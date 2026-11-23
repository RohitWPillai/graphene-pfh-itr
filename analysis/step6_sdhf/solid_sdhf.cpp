#include "utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <dirent.h>
#include <algorithm>
#include <fftw3.h>

struct TempAtom {
    int id;
    double vx, vy, vz, fx, fy, fz;
};

// Process atom lines into atom data for one frame
std::vector<AtomData> processFrameToAtoms(const std::vector<std::string>& atomLines, const std::map<int, int>& idToIdx) {
    std::vector<AtomData> frameData(idToIdx.size());
    for (const auto& line : atomLines) {
        std::istringstream iss(line);
        TempAtom atom;
        int mol, type, ix, iy, iz;
        double x, y, z;
        iss >> atom.id >> type >> x >> y >> z >> atom.vx >> atom.vy >> atom.vz >> atom.fx >> atom.fy >> atom.fz;
        auto it = idToIdx.find( atom.id );
        if (it != idToIdx.end()) {
            frameData[it->second] = {atom.vx, atom.vy, atom.vz, atom.fx, atom.fy, atom.fz};
        }
    }
    return frameData;
}

// Load exactly L sampled frames starting at global sampled index `start`
void loadWindow(int start, int L, const std::vector<std::string>& files, int framesPerFileRaw, int sampling, const std::map<int, int>& idToIdx, std::vector<std::vector<AtomData>>& atomTS) {
    std::cout << "Loading window starting at sampled index " << start << ", length " << L << std::endl;
    int framesPerSample = framesPerFileRaw / sampling;
    int raw_start = 1 + start * sampling;  // Skip the first raw frame (index 0), start from 1 + start * sampling
    int fileIdx = raw_start / framesPerFileRaw;
    int offsetInFileRaw = raw_start % framesPerFileRaw;
    int got = 0;
    int natoms = idToIdx.size();

    while (got < L) {
        std::cout << "Opening file: " << files[fileIdx] << ", raw offset: " << offsetInFileRaw << std::endl;
        std::ifstream fin(files[fileIdx]);
        if (!fin.is_open()) {
            std::cerr << "Error: Could not open file: " << files[fileIdx] << std::endl;
            exit(1);
        }
        // Skip to the offsetInFileRaw-th raw frame
        std::string line;
        for (int f = 0; f < offsetInFileRaw; ++f) {
            while (std::getline(fin, line) && line.find("ITEM: ATOMS") != 0);
            for (int i = 0; i < natoms; ++i) std::getline(fin, line);
        }

        // Determine how many sampled frames to get from this file
        int avail_raw = framesPerFileRaw - offsetInFileRaw;
        int avail_sampled = (avail_raw - 1) / sampling + 1;
        int toGet = std::min(L - got, avail_sampled);

        if (toGet > 0) {
            int first_local_raw = offsetInFileRaw;
            int last_local_raw = offsetInFileRaw + (toGet - 1) * sampling;
            std::cout << "Reading " << toGet << " sampled frames from this file, local raw frames "
                      << first_local_raw << " to " << last_local_raw << " every " << sampling << std::endl;
        }

        for (int ss = 0; ss < toGet; ++ss) {
            while (std::getline(fin, line) && line.find("ITEM: ATOMS") != 0);
            std::vector<std::string> atomLines(natoms);
            for (int i = 0; i < natoms; ++i) std::getline(fin, atomLines[i]);
            std::vector<AtomData> frameData = processFrameToAtoms(atomLines, idToIdx);
            for (int a = 0; a < natoms; ++a) {
                atomTS[a][got] = frameData[a];
            }
            ++got;
            // Skip sampling - 1 raw frames
            for (int rr = 0; rr < sampling - 1; ++rr) {
                while (std::getline(fin, line) && line.find("ITEM: ATOMS") != 0);
                for (int i = 0; i < natoms; ++i) std::getline(fin, line);
            }
        }
        offsetInFileRaw = 0;
        ++fileIdx;
    }
    std::cout << "Loaded " << got << " frames for the window." << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 8) {
        std::cerr << "Usage: " << argv[0] << " <folder> <suffix> <outprefix> <win_len> <over_pct> <sampling> <segment_id>\n";
        std::cerr << "Note: sampling is required, segment_id is 0-based.\n";
        return 1;
    }

    std::string folder = argv[1];
    std::string suffix = argv[2];
    std::string outPref = argv[3];
    int L = std::atoi(argv[4]);
    double overPct = std::atof(argv[5]);
    int sampling = std::atoi(argv[6]);
    int s = std::atoi(argv[7]);

    int overlap = static_cast<int>(std::round(L * overPct / 100.0));
    int step = L - overlap;

    // Gather and sort dump files
    std::vector<std::pair<int, std::string>> ent;
    DIR* dp = opendir(folder.c_str());
    if (!dp) {
        std::cerr << "Error: Could not open directory: " << folder << std::endl;
        return 1;
    }
    struct dirent* de;
    while ((de = readdir(dp)) != nullptr) {
        if (de->d_type != DT_REG) continue;
        std::string nm(de->d_name);
        if (nm.rfind("dump", 0) != 0) continue;
        size_t dot = nm.find('.');
        if (dot == std::string::npos) continue;
        if (nm.substr(dot + 1) != suffix) continue;
        int idx = std::stoi(nm.substr(4, dot - 4));
        ent.emplace_back(idx, folder + "/" + nm);
    }
    closedir(dp);
    std::sort(ent.begin(), ent.end());
    std::vector<std::string> files;
    for (auto& e : ent) files.push_back(e.second);

    if (files.empty()) {
        std::cerr << "Error: No dump files found.\n";
        return 1;
    }

    std::cout << "Found " << files.size() << " dump files." << std::endl;

    // Read first file to get natoms and framesPerFileRaw
    std::ifstream f1(files[0]);
    std::string line;
    int natoms = 0;
    int framesPerFileRaw = 0;
    while (std::getline(f1, line)) {
        if (line == "ITEM: NUMBER OF ATOMS") {
            std::getline(f1, line);
            natoms = std::stoi(line);
        } else if (line == "ITEM: TIMESTEP") {
            ++framesPerFileRaw;
        }
    }
    f1.close();

    std::cout << "natoms = " << natoms << ", framesPerFileRaw = " << framesPerFileRaw << std::endl;

    // Read the first frame of the first file to get atom IDs and create mapping
    std::map<int, int> idToIdx;
    {
        std::ifstream fin(files[0]);
        while (std::getline(fin, line) && line.find("ITEM: ATOMS") != 0);
        std::set<int> atomIDs;
        for (int i = 0; i < natoms; ++i) {
            std::getline(fin, line);
            std::istringstream iss(line);
            int id, type;
            iss >> id >> type;
            atomIDs.insert(id);
        }
        if (atomIDs.size() != static_cast<size_t>(natoms)) {
            std::cerr << "Error: Number of unique atom IDs does not match natoms." << std::endl;
            return 1;
        }
        int idx = 0;
        int minID = *atomIDs.begin();
        int maxID = *atomIDs.rbegin();
        for (int id : atomIDs) {
            idToIdx[id] = idx++;
        }
        std::cout << "Atom IDs range from " << minID << " to " << maxID << std::endl;
    }

    // Compute maximum number of segments
    long long total_raw = (long long)files.size() * framesPerFileRaw;
    long long max_raw_index = total_raw - 1;
    long long max_k = (max_raw_index - 1) / sampling;
    long long total_sampled = max_k + 1;
    long long max_start = total_sampled - L;
    if (max_start < 0) {
        std::cerr << "Error: Window length L is too large for the available data." << std::endl;
        return 1;
    }
    long long nSeg = (max_start / step) + 1;
    std::cout << "Maximum number of segments: " << nSeg << std::endl;

    if (s >= nSeg || s < 0) {
        std::cerr << "Error: Segment " << s << " is out of range. Must be between 0 and " << (nSeg - 1) << " inclusive." << std::endl;
        return 1;
    }

    int start = s * step;

    std::cout << "Processing segment " << s << ", starting sampled index = " << start << std::endl;

    // Build Welch window
    std::vector<double> W(L);
    double U = 0.0;
    double half = (L - 1) / 2.0;
    for (int i = 0; i < L; ++i) {
        double t = (i - half) / half;
        W[i] = 1 - t * t;
        U += W[i] * W[i];
    }
    U /= L;

    std::cout << "Value of U is : " << U << std::endl;

    // Allocate FFTW arrays
    fftw_complex* fft_in = alloc_fftw_complex(L);
    fftw_complex* fft_out = alloc_fftw_complex(L);
    fftw_complex* fft_res = alloc_fftw_complex(L);
    fftw_plan plan = fftw_plan_dft_1d(L, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE);

    // Prepare segment accumulators
    int outSize = L / 2 + 1;
    std::vector<double> seg_x(outSize, 0.0), seg_y(outSize, 0.0), seg_z(outSize, 0.0);

    // Load the window
    std::vector<std::vector<AtomData>> atomTS(natoms, std::vector<AtomData>(L));
    loadWindow(start, L, files, framesPerFileRaw, sampling, idToIdx, atomTS);

    // Process each atom
    for (int a = 0; a < natoms; ++a) {
        computeAtomContribution(atomTS[a], L, fft_in, fft_out, fft_res, plan, seg_x, seg_y, seg_z, outSize, W);
    }

    // Output the segment
    double timestep = 0.001 * 1e-12;  // ps to s
    double dump_freq = 5.0;           // dump every 5 timesteps
    double dt = timestep * dump_freq * sampling;
    double conv = (2.0 * dt / L) / U;

    std::ofstream souts(outPref + "_seg" + std::to_string(s) + ".txt");
    souts << "# Segment: " << s << "\n";
    souts << "# Frequency[Hz]    Sx    Sy    Sz\n"
          << std::scientific << std::setw(20) << std::setprecision(7);
    for (int k = 0; k < outSize; ++k) {
        double freq = k / (dt * L);
        souts << std::setw(20) << freq
              << std::setw(20) << seg_x[k] * conv
              << std::setw(20) << seg_y[k] * conv
              << std::setw(20) << seg_z[k] * conv << "\n";
    }
    souts.close();

    // Cleanup
    free_fftw_complex(fft_in);
    free_fftw_complex(fft_out);
    free_fftw_complex(fft_res);
    fftw_destroy_plan(plan);
    fftw_cleanup();

    return 0;
}
