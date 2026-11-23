#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <array>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include <fftw3.h>
// Helper: Next power of two
static unsigned int nextPowerOfTwo(unsigned int x)
{
    if (x == 0) return 1;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
}
// FFT routine for 1D time series
void calcFFT_1D(unsigned int fft_len, std::vector<double>& fft_conj, const std::vector<double>& time_series)
{
    unsigned int n = time_series.size();
    fftw_complex* fft_inp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_len);
    fftw_complex* fft_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_len);
    fftw_plan fft_plan = fftw_plan_dft_1d(fft_len, fft_inp, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
    for (unsigned int i = 0; i < n; ++i) {
        fft_inp[i][0] = time_series[i];
        fft_inp[i][1] = 0.0;
    }
    for (unsigned int i = n; i < fft_len; ++i) {
        fft_inp[i][0] = 0.0;
        fft_inp[i][1] = 0.0;
    }
    fftw_execute(fft_plan);
    for (unsigned int i = 0; i < fft_len; ++i) {
        double re = fft_out[i][0];
        double im = fft_out[i][1];
        fft_conj[i] = re * re + im * im;
    }
    fftw_destroy_plan(fft_plan);
    fftw_free(fft_inp);
    fftw_free(fft_out);
}
// Helper: Intersect two sets
template <typename T>
std::set<T> intersectSets(const std::set<T>& s1, const std::set<T>& s2) {
    std::set<T> result;
    std::set_intersection(
        s1.begin(), s1.end(),
        s2.begin(), s2.end(),
        std::inserter(result, result.begin())
    );
    return result;
}
// Structure to hold per-atom data from the dump file
struct AtomData {
    int id;
    int mol;
    int type;
    double x, y, z;
    double vx, vy, vz;
    int ix, iy, iz;
};
// Structure to hold data for one frame
struct FrameData {
    int timeStep;
    std::map<int, std::vector<std::array<double, 3>>> atomicVelocities; // molID -> vector of {vx, vy, vz} for ordered atoms
};
// Get the mass for an atom based on its type (in g/mol)
double getMass(int type) {
    if (type == 1 || type == 2) return 12.011; // carbon
    if (type == 3) return 18.998; // fluorine
    return 1.0; // default
}
// Function to count total number of frames in the dump file
int countTotalFrames(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open()) return 0;
    std::string line;
    int count = 0;
    while (std::getline(fin, line)) {
        if (line.find("ITEM: TIMESTEP") == 0) {
            ++count;
        }
    }
    fin.close();
    return count;
}
// Function to read the next frame from the dump file
bool readNextFrame(std::ifstream& fin, int& timeStep, int& nAtoms, std::map<int, std::vector<AtomData>>& atomsByMol) {
    atomsByMol.clear();
    std::string line;
    // Find ITEM: TIMESTEP
    while (std::getline(fin, line)) {
        if (line.find("ITEM: TIMESTEP") == 0) {
            break;
        }
    }
    if (!fin.good()) return false;
    // Read timestep value
    if (!std::getline(fin, line)) return false;
    std::istringstream iss_ts(line);
    iss_ts >> timeStep;
    // Find ITEM: NUMBER OF ATOMS
    while (std::getline(fin, line)) {
        if (line.find("ITEM: NUMBER OF ATOMS") == 0) {
            break;
        }
    }
    // Read natoms value
    if (!std::getline(fin, line)) return false;
    std::istringstream iss_na(line);
    iss_na >> nAtoms;
    // Find ITEM: BOX BOUNDS
    while (std::getline(fin, line)) {
        if (line.find("ITEM: BOX BOUNDS") == 0) {
            break;
        }
    }
    // Read 3 box lines
    for (int i = 0; i < 3; ++i) {
        if (!std::getline(fin, line)) return false;
    }
    // Find ITEM: ATOMS
    while (std::getline(fin, line)) {
        if (line.find("ITEM: ATOMS") == 0) {
            break;
        }
    }
    // Read atom lines
    for (int i = 0; i < nAtoms; ++i) {
        if (!std::getline(fin, line)) {
            std::cerr << "Error: Unexpected end of file or read failure while reading atom line " << (i + 1) << " of " << nAtoms << " at timestep " << timeStep << "\n";
            return false;
        }
        if (line.empty() || line[0] == '#') {
            --i;
            continue;
        }
        std::istringstream iss(line);
        AtomData atom;
        if (!(iss >> atom.id >> atom.mol >> atom.type >> atom.x >> atom.y >> atom.z >> atom.ix >> atom.iy >> atom.iz >> atom.vx >> atom.vy >> atom.vz)) {
            std::cerr << "Warning: incomplete atom data in line: " << line << "\n";
            --i;
            continue;
        }
        atomsByMol[atom.mol].push_back(atom);
    }
    return true;
}
int main(int argc, char* argv[])
{
    if (argc < 8) {
        std::cerr << "Usage: " << argv[0] << " <input_dump_file> <output_file> <z_lim_min> <z_lim_max> <dt> <segment_size> <segment_number>\n\n"
                  << "This program processes a single LAMMPS dump file to compute the vibrational density of states (VDOS)\n"
                  << "for interfacial PFH molecules that stay within [z_lim_min, z_lim_max] throughout the specified segment.\n"
                  << "Skips the first frame (frame 0), then processes the segment_number-th segment of segment_size frames\n"
                  << "(frame numbering starts at 1 after skip). If the segment is incompatible (e.g., starts beyond available frames),\n"
                  << "terminates with a warning. Processes partial last segments if applicable.\n"
                  << "Velocities in Ang/ps, masses in g/mol.\n"
                  << "Output: Frequency, PowerX, PowerY, PowerZ, Sum (X+Y+Z).\n";
        return 1;
    }
    std::string inputFile = argv[1];
    std::string outFile = argv[2];
    double zLimMin = std::stod(argv[3]);
    double zLimMax = std::stod(argv[4]);
    double dt = std::stod(argv[5]);
    int segment_size = std::stoi(argv[6]);
    int segment_number = std::stoi(argv[7]);
    if (dt <= 0) {
        std::cerr << "Error: dt must be positive.\n";
        return 1;
    }
    if (segment_size <= 0) {
        std::cerr << "Error: segment_size must be positive.\n";
        return 1;
    }
    if (segment_number < 1) {
        std::cerr << "Error: segment_number must be at least 1.\n";
        return 1;
    }
    if (zLimMin > zLimMax) {
        std::cerr << "Error: z_lim_min should be less than or equal to z_lim_max.\n";
        return 1;
    }
    // Count total frames
    int total_frames = countTotalFrames(inputFile);
    if (total_frames <= 1) {
        std::cerr << "Error: Insufficient frames in dump file (need at least 2: one to skip, one to process).\n";
        return 1;
    }
    int total_usable_frames = total_frames - 1;
    // Calculate start and end
    int start_frame = (segment_number - 1) * segment_size + 1;
    int end_frame = start_frame + segment_size - 1;
    if (start_frame > total_usable_frames) {
        std::cerr << "Warning: Segment " << segment_number << " starts beyond available frames (" << total_usable_frames << "). Terminating.\n";
        return 0;
    }
    if (end_frame > total_usable_frames) {
        end_frame = total_usable_frames;
    }
    std::ifstream fin(inputFile);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open input file " << inputFile << "\n";
        return 1;
    }
    // Skip the first frame (frame 0)
    int dummy_ts, dummy_na;
    std::map<int, std::vector<AtomData>> dummy_atoms;
    if (!readNextFrame(fin, dummy_ts, dummy_na, dummy_atoms)) {
        std::cerr << "Error: Failed to skip first frame.\n";
        return 1;
    }
    // Skip additional frames up to start_frame - 1
    for (int cf = 1; cf < start_frame; ++cf) {
        if (!readNextFrame(fin, dummy_ts, dummy_na, dummy_atoms)) {
            std::cerr << "Error: File ended before reaching start_frame " << start_frame << ".\n";
            return 1;
        }
    }
    // Now process the single segment from start_frame to end_frame
    std::cout << "Processing segment " << segment_number << " from frame " << start_frame << " to " << end_frame << "\n";
    std::vector<FrameData> segmentFrames;
    std::set<int> segmentIntersection;
    bool firstInSegment = true;
    std::map<int, std::vector<int>> orderedAtomIDs;
    int framesInSegment = 0;
    int timeStep, nAtoms;
    std::map<int, std::vector<AtomData>> atomsByMol;
    int current_frame = start_frame;
    while (current_frame <= end_frame) {
        if (!readNextFrame(fin, timeStep, nAtoms, atomsByMol)) break;
        std::cout << "nAtoms equals = " << nAtoms << std::endl;
        // Determine ordered atom IDs in the first frame of the segment
        if (firstInSegment) {
            for (const auto& pair : atomsByMol) {
                int molID = pair.first;
                const auto& atomList = pair.second;
                std::vector<int> carbonIDs;
                std::vector<int> fluorineIDs;
                for (const auto& atom : atomList) {
                    if (atom.type == 1 || atom.type == 2) {
                        carbonIDs.push_back(atom.id);
                    } else if (atom.type == 3) {
                        fluorineIDs.push_back(atom.id);
                    }
                }
                std::sort(carbonIDs.begin(), carbonIDs.end());
                std::sort(fluorineIDs.begin(), fluorineIDs.end());
                std::vector<int> orderedIDs = carbonIDs;
                orderedIDs.insert(orderedIDs.end(), fluorineIDs.begin(), fluorineIDs.end());
                orderedAtomIDs[molID] = orderedIDs;
            }
            firstInSegment = false;
        }
        // Compute qualifying molecules based on COM z
        std::set<int> currentQualifying;
        for (const auto& pair : atomsByMol) {
            int molID = pair.first;
            const auto& atomList = pair.second;
            double totalMass = 0.0;
            double comZ = 0.0;
            for (const auto& atom : atomList) {
                double m = getMass(atom.type);
                totalMass += m;
                comZ += m * atom.z;
            }
            if (totalMass > 0.0) {
                comZ /= totalMass;
                if (comZ >= zLimMin && comZ <= zLimMax) {
                    currentQualifying.insert(molID);
                }
            }
        }
        // Update intersection
        if (framesInSegment == 0) {
            segmentIntersection = currentQualifying;
        } else {
            segmentIntersection = intersectSets(segmentIntersection, currentQualifying);
        }
        // Store frame data
        FrameData frameData;
        frameData.timeStep = timeStep;
        for (const auto& pair : orderedAtomIDs) {
            int molID = pair.first;
            const auto& orderedIDs = pair.second;
            if (atomsByMol.count(molID) == 0) continue;
            const auto& atomList = atomsByMol[molID];
            std::map<int, AtomData> atomMap;
            for (const auto& atom : atomList) {
                atomMap[atom.id] = atom;
            }
            std::vector<std::array<double, 3>> velocities;
            for (int id : orderedIDs) {
                auto it = atomMap.find(id);
                if (it != atomMap.end()) {
                    velocities.push_back({it->second.vx, it->second.vy, it->second.vz});
                } else {
                    velocities.push_back({0.0, 0.0, 0.0});
                }
            }
            frameData.atomicVelocities[molID] = velocities;
        }
        segmentFrames.push_back(std::move(frameData));
        framesInSegment++;
        current_frame++;
    }
    fin.close();
    if (framesInSegment == 0) {
        std::cerr << "Error: No frames found in the specified segment.\n";
        return 1;
    }
    std::cout << "Segment contains " << framesInSegment
              << " frames and " << segmentIntersection.size() << " qualifying molecules.\n";
    if (segmentIntersection.empty()) {
        std::cerr << "Skipping segment due to no qualifying molecules.\n";
        return 0; // Or output empty file if preferred
    }
    unsigned int N = framesInSegment;
    unsigned int fft_len = nextPowerOfTwo(2 * N);
    std::vector<std::vector<double>> powerSum(fft_len, std::vector<double>(3, 0.0));
    int validMolCount = 0;
    for (int molID : segmentIntersection) {
        std::vector<std::vector<double>> samples(N, std::vector<double>(60, 0.0)); // 20 atoms * 3 dims
        bool complete = true;
        for (unsigned int i = 0; i < N; ++i) {
            const auto& frame = segmentFrames[i];
            auto it = frame.atomicVelocities.find(molID);
            if (it == frame.atomicVelocities.end()) {
                complete = false;
                break;
            }
            const auto& vels = it->second;
            if (vels.size() != 20) {
                complete = false;
                break;
            }
            for (int a = 0; a < 20; ++a) {
                samples[i][a * 3 + 0] = vels[a][0];
                samples[i][a * 3 + 1] = vels[a][1];
                samples[i][a * 3 + 2] = vels[a][2];
            }
        }
        if (!complete) continue;
        validMolCount++;
        std::vector<std::vector<double>> total_power(fft_len, std::vector<double>(3, 0.0));
        for (int a = 0; a < 20; ++a) {
            double mass = (a < 6) ? 12.011 : 18.998;
            for (int dim = 0; dim < 3; ++dim) {
                std::vector<double> time_series(N);
                for (unsigned int i = 0; i < N; ++i) {
                    time_series[i] = samples[i][a * 3 + dim];
                }
                std::vector<double> fft_conj(fft_len, 0.0);
                calcFFT_1D(fft_len, fft_conj, time_series);
                for (unsigned int bin = 0; bin < fft_len; ++bin) {
                    total_power[bin][dim] += mass * fft_conj[bin];
                }
            }
        }
        for (unsigned int i = 0; i < fft_len; ++i) {
            powerSum[i][0] += total_power[i][0];
            powerSum[i][1] += total_power[i][1];
            powerSum[i][2] += total_power[i][2];
        }
    }
    if (validMolCount == 0) {
        std::cerr << "Segment has 0 valid molecules. Skipping.\n";
        return 0; // Or output empty file if preferred
    }
    std::cout << "Segment: Valid molecules = " << validMolCount << "\n";
    std::ofstream fout(outFile);
    if (!fout.is_open()) {
        std::cerr << "Error: could not open output file " << outFile << "\n";
        return 1;
    }
    unsigned int nFreq = fft_len / 2 + 1;
    fout << "# zero padded to fft_len = " << fft_len << " (padded by a factor of 2)\n";
    fout << "# Columns: Frequency\tPowerX\tPowerY\tPowerZ\tSumPower\n";
    double df = 1.0 / (fft_len * dt);
    for (unsigned int i = 0; i < nFreq; ++i) {
        double freq = i * df;
        double avgX = powerSum[i][0] / validMolCount;
        double avgY = powerSum[i][1] / validMolCount;
        double avgZ = powerSum[i][2] / validMolCount;
        double sumPower = avgX + avgY + avgZ;
        fout << freq << "\t" << avgX << "\t" << avgY << "\t" << avgZ << "\t" << sumPower << "\n";
    }
    fout.close();
    std::cout << "Processed segment (" << framesInSegment << " frames).\n";
    std::cout << "Power spectrum written to " << outFile << "\n";
    return 0;
}
