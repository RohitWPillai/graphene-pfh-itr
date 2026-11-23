#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cmath>

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
    double ix, iy, iz;
};

// Structure to hold data for one frame
struct FrameData {
    int timeStep;
    std::map<int, std::vector<std::array<double, 3>>> atomicVelocities; // molID -> vector of {vx, vy, vz} for ordered atoms
};

// Get the mass for an atom based on its type
double getMass(int type) {
    if (type == 1) return 12.011;   // carbonTip
    if (type == 2) return 12.011;   // carbonBackbone
    if (type == 3) return 18.998;   // fluorine
    return 1.0; // default mass
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_dump_file> <output_timeSeries_file> <zthreshold> <blockSize>\n\n"
                  << "This program:\n"
                  << "  1) Reads a LAMMPS dump file and discards the first frame.\n"
                  << "  2) Processes frames in blocks of <blockSize>.\n"
                  << "  3) Filters PFH molecules based on COM z <= zthreshold.\n"
                  << "  4) Outputs atomic velocities for VDOS:\n"
                  << "       timeStep molID [vx vy vz for each of 20 atoms]\n";
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    double zThreshold = std::stod(argv[3]);
    int blockSize = std::stoi(argv[4]);

    std::ifstream fin(inputFile);
    if (!fin.is_open()) {
        std::cerr << "Error: Could not open input file: " << inputFile << "\n";
        return 1;
    }

    // Skip the first frame
    {
        std::string line;
        while (std::getline(fin, line)) {
            if (line.find("ITEM: TIMESTEP") != std::string::npos) {
                std::getline(fin, line); // timestep
                std::getline(fin, line); // NUMBER OF ATOMS header
                int nAtoms;
                std::getline(fin, line); std::istringstream(line) >> nAtoms;
                std::getline(fin, line); // box bounds header
                for (int i = 0; i < 3; i++) std::getline(fin, line);
                std::getline(fin, line); // ATOMS header
                for (int i = 0; i < nAtoms; i++) std::getline(fin, line);
                break;
            }
        }
    }

    const int EXPECTED_FRAMES = 250000;
    int totalFrameCount = 0;
    int totalBlockCount = 0;

    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        std::cerr << "Error: Could not create output file: " << outputFile << "\n";
        return 1;
    }

    while (totalFrameCount < EXPECTED_FRAMES) {
        std::vector<FrameData> blockFrames;
        std::set<int> blockGlobalIntersection;
        bool blockFirstFrame = true;
        int blockFrameCount = 0;
        std::map<int, std::vector<int>> orderedAtomIDs; // molID -> ordered atom IDs

        while (blockFrameCount < blockSize && totalFrameCount < EXPECTED_FRAMES) {
            std::string line;
            while (std::getline(fin, line)) {
                if (line.find("ITEM: TIMESTEP") != std::string::npos) break;
            }
            if (fin.eof()) break;

            int timeStep;
            std::getline(fin, line); std::istringstream(line) >> timeStep;

            std::getline(fin, line); // NUMBER OF ATOMS
            int nAtoms;
            std::getline(fin, line); std::istringstream(line) >> nAtoms;

            std::getline(fin, line); // box bounds header
            for (int i = 0; i < 3; i++) std::getline(fin, line);

            std::getline(fin, line); // ATOMS header

            std::map<int, std::vector<AtomData>> atomsByMol;
            for (int i = 0; i < nAtoms; i++) {
                if (!std::getline(fin, line)) break;
                if (line.empty() || line[0] == '#') { i--; continue; }
                std::istringstream iss(line);
                AtomData atom;
                if (iss >> atom.id >> atom.type >> atom.mol 
                        >> atom.x >> atom.y >> atom.z 
                        >> atom.vx >> atom.vy >> atom.vz) {
                    atomsByMol[atom.mol].push_back(atom);
                }
            }

            // In first frame, determine ordered atom IDs
            if (blockFirstFrame) {
                for (const auto &pair : atomsByMol) {
                    int molID = pair.first;
                    const auto &atomList = pair.second;
                    std::vector<int> carbonIDs;
                    std::vector<int> fluorineIDs;
                    for (const auto &atom : atomList) {
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
                blockFirstFrame = false;
            }

            FrameData frameData;
            frameData.timeStep = timeStep;

            // Store atomic velocities
            for (const auto &pair : orderedAtomIDs) {
                int molID = pair.first;
                const auto &orderedIDs = pair.second;
                if (atomsByMol.count(molID) == 0) continue;
                const auto &atomList = atomsByMol[molID];
                std::map<int, AtomData> atomMap;
                for (const auto &atom : atomList) {
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

            // Filter based on COM position
            std::set<int> currentQualifyingIDs;
            for (const auto &pair : atomsByMol) {
                int molID = pair.first;
                const auto &atomList = pair.second;
                double totalMass = 0.0;
                double com_z = 0.0;
                for (const auto &atom : atomList) {
                    double m = getMass(atom.type);
                    totalMass += m;
                    com_z += m * atom.z;
                }
                if (totalMass > 0.0) {
                    com_z /= totalMass;
                    if (com_z <= zThreshold) {
                        currentQualifyingIDs.insert(molID);
                    }
                }
            }

            if (blockFrameCount == 0) {
                blockGlobalIntersection = currentQualifyingIDs;
            } else {
                blockGlobalIntersection = intersectSets(blockGlobalIntersection, currentQualifyingIDs);
            }

            blockFrames.push_back(std::move(frameData));
            blockFrameCount++;
            totalFrameCount++;
        }

        if (blockFrameCount == 0) break;

        totalBlockCount++;
        int numFramesInBlock = blockFrames.size();
        int numMolsInBlock = blockGlobalIntersection.size();
        std::cout << "Block " << totalBlockCount << ": Processed " << numFramesInBlock
                  << " frames; found " << numMolsInBlock << " qualifying molecules.\n";

        long long totalDataLines = static_cast<long long>(numFramesInBlock) * numMolsInBlock;
        fout << "# Block " << totalBlockCount << "\n";
        fout << "# NumMolecules " << numMolsInBlock << "\n";
        fout << "# NumLines " << totalDataLines << "\n";
        fout << "# Columns: timeStep molID";
        for (int i = 1; i <= 6; i++) fout << " vx_c" << i << " vy_c" << i << " vz_c" << i;
        for (int i = 1; i <= 14; i++) fout << " vx_f" << i << " vy_f" << i << " vz_f" << i;
        fout << "\n";

        for (const auto &frame : blockFrames) {
            int t = frame.timeStep;
            for (int molID : blockGlobalIntersection) {
                fout << t << " " << molID;
                auto it = frame.atomicVelocities.find(molID);
                if (it != frame.atomicVelocities.end()) {
                    const auto &velocities = it->second;
                    for (const auto &vel : velocities) {
                        fout << " " << vel[0] << " " << vel[1] << " " << vel[2];
                    }
                } else {
                    for (int i = 0; i < 20; i++) fout << " 0 0 0";
                }
                fout << "\n";
            }
        }
        fout << "\n";
    }

    fin.close();
    fout.close();

    std::cout << "Processed " << totalBlockCount << " blocks (" << totalFrameCount << " frames).\n";
    std::cout << "Output written to " << outputFile << "\n";

    return 0;
}
