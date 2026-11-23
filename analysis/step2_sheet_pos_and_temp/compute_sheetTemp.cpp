#include <iostream>
#include <iomanip>
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

// Structure to hold per-atom data from the dump file
struct AtomData {
    int id;
    int mol;
    int type;
    double x, y, z;
    double vx, vy, vz;
    double ix, iy, iz;
};

// Get the mass for an atom based on its type (carbon for graphene)
double getMass(int type) {
    return 12.011; // g/mol
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_dump_file> <output_file> <zLo> <zHi> <blockSize>\n\n"
                  << "This program:\n"
                  << "  1) Reads a LAMMPS dump file and discards the first frame.\n"
                  << "  2) Processes frames in blocks of <blockSize>.\n"
                  << "  3) Computes temperature and COM for each graphene sheet.\n"
                  << "  4) Outputs sheet type, COM, and temperature per block.\n";
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    double zLo = std::stod(argv[3]);
    double zHi = std::stod(argv[4]);
    int blockSize = std::stoi(argv[5]);

    /*
        std::cout << " Input file is : "    << inputFile  << std::endl;
        std::cout << " Output file is : "   << outputFile << std::endl;
        std::cout << " zLo file is : "      << zLo       << std::endl;
        std::cout << " zHi file is : "      << zHi       << std::endl;
        std::cout << " blockSize is : "     << blockSize << std::endl;
    */

    double nAv = 6.02214076e23;         // Avogadro's number, 1/mol
    double kB = 1.380649e-23;           // Boltzmann constant, J/K
    double gmol_to_kg = 1e-3 / nAv;     // g/mol to kg/atom
    double ang2_ps2_to_m2_s2 = 1e4;     // (Å/ps)^2 to (m/s)^2
    double to_joules = ang2_ps2_to_m2_s2 * gmol_to_kg; // Energy to joules
    double tempConvertSI = to_joules / (3.0 * kB);     // Conversion to K

    double boxLen[3];

    // Maps to accumulate data per sheet (keyed by atom type)
    std::map<int, double> ke;        // Kinetic energy per sheet
    std::map<int, int> atom_count;   // Atom count per sheet
    std::map<int, double> sum_z;     // Sum of z-positions per sheet

    // Open input file
    std::ifstream fin(inputFile);
    if (!fin.is_open()) {
        std::cerr << "Error: Could not open input file: " << inputFile << "\n";
        return 1;
    }

    // Open output file
    std::ofstream fout(outputFile);
    if (!fout.is_open()) {
        std::cerr << "Error: Could not open output file: " << outputFile << "\n";
        fin.close();
        return 1;
    }

    // Skip the first frame
    {
        std::string line;
        while (std::getline(fin, line)) {
            if (line.find("ITEM: TIMESTEP") != std::string::npos) {
                std::getline(fin, line); // timestep
                std::getline(fin, line); // NUMBER OF ATOMS
                int nAtoms;
                std::getline(fin, line); std::istringstream(line) >> nAtoms;
                std::getline(fin, line); // box bounds header
                for (int i = 0; i < 3; i++) {
                    std::getline(fin, line);
                    std::istringstream iss(line);
                    double Lo, Hi;
                    iss >> Lo >> Hi;
                    boxLen[i] = Hi - Lo;
                }
                std::getline(fin, line); // ATOMS header
                for (int i = 0; i < nAtoms; i++) std::getline(fin, line);
                break;
            }
        }
    }

    std::cout << "Temperature conversion factor: " << tempConvertSI << "\n";
    std::cout << "Box lengths: x=" << boxLen[0] << ", y=" << boxLen[1] << ", z=" << boxLen[2] << "\n";

    const int EXPECTED_FRAMES = 250000;
    int totalFrameCount = 0;
    int totalBlockCount = 0;

    while (totalFrameCount < EXPECTED_FRAMES) {
        int blockFrameCount = 0;

        // Clear maps for each block
        ke.clear();
        atom_count.clear();
        sum_z.clear();

        while (blockFrameCount < blockSize && totalFrameCount < EXPECTED_FRAMES) {
            std::string line;
            while (std::getline(fin, line)) {if (line.find("ITEM: TIMESTEP") != std::string::npos) break;}

            if (fin.eof()) break;

            std::getline(fin, line); // timestep
            int timeStep;
            std::istringstream(line) >> timeStep;

            std::getline(fin, line); // NUMBER OF ATOMS
            int nAtoms;
            std::getline(fin, line); std::istringstream(line) >> nAtoms;

            std::getline(fin, line); // box bounds header
            for (int i = 0; i < 3; i++) std::getline(fin, line);

            std::getline(fin, line); // ATOMS header

            for (int i = 0; i < nAtoms; i++) {
                if (!std::getline(fin, line)) break;
                if (line.empty() || line[0] == '#') { i--; continue; }

                std::istringstream iss(line);
                AtomData atom;
                if ( iss >> atom.id >> atom.type
                         >> atom.x >> atom.y >> atom.z 
                         >> atom.vx >> atom.vy >> atom.vz ) {
                    if ( atom.z >= zLo && atom.z <= zHi ) {
                        int sheet = atom.type; // Sheet identified by atom type
                        double mass = getMass(atom.type);
                        double v2 = atom.vx * atom.vx + atom.vy * atom.vy + atom.vz * atom.vz;
                        double Two_ke = mass * v2;

                        ke[sheet] += Two_ke;
                        atom_count[sheet] += 1;
                        sum_z[sheet] += atom.z;
                    }
                }
            }

            blockFrameCount++;
            totalFrameCount++;

            if (totalFrameCount % 100 == 0) {
                std::cout << totalFrameCount << " frames completed" << std::endl;
            }
        }

        if (blockFrameCount == 0) break;

        totalBlockCount++;
        std::cout << "Block " << totalBlockCount << ": Processed " << blockFrameCount << " frames\n";

        fout << "# Graphene sheet temperature profile\n"
             << "# Block: " << totalBlockCount << "\n"
             << "# z range: " << zLo << " to " << zHi << " Å\n"
             << "# Columns: sheet_type COM_z [Å] temperature [K]\n";

        for (const auto& [sheet, ke_sum] : ke) {
            if (atom_count[sheet] > 0) {
                double keAv = ke_sum / static_cast<double>(atom_count[sheet]);
                double temp = keAv * tempConvertSI;
                double com_z = sum_z[sheet] / atom_count[sheet];

                fout << std::setprecision(4)
                     << std::setw(15) << sheet
                     << std::setw(15) << com_z
                     << std::setw(15) << temp << "\n";
            }
        }
    }

    fin.close();
    fout.close();

    std::cout << "Processed " << totalBlockCount << " blocks (" << totalFrameCount << " frames).\n";
    std::cout << "Output written to " << outputFile << "\n";

    return 0;
}
