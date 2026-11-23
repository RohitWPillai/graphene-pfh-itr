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

// Get the mass for an atom based on its type
double getMass( int type ) {
    if ( type == 1 ) return 12.017;   // carbonTip
    if ( type == 2 ) return 12.017;   // carbonBackbone
    if ( type == 3 ) return 18.998;   // fluorine
    return 1.0; // default mass
}

int main( int argc, char* argv[] ) {
    if ( argc < 5 ) {
        std::cerr << "Usage: " << argv[0]
            << " <input_dump_file> <output_density_file> <zHot> <zCold> <nBins> <blockSize>\n\n"
            << "This program:\n"
            << "  1) Reads a LAMMPS dump file and discards the first frame.\n"
            << "  2) Processes frames in blocks of <blockSize>.\n"
            << "  3) Computes density in SI units for PFH between zHot and zCold.\n"
            << "  4) Outputs densities per bin for each block.\n";
        return 1;
    }

    std::string inputFile       = argv[1];
    std::string outputFile      = argv[2];
    double      zHot            = std::stod( argv[3] );
    double      zCold           = std::stod( argv[4] );
    double      nBins           = std::stoi( argv[5] );
    int         blockSize       = std::stoi( argv[6] );
    std::string domainSide      = argv[7];

    double      wallSep             = std::abs(std::abs(zCold) - std::abs(zHot));
    double      binw                = wallSep / nBins;          // in Ang.
    double      Ang3_to_m3          = 1e-30;                    // conversion from Ang.^3 to m^3
    double      nAv                 = 6.02214076e23;            // in units of 1/mol
    double      gmol_to_kg          = 1e-3/nAv;                 // conversion from g/mol to kg                                                         
    double      densConvertSI       = gmol_to_kg / Ang3_to_m3;  // convert dens from g/mol/Ang^3 to kg/m^3

    double      boxLen[3];

    std::vector<double> mass( nBins, 0. );  // mass vector with all elements set to zero

    // Open file stream for input file 
    std::ifstream fin( inputFile );
    if ( !fin.is_open() ) {
        std::cerr << "Error: Could not open input file: " << inputFile << "\n";
        return 1;
    }

    // Open file stream for output file 
    std::ofstream fout( outputFile );
    if ( !fout.is_open() ) {
        std::cerr << "Error: Could not open output file: " << outputFile << std::endl;
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

    std::cout << "The conversion factor is = " << densConvertSI << std::endl;
    std::cout << "The box length in x dir is = " << boxLen[0] << std::endl;
    std::cout << "The box length in y dir is = " << boxLen[1] << std::endl;
    std::cout << "The box length in z dir is = " << boxLen[2] << std::endl;

    const int EXPECTED_FRAMES = 250000;
    int totalFrameCount = 0;
    int totalBlockCount = 0;

    while ( totalFrameCount < EXPECTED_FRAMES ) {

        bool                    blockFirstFrame = true;
        int                     blockFrameCount = 0;

        while ( blockFrameCount < blockSize && totalFrameCount < EXPECTED_FRAMES ) 
        {
            std::string line;
            while (std::getline(fin, line)) 
            {
                if (line.find("ITEM: TIMESTEP") != std::string::npos) break;
            }

            if ( fin.eof() ) break;

            int timeStep;
            std::getline(fin, line); std::istringstream(line) >> timeStep;

            std::getline(fin, line); // NUMBER OF ATOMS
            int nAtoms;
            std::getline(fin, line); std::istringstream(line) >> nAtoms;

            std::getline(fin, line); // box bounds header
            for (int i = 0; i < 3; i++) std::getline(fin, line);

            std::getline(fin, line); // ATOMS header

            std::map<int, std::vector<AtomData>> atomsByMol;
            for ( int i = 0; i < nAtoms; i++ ) 
            {
                if ( !std::getline(fin, line) ) break;
                if ( line.empty() || line[0] == '#' ) 
                { i--; continue; }

                std::istringstream iss(line);
                AtomData atom;
                if (iss >> atom.id >> atom.mol >> atom.type
                        >> atom.x >> atom.y >> atom.z 
                        >> atom.ix >> atom.iy >> atom.iz) 
                {

                    int bin = static_cast<int>(std::abs( std::abs(atom.z) - std::abs(zHot) ) / binw );
                    if ( bin >= 0 && bin < nBins )
                    {
                        mass[bin] += getMass( atom.type );
                    }
                    else
                    {
                        std::cout << "atom.type = " << atom.type << "\n";
                        std::cout << "atom.z    = " << atom.z    << "\n";
                        std::cout << "zlo       = " << zHot       << "\n";
                        std::cout << "binw      = " << binw      << "\n";
                        std::cout << "bin       = " << bin       << "\n";
                        std::cout << " ABORTED !!! Incorrect assignment " << std::endl;
                        abort();
                    }
                }
            }

            blockFrameCount++;
            totalFrameCount++;

            if ( totalFrameCount % 1000 == 0 ) { std::cout << totalFrameCount << " completed " << std::endl; }
        }

        if ( blockFrameCount == 0 ) break;

        double  binVolume   = boxLen[0] * boxLen[1] * binw;         // in Ang.^3
        totalBlockCount++;
        int numFramesInBlock = blockFrameCount;   
        std::cout << "Block " << totalBlockCount << ": Processed " << numFramesInBlock << " frames; \n";

        fout << "# Density profile (individual atoms)\n"
             << "# Block           : " << totalBlockCount << "\n"
             << "# zlim_min        : " << zHot << " zlim_max: " << zCold << " nbins: " << nBins << "\n"
             << "# Bin volume      : " << binVolume << " (A^3, fixed xy_area: " << boxLen[0]*boxLen[1] << ")\n"
             << "# Columns         : bin_position [Ang] z_center density [kg/m^3] \n";

        for ( int bin = 0; bin < nBins; ++bin ){    // write in the direction of hot to cold, easy for post-processing

            double binPos;

            if      ( domainSide.compare("L") == 0 ) { binPos = zHot - ( bin + 0.5 ) * binw; }
            else if ( domainSide.compare("R") == 0 ) { binPos = zHot + ( bin + 0.5 ) * binw; }

            double dens   = mass[bin] / ( numFramesInBlock * binVolume );   // in g/mol/Ang^3

            double dens_si = dens * densConvertSI;

            // Set persistent formatting for the stream
            fout << std::scientific << std::setprecision(6);

            fout    << std::setw( 15 ) << binPos 
                    << std::setw( 15 ) << dens_si << std::endl;

            mass[bin] = 0.;     // reset mass[bin] to zero since this is not a running average
        }
    }

    fin.close();
    fout.close();

    std::cout << "Processed " << totalBlockCount << " blocks (" << totalFrameCount << " frames).\n";
    std::cout << "Output written to " << outputFile << "\n";

    return 0;
}
