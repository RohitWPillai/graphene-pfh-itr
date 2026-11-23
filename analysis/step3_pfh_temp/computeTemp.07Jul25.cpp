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
    double      nAv                 = 6.02214076e23;            // in units of 1/mol
    double      kB                  = 1.380649e-23;             // in units of J/K                                                            
    double      gmol_to_kg          = 1e-3/nAv;                 // conversion from g/mol to kg                                                         

    double      ang2_ps2_to_m2_s2   =   1e4;                                // ang^2 / fs^2 to m^2 / s^2
    double      to_joules           =   ang2_ps2_to_m2_s2 * gmol_to_kg;     // energy to joules 
    double      tempConvertSI       =   to_joules / ( 3. * kB );                                                                                

    double      boxLen[3];

    std::vector<double> ke( nBins, 0. );          // kinetic energy vector with all elements set to zero
    std::vector<int>    atom_count( nBins, 0);      // count all atoms in a bin                                           

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

    std::cout << "The factor to convert energy to temperature = " << tempConvertSI << std::endl;
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
                if ( iss >> atom.id >> atom.mol  >> atom.type
                         >> atom.x  >> atom.y    >> atom.z 
                         >> atom.ix >> atom.iy   >> atom.iz 
                         >> atom.vx >> atom.vy   >> atom.vz  ) 
                {

                    int bin = static_cast<int>(std::abs( std::abs(atom.z) - std::abs(zHot) ) / binw );
                    if ( bin >= 0 && bin < nBins )
                    {
                        double mass = getMass( atom.type );
                        double v2   =   atom.vx*atom.vx + 
                                        atom.vy*atom.vy + 
                                        atom.vz*atom.vz ;

                        double  Two_ke      = mass * v2; 
                        double  atom_temp   = Two_ke;

                        ke[bin] += Two_ke;
                        atom_count[bin] += 1;
                    }
                    else
                    {
                        std::cout << "atom.type = " << atom.type << "\n";
                        std::cout << "atom.z    = " << atom.z    << "\n";
                        std::cout << "zHot      = " << zHot       << "\n";
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

        totalBlockCount++;
        int numFramesInBlock = blockFrameCount;   
        std::cout << "Block " << totalBlockCount << ": Processed " << numFramesInBlock << " frames; \n";

        fout << "# Temperature profile (individual atoms)\n"
             << "# Block           : " << totalBlockCount << "\n"
             << "# zlim_min        : " << zHot << " zlim_max: " << zCold << " nbins: " << nBins << "\n"
             << "# Columns         : bin_position [Ang] z_center density [kg/m^3] \n";

        for ( int bin = 0; bin < nBins; ++bin ){ 

            double binPos;
            if      ( domainSide.compare("L") == 0 ) { binPos = zHot - ( bin + 0.5 ) * binw; }
            else if ( domainSide.compare("R") == 0 ) { binPos = zHot + ( bin + 0.5 ) * binw; }
            double temp = 0.;

            if ( atom_count[bin] > 0 ){

                double keAv = ke[bin] / (1. * atom_count[bin]);   
                temp = keAv * tempConvertSI;     // in K
            }

            // Set persistent formatting for the stream
            fout << std::setprecision( 4 );

            fout << std::setw( 15 ) << binPos 
                 << std::setw( 15 ) << temp << "\n";

            ke[bin]         = 0.;     // reset mass[bin] to zero since this is not a running average
            atom_count[bin] = 0;      // set atom count to zero 
        }
    }

    fin.close();
    fout.close();

    std::cout << "Processed " << totalBlockCount << " blocks (" << totalFrameCount << " frames).\n";
    std::cout << "Output written to " << outputFile << "\n";

    return 0;
}
