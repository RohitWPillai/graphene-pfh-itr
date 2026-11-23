#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>

int main(int argc, char** argv) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <folder> <suffixes> <winlen> <startSeg> <endSeg> <outFile>\n";
        std::cerr << "<suffixes> is comma-separated, e.g., suffix1,suffix2\n";
        return 1;
    }
    std::string folder = argv[1];
    std::string suffixes_str = argv[2];
    int winlen = std::atoi( argv[3] );
    int startSeg = std::atoi( argv[4] );
    int endSeg = std::atoi( argv[5] );
    std::string outFile = argv[6];
    // Parse suffixes
    std::vector<std::string> suffixes;
    std::stringstream ss(suffixes_str);
    std::string suf;
    while (std::getline(ss, suf, ',')) {
        suffixes.push_back(suf);
    }
    unsigned int outSize= winlen/2 + 1;
    std::string fname, line;
    std::vector<double> freq( outSize,0.);
    std::vector<double> pxcum(outSize, 0.);
    std::vector<double> pycum(outSize, 0.);
    std::vector<double> pzcum(outSize, 0.);
    unsigned int nSeg=0;
    for (const auto& suffix : suffixes) {
        for( int i=startSeg; i<=endSeg; i++){
            std::cout << "Reading seg " << i << " for suffix " << suffix << std::endl;
            fname= folder + "/" + suffix + "_seg" + std::to_string(i) + ".txt";
            std::ifstream fin(fname);
            if (!fin.is_open()) {
                std::cerr << "Error: Could not open file: " << fname << std::endl;
                exit(1);
            }
            unsigned lnum=0;
            double px, py, pz;
            while ( std::getline(fin, line) ){
                if ( line.find("#") == 0 ) continue;
                std::istringstream iss(line);
                iss >> freq[lnum] >> px >> py >> pz;
                // std::cout << freq[lnum] << ", " << px << ", " << py << ", " << pz << std::endl;
                pxcum[lnum] += px;
                pycum[lnum] += py;
                pzcum[lnum] += pz;
                lnum++;
            }
            if ( lnum != outSize ){ std::cout << "Total lines " << lnum << " does not match the fft out size" << outSize << std::endl; abort(); }
            nSeg += 1;
            fin.close();
            /*
            while (std::getline(fin, line) && line.find("#") ==0 ) continue;
            {
                std::cout << line << std::endl;
            }
            */
        }
    }
    std::cout << "Total Segments processed = " << nSeg << std::endl;
    // Final output
    std::ofstream outs( outFile );
    outs << "# Frequency [Hz] , Px [eV/s], Py [eV/s], Pz[eV/s], Ptot[eV/s]" << std::endl;
    for ( int i=0; i < outSize; i++ ){
        outs << std::setw( 15 ) << std::scientific << freq[i] <<
                    std::setw( 15 ) << std::scientific << pxcum[i] / nSeg <<
                    std::setw( 15 ) << std::scientific << pycum[i] / nSeg <<
                    std::setw( 15 ) << std::scientific << pzcum[i] / nSeg <<
                    std::setw( 15 ) << std::scientific << (pxcum[i] + pycum[i] + pzcum[i]) / nSeg << std::endl;
    }
    return 0;
}
