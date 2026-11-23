#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
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

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <timeSeries_file> <overallOutputFile> [dt]\n\n"
                  << "This program reads a time-series file containing PFH molecule atomic velocities\n"
                  << "in the format:\n"
                  << "   timeStep molID vx_c1 vy_c1 vz_c1 ... vx_f14 vy_f14 vz_f14 (62 columns)\n\n"
                  << "Each block begins with '# Block' and ends with a blank line.\n"
                  << "Output: Frequency, AvgPowerX, AvgPowerY, AvgPowerZ, Sum (X+Y+Z).\n";
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outFileAvg = argv[2];
    double dt = (argc >= 4) ? std::atof(argv[3]) : 1.0;
    if (dt <= 0) {
        std::cerr << "Error: dt must be positive.\n";
        return 1;
    }

    std::ifstream fin(inputFile);
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open input file " << inputFile << "\n";
        return 1;
    }

    std::vector<std::vector<double>> overallSum;
    int validBlockCount = 0;
    unsigned int global_fft_len = 0;
    std::string line;
    int blockNumber = 0;

    while (std::getline(fin, line)) {
        if (line.substr(0, 7) != "# Block") continue;

        blockNumber++;
        std::cout << "Processing Block " << blockNumber << "...\n";

        std::map<int, std::vector<std::vector<double>>> velocityByMol;
        int blockFrameCount = 0;
        while (std::getline(fin, line)) {
            if (line.empty()) break;
            if (line[0] == '#') continue;
            std::istringstream iss(line);
            int timeStep, molID;
            std::vector<double> velocities(60);
            if (!(iss >> timeStep >> molID)) {
                std::cerr << "Warning: could not parse timeStep or molID in line: " << line << "\n";
                continue;
            }
            for (int i = 0; i < 60; ++i) {
                if (!(iss >> velocities[i])) {
                    std::cerr << "Warning: incomplete velocity data at column " << i + 3 << " in line: " << line << "\n";
                    break;
                }
            }
            velocityByMol[molID].push_back(velocities);
            blockFrameCount++;
        }
        std::cout << "  Block " << blockNumber << " contains " << blockFrameCount
                  << " frames and " << velocityByMol.size() << " molecules.\n";

        if (velocityByMol.empty()) {
            std::cerr << "Block " << blockNumber << " contains no data. Skipping.\n";
            continue;
        }

        unsigned int N = 0;
        for (const auto &p : velocityByMol)
            N = std::max(N, static_cast<unsigned int>(p.second.size()));
        if (N == 0) continue;

        unsigned int fft_len = nextPowerOfTwo(2 * N);
        if (global_fft_len == 0) {
            global_fft_len = fft_len;
            overallSum.resize(global_fft_len, std::vector<double>(3, 0.0));
        } else if (fft_len != global_fft_len) {
            std::cerr << "Warning: Block " << blockNumber << " has fft_len " << fft_len
                      << " differing from global_fft_len " << global_fft_len << ". Skipping.\n";
            continue;
        }

        std::vector<std::vector<double>> blockSum(fft_len, std::vector<double>(3, 0.0));
        int validMolCount = 0;
        for (const auto &p : velocityByMol) {
            const auto &samples = p.second;
            if (samples.size() != N) continue;
            validMolCount++;

            std::vector<std::vector<double>> total_power(fft_len, std::vector<double>(3, 0.0));
            for (int a = 0; a < 20; ++a) {
                double mass = (a < 6) ? 12.017 : 18.998; // Carbon: 12.017, Fluorine: 18.998
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
                blockSum[i][0] += total_power[i][0];
                blockSum[i][1] += total_power[i][1];
                blockSum[i][2] += total_power[i][2];
            }
        }

        if (validMolCount == 0) {
            std::cerr << "Block " << blockNumber << " has 0 valid molecules. Skipping.\n";
            continue;
        }

        std::vector<std::vector<double>> blockAvg(fft_len, std::vector<double>(3, 0.0));
        for (unsigned int i = 0; i < fft_len; ++i) {
            blockAvg[i][0] = blockSum[i][0] / validMolCount;
            blockAvg[i][1] = blockSum[i][1] / validMolCount;
            blockAvg[i][2] = blockSum[i][2] / validMolCount;
        }

        for (unsigned int i = 0; i < fft_len; ++i) {
            overallSum[i][0] += blockAvg[i][0];
            overallSum[i][1] += blockAvg[i][1];
            overallSum[i][2] += blockAvg[i][2];
        }
        validBlockCount++;
        std::cout << "  Block " << blockNumber << ": Valid molecules = " << validMolCount << "\n";
    }

    fin.close();

    std::ofstream foutAvg(outFileAvg);
    if (!foutAvg.is_open()) {
        std::cerr << "Error: could not open output file " << outFileAvg << "\n";
        return 1;
    }

    unsigned int nFreq = global_fft_len / 2 + 1;
    foutAvg << "# zero padded to fft_len = " << global_fft_len << " (padded by a factor of 2)\n";
    foutAvg << "# Columns: Frequency\tAvgPowerX\tAvgPowerY\tAvgPowerZ\tSumPower\n";
    if (validBlockCount > 0) {
        double df = 1.0 / (global_fft_len * dt);
        for (unsigned int i = 0; i < nFreq; ++i) {
            double freq = i * df;
            double avgX = overallSum[i][0] / validBlockCount;
            double avgY = overallSum[i][1] / validBlockCount;
            double avgZ = overallSum[i][2] / validBlockCount;
            double sumPower = avgX + avgY + avgZ;
            foutAvg << freq << "\t" << avgX << "\t" << avgY << "\t" << avgZ << "\t" << sumPower << "\n";
        }
    } else {
        foutAvg << "# No valid blocks were processed.\n";
    }
    foutAvg.close();

    std::cout << "Processed " << validBlockCount << " blocks with valid data.\n";
    std::cout << "Overall average power spectrum written to " << outFileAvg << "\n";

    return 0;
}
