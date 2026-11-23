#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <fftw3.h>

struct AtomData {
    double vx, vy, vz, fx, fy, fz;
};

fftw_complex* alloc_fftw_complex(size_t n);
void free_fftw_complex(fftw_complex* p);

void computeAtomContribution(const std::vector<AtomData>& atomTS, int L,
                             fftw_complex* fft_in, fftw_complex* fft_out, fftw_complex* fft_res,
                             fftw_plan plan, std::vector<double>& seg_x,
                             std::vector<double>& seg_y, std::vector<double>& seg_z, int outSize,
                             const std::vector<double>& W);

#endif // UTILS_H
