#include "utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

// Alloc for FFT
fftw_complex* alloc_fftw_complex(size_t n) {
    fftw_complex* p = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
    if (!p) {
        std::cerr << "Out of memory\n";
        std::exit(1);
    }
    return p;
}

// Free for FFT
void free_fftw_complex(fftw_complex* p) {
    fftw_free(p);
}

/*-------------------------------------
  Compute contribution from each atom
-------------------------------------*/
void computeAtomContribution(const std::vector<AtomData>& atomTS,
                             int L,
                             fftw_complex* fft_in,
                             fftw_complex* fft_out,
                             fftw_complex* fft_res,
                             fftw_plan plan,
                             std::vector<double>& seg_x,
                             std::vector<double>& seg_y,
                             std::vector<double>& seg_z,
                             int outSize,
                             const std::vector<double>& W) {
    // X component
    for (int t = 0; t < L; ++t) {
        fft_in[t][0] = atomTS[t].fx * W[t];
        fft_in[t][1] = 0.0;
    }
    fftw_execute(plan);
    for (int t = 0; t < L; ++t) {
        fft_res[t][0] = fft_out[t][0];
        fft_res[t][1] = fft_out[t][1];
    }
    for (int t = 0; t < L; ++t) {
        fft_in[t][0] = atomTS[t].vx * W[t];
        fft_in[t][1] = 0.0;
    }
    fftw_execute(plan);
    for (int k = 0; k < outSize; ++k) {
        seg_x[k] += fft_res[k][0] * fft_out[k][0] + fft_res[k][1] * fft_out[k][1];
    }

    // Y component
    for (int t = 0; t < L; ++t) {
        fft_in[t][0] = atomTS[t].fy * W[t];
        fft_in[t][1] = 0.0;
    }
    fftw_execute(plan);
    for (int t = 0; t < L; ++t) {
        fft_res[t][0] = fft_out[t][0];
        fft_res[t][1] = fft_out[t][1];
    }
    for (int t = 0; t < L; ++t) {
        fft_in[t][0] = atomTS[t].vy * W[t];
        fft_in[t][1] = 0.0;
    }
    fftw_execute(plan);
    for (int k = 0; k < outSize; ++k) {
        seg_y[k] += fft_res[k][0] * fft_out[k][0] + fft_res[k][1] * fft_out[k][1];
    }

    // Z component
    for (int t = 0; t < L; ++t) {
        fft_in[t][0] = atomTS[t].fz * W[t];
        fft_in[t][1] = 0.0;
    }
    fftw_execute(plan);
    for (int t = 0; t < L; ++t) {
        fft_res[t][0] = fft_out[t][0];
        fft_res[t][1] = fft_out[t][1];
    }
    for (int t = 0; t < L; ++t) {
        fft_in[t][0] = atomTS[t].vz * W[t];
        fft_in[t][1] = 0.0;
    }
    fftw_execute(plan);
    for (int k = 0; k < outSize; ++k) {
        seg_z[k] += fft_res[k][0] * fft_out[k][0] + fft_res[k][1] * fft_out[k][1];
    }
}
