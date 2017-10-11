//
// The code is free for any use
//

#ifndef FFT_FAST_H
#define FFT_FAST_H

#include <stdlib.h>
#include <math.h>
#include <omp.h>

void fft_reverse(size_t N, double *xr, double *xi);
void FFT_DIT_RN(const size_t N, double *xr, double *xi);
void FFT_DIT_RN_2(const size_t N, double *xr, double *xi);

#endif //FFT_FAST_H
