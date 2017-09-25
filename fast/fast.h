//
// Created by ysi on 9/25/17.
//

#ifndef FFT_FAST_H
#define FFT_FAST_H

#include <stdlib.h>
#include <math.h>
#include <omp.h>

int fft(size_t N, double *re, double *im);

#endif //FFT_FAST_H
