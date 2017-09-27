#include <stdio.h>
//#include <time.h>
#include "fast.h"

int reverse(int N, int n) {
    int j, p = 0;
    for (j = 1; j <= log2(N); ++j) {
        if (n & (1 << ((int) log2(N) - j)))
            p |= 1 << (j - 1);
    }
    return p;
}

/******************************************************************************/

void fft_reverse(size_t N, double *re, double *im) {
    double tmp;
    int i;
    const size_t N2 = N >> 1;
    const int log2N = (int)log2(N);

//#pragma omp parallel shared(re, im) private (i, tmp)
//#pragma omp for
    for (int i = 0; i < N2; ++i) {

        //int reverted = reverse(N, i);
        int j, reverted = 0;
        for (j = 1; j <= log2N; ++j) {
            if (i & (1 << (log2N - j)))
                reverted |= 1 << (j - 1);
        }

        /////
        tmp = re[i];
        re[i] = re[reverted];
        re[reverted] = tmp;
        tmp = im[i];
        im[i] = im[reverted];
        im[reverted] = tmp;
    }
}


int fft(const size_t N, double *re, double *im) {

    double *wre, *wim;
    size_t b;
    size_t a;

    size_t i, j;

    b = 1;
    a = N >> 1;

    wre = (double *) malloc(a * sizeof(double));
    wim = (double *) malloc(a * sizeof(double));

    const int log2N = (int)log2(N);
    double phi = -M_2_PI / N;

    // предвычисление множителей
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < N / 2; ++i) {
        double alpha = phi * i;
        wre[i] = cos(alpha);
        wim[i] = sin(alpha);

        //wre[1] * wre[i - 1] - wim[1] * wim[i - 1];
        //wim[i] = wre[1] * wim[i - 1] - wim[1] * wre[i - 1];
    }

    //clock_t start = clock();

    fft_reverse(N, re, im);

    for (j = 0; j < log2N; ++j) {

#pragma omp parallel
#pragma omp for
        for (i = 0; i < N; ++i) {
            if (!(i & b)) {
                int ib = i + b;
                double r1 = re[i];
                double i1 = im[i];
                int index = (i * a) % (b * a);
                double r2 = wre[index] * re[ib] - wim[index] * im[ib];
                double i2 = wre[index] * im[ib] + wim[index] * re[ib];

                re[i] = r1 + r2;
                im[i] = i1 + i2;
                re[ib] = r1 - r2;
                im[ib] = i1 - i2;
            }
        }


        b << 1;
        a >> 1;
    }

    free(wre);
    free(wim);

    return 0;
}



