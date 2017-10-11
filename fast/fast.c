//
// The code is free for any use
//

#include <stdio.h>
#include "fast.h"

int reverse(int N, int n) {
    int log2N = (int) log2(N);
    int j, p = 0;
    for (j = 1; j <= log2N; ++j) {
        if (n & (1 << (log2N - j)))
            p |= 1 << (j - 1);
    }
    return p;
}


void fft_reverse(size_t N, double *xr, double *xi) {
    const size_t N2 = N >> 1;

#pragma omp parallel for
    for (int i = 0; i < N2; ++i) {
        int reverted = reverse(N, i);
        double tmp = xr[i];
        xr[i] = xr[reverted];
        xr[reverted] = tmp;
        tmp = xi[i];
        xi[i] = xi[reverted];
        xi[reverted] = tmp;
    }
}


void fft_twiddle(const size_t N_2, double *wre, double *wim) {
    const double phi = -M_PI / (N_2);
#pragma omp parallel for
    for (size_t i = 0; i < N_2; ++i) {
        double alpha = phi * i;
        wre[i] = cos(alpha);
        wim[i] = sin(alpha);
    }
}

/*
 * Chu, E., George, A. Inside the FFT black box
 */
void FFT_DIT_RN(const size_t N, double *xr, double *xi) {
    const size_t N_2 = N >> 1;
    double *wre = (double *) malloc(N_2 * sizeof(double));
    double *wim = (double *) malloc(N_2 * sizeof(double));

    fft_twiddle(N_2, wre, wim);

    int pairs = N >> 1; // in a group
    int groups = 1;

    while (groups < N) {
        int gapToNextPair = groups << 1;
        int gapToLastPair = gapToNextPair * (pairs - 1);

#pragma omp parallel for
        for (int k = 0; k < groups; ++k) {
            int j = k;
            int jlast = k + gapToLastPair;
            int jtwiddle = k * pairs;

            double wr = wre[jtwiddle];
            double wi = wim[jtwiddle];
            while (j <= jlast) {
                int jd = j + groups;
                double tr = xr[jd];
                double ti = xi[jd];
                double vr = wr * tr - wi * ti;
                double vi = wr * ti + wi * tr;
                xr[jd] = xr[j] - vr;
                xi[jd] = xi[j] - vi;
                xr[j] += vr;
                xi[j] += vi;
                j += gapToNextPair;
            }
        }
        pairs >>= 1;
        groups <<= 1;
    }
}

void FFT_DIT_RN_2(const size_t N, double *xr, double *xi) {
    const size_t N_2 = N >> 1;
    double *wre = (double *) malloc(N_2 * sizeof(double));
    double *wim = (double *) malloc(N_2 * sizeof(double));

    fft_twiddle(N_2, wre, wim);

    size_t groups = 1;
    size_t pairs = N_2;
    const int log2N = (int) log2(N);

    for (int j = 0; j < log2N; ++j) {

#pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            if (!(i & groups)) {

                // for a half of i's
                // mask1: 10101010, groups=1=0001, pairs=N/2
                // mask2: 11001100, groups=2=0010, pairs=N/4
                // mask3: 11110000, groups=4=0100, pairs=N/8
                //   ...
                // for every pair (i, id) apply the butterfly
                // u = x[i]
                // v = W*x[id]
                // x[i] = u + v
                // x[id] = u - v

                double ur = xr[i];
                double ui = xi[i];
                int itwiddle = (i * pairs) % N_2;

                double wr = wre[itwiddle];
                double wi = wim[itwiddle];

                int id = i + groups;

                double tr = xr[id];
                double ti = xi[id];
                double vr = wr * tr - wi * ti;
                double vi = wr * ti + wi * tr;

                xr[i] = ur + vr;
                xi[i] = ui + vi;
                xr[id] = ur - vr;
                xi[id] = ui - vi;
            }
        }
        groups <<= 1;
        pairs >>= 1;
    }
    free(wre);
    free(wim);
}
