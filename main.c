//
// The code is free for any use
//
// Two realizations of DIT RN FFT
//

#include <stdio.h>
#include "fast.h"

#define bits 22
#define n (1 << bits)

double xr[n];
double xi[n];

void print_x(const size_t N, double *ar, double *ai) {
    for (int i = 0; i < N; ++i) {
        printf("\n%f  %f\n", ar[i], ai[i]);
    }
}

int test() {
    omp_set_num_threads(1);

    xr[0] = 1;
    xr[1] = 2;
    xr[2] = 0;
    xr[3] = -1;
    xi[0] = 0;
    xi[1] = -1;
    xi[2] = -1;
    xi[3] = 2;

    fft_reverse(4, xr, xi);
    FFT_DIT_RN(4, xr, xi);
    print_x(4, xr, xi);

    xr[0] = 1;
    xr[1] = 2;
    xr[2] = 0;
    xr[3] = -1;
    xi[0] = 0;
    xi[1] = -1;
    xi[2] = -1;
    xi[3] = 2;

    fft_reverse(4, xr, xi);
    FFT_DIT_RN_2(4, xr, xi);
    print_x(4, xr, xi);

    return 0;
}

int main() {

    //test();

    double start;

    for (int t = 1; t <= 4; ++t) {
        omp_set_num_threads(t);

        for (int k = 0; k < 1; k++) {

            //srand(time(NULL));
            srand(12345);

            for (int i = 0; i < n; ++i) {
                xr[i] = (double) rand() / RAND_MAX;
                xi[i] = (double) rand() / RAND_MAX;
            }

            fft_reverse(n, xr, xi);

            start = omp_get_wtime();
            FFT_DIT_RN(n, xr, xi);
            double time1 = omp_get_wtime() - start;
            //print_x(1, xr, xi);

            srand(12345);

            for (int i = 0; i < n; ++i) {
                xr[i] = (double) rand() / RAND_MAX;
                xi[i] = (double) rand() / RAND_MAX;
            }

            fft_reverse(n, xr, xi);

            start = omp_get_wtime();
            FFT_DIT_RN_2(n, xr, xi);
            double time2 = omp_get_wtime() - start;
            //print_x(1, xr, xi);

            printf("threads=%d  time1=%f  time2=%f\n", t, time1, time2);

        }
    }

    return 0;
}