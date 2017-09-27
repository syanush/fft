//
// Created by ysi on 9/25/17.
//

#include <stdio.h>
#include "fast.h"

#define bits 22
#define n (1 << bits)

double input_real[n]; // = {1, 2, 0, -1};
double input_imag[n]; // = {0, -1, -1, 2};

int main() {

    int error_code;

    for (int t = 1; t <= 6; ++t) {
        omp_set_num_threads(t);
        printf("threads=%d\n", t);
        for (int k = 0; k < 5; k++) {

            //srand(time(NULL));
            srand(12345);

            for (int i = 0; i < n; ++i) {
                input_real[i] = (double) rand() / RAND_MAX;
                input_imag[i] = (double) rand() / RAND_MAX;
            }

            double wall_start = omp_get_wtime();

            //tut2();
            error_code = fft(n, input_real, input_imag);

            double wall_end = omp_get_wtime();
            double sec = wall_end - wall_start;
            printf("   %f sec\n", sec);


            //clock_t start = clock();
            //clock_t end = clock();
            //float seconds = (float) (end - start) / CLOCKS_PER_SEC;
            //printf("   %f sec\n", seconds);
        }
    }


    return 0;
}