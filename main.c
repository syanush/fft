//
// Created by ysi on 9/25/17.
//

#include <stdio.h>
#include <time.h>

#include "fast.h"

//#define max_bits 25
//#define N (1 << max_bits)

#define bits 22
#define n (1 << bits)

double input_real[n]; // = {1, 2, 0, -1};
double input_imag[n]; // = {0, -1, -1, 2};

//#define M 40000000
//double a[M];
//
//void test1() {
//#pragma omp parallel for
//    for(int i=0;i<M;i++) {
//        a[i] = cos(sin(i));
//    }
//}

int main() {
    int error_code;

    printf ( "  Number of processors available = %d\n", omp_get_num_procs () );
    printf ( "  Number of threads =              %d\n", omp_get_max_threads () );
    omp_set_num_threads(1);
//#define M 40000000
//double a[M];
//
//void test1() {
//#pragma omp parallel for
//    for(int i=0;i<M;i++) {
//        a[i] = cos(sin(i));
//    }
//}

void tut1() {
    int tsum=0;
#pragma omp parallel
    tsum += omp_get_thread_num();
    printf("Sum is %d   ",tsum);
}

#define nmax 500000000
void tut2() {
    const int n=500000000;
    double dx = 1.0/nmax;
    double sum = 0;

    int i;

#pragma omp parallel for reduction(+:sum)
    for(i=0;i<nmax;++i) {
        double x=i*dx;
        sum += dx*sqrt(1-x*x);
    }

    printf("delta = %e    ", abs(M_PI/4-sum));
}

int main() {

    for (int t = 1; t <= 4; ++t) {
        omp_set_num_threads(t);

        clock_t start = clock();

        tut2();

        clock_t end = clock();
        float seconds = (float) (end - start) / CLOCKS_PER_SEC;
        printf("   %f sec\n", seconds);
    }




    //int error_code;

//    printf ( "  Number of processors available = %d\n", omp_get_num_procs () );
//    printf ( "  Number of threads =              %d\n", omp_get_max_threads () );
//    omp_set_num_threads(1);
//#pragma omp parallel
//    {
//#pragma omp master
//        printf("  Using %i threads.\n", omp_get_num_threads());
//    }

    //srand(time(NULL));

//    for (int t = 1; t <= 4; ++t) {
//        omp_set_num_threads(t);
//        srand(12345);
//
//        for (int i = 0; i < n; ++i) {
//            input_real[i] = (double) rand() / RAND_MAX;
//            input_imag[i] = (double) rand() / RAND_MAX;
//        }
//
//        clock_t start = clock();
//
//        error_code = fft(n, input_real, input_imag);
//        //test1();
//
//        clock_t end = clock();
//        float seconds = (float) (end - start) / CLOCKS_PER_SEC;
//
//        printf("%f sec\n", seconds);
//    }

    //for(int i=0; i<n; ++i) {
    //    printf("%8.2f %8.2f\n", input_real[i], input_imag[i]);
    //}

    //if (error_code) {
    //    printf("Error code: %d\n", error_code);
    //}

    return 0;
}