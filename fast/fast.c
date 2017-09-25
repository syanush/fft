#include <stdio.h>
//#include <time.h>
//#include <pthread.h>
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

//void fft_work(size_t m1, size_t m2, double *re, double *im, double *wre, double *wim, size_t b, size_t a) {
//    for (size_t i = m1; i < m2; ++i) {
//        if (!(i & b)) {
//            int in = i + b;
//            double r1 = re[i];
//            double i1 = im[i];
//            int index = (i * a) % (b * a);
//            double r2 = wre[index] * re[in] - wim[index] * im[in];
//            double i2 = wre[index] * im[in] + wim[index] * re[in];
//            re[i] = r1 + r2;
//            im[i] = i1 + i2;
//            re[in] = r1 - r2;
//            im[in] = i1 - i2;
//        }
//    }
//}




//int thread_count;

//struct params {
//    size_t m1;
//    size_t m2;
//    double *re;
//    double *im;
//    double *wre;
//    double *wim;
//    size_t b;
//    size_t a;
//};

//void* fft_pthreads(void *args) {
//    struct params *fargs = (struct params*)args;
//    size_t m1 = fargs->m1;
//    size_t m2 = fargs->m2;
//    double *re = fargs->re;
//    double *im = fargs->im;
//    double *wre = fargs->wre;
//    double *wim = fargs->wim;
//    size_t b = fargs->b;
//    size_t a = fargs->a;
//
//    for (size_t i = m1; i < m2; ++i) {
//        if (!(i & b)) {
//            int in = i + b;
//            double r1 = re[i];
//            double i1 = im[i];
//            int index = (i * a) % (b * a);
//            double r2 = wre[index] * re[in] - wim[index] * im[in];
//            double i2 = wre[index] * im[in] + wim[index] * re[in];
//            re[i] = r1 + r2;
//            im[i] = i1 + i2;
//            re[in] = r1 - r2;
//            im[in] = i1 - i2;
//        }
//    }
//
//
//    return NULL;
//}

int fft(const size_t N, double *re, double *im) {

    double *wre, *wim;
    size_t b;
    size_t a;

//    long thread;
//    pthread_t *thread_handles;
//    thread_count = 2;
//    thread_handles = malloc (thread_count*sizeof(pthread_t));

//    struct params args0;
//    args0.re = re;
//    args0.im = im;
//    args0.wre = wre;
//    args0.wim = wim;
//
//    struct params args1;
//    args1.re = re;
//    args1.im = im;
//    args1.wre = wre;
//    args1.wim = wim;


    size_t i, j;

    b = 1;
    a = N >> 1;

    wre = (double *) malloc(a * sizeof(double));
    wim = (double *) malloc(a * sizeof(double));

    const int log2N = (int)log2(N);
    double phi = -M_2_PI / N;

    // предвычисление множителей
//#pragma omp parallel shared ( wre, wim, phi ) private ( i )
//#pragma omp for
    for (i = 0; i < N / 2; ++i) {
        double alpha = phi * i;
        wre[i] = cos(alpha);
        wim[i] = sin(alpha);

        //wre[1] * wre[i - 1] - wim[1] * wim[i - 1];
        //wim[i] = wre[1] * wim[i - 1] - wim[1] * wre[i - 1];
    }

    //clock_t start = clock();

    fft_reverse(N, re, im);

    //clock_t end = clock();
    //float seconds = (float) (end - start) / CLOCKS_PER_SEC;

    //printf("%f sec\n", seconds);




    for (j = 0; j < log2N; ++j) {
        // Main loop paralelization

//        struct params args;
//
//#pragma omp parallel shared(args, re, im, wre, wim, b, a)
//        {
//
//            args.re = re;
//            args.im = im;
//            args.wre = wre;
//            args.wim = wim;
//            args.b = b;
//            args.a = a;
//
//#pragma omp sections
//            {
//#pragma omp section
//                {
//                    //fft_work(0, N / 2, re, im, wre, wim, b, a);
//                    args.m1 = 0;
//                    args.m2 = N/2;
//                    fft_pthreads(&args);
//                }
//#pragma omp section
//                {
//                    //fft_work(N / 2, N, re, im, wre, wim, b, a);
//                    args.m1 = N/2;
//                    args.m2 = N;
//                    fft_pthreads(&args);
//                }
//            }
//        }

        ////



//        args0.b = b;
//        args0.a = a;
//        args0.m1 = 0;
//        args0.m2 = N;
//
//        args1.b = b;
//        args1.a = a;
//        args1.m1 = N/2;
//        args1.m2 = N;
//
//        fft_pthreads(&args0);

        //pthread_create(&thread_handles[0], NULL, fft_pthreads, &args0);
        //pthread_create(&thread_handles[1], NULL, fft_pthreads, &args1);
        //pthread_join(thread_handles[0], NULL);
        //pthread_join(thread_handles[1], NULL);



        ////

#pragma omp parallel shared ( re, im, wre, wim, b, a, j ) private (i)
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

    //free(thread_handles);

    return 0;
}



