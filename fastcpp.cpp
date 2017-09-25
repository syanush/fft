#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <random>
//#include <ctime>
#include <climits>
#include <cassert>
#include <omp.h>

using namespace std;

template<typename T>
T reverse(T n, size_t b = sizeof(T) * CHAR_BIT) {
    assert(b <= std::numeric_limits<T>::digits);

    T rv = 0;

    for (size_t i = 0; i < b; ++i, n >>= 1) {
        rv = (rv << 1) | (n & 0x01);
    }

    return rv;
}

void print_vec(vector<complex<double>> a) {
    for (int i = 0; i < a.size(); ++i) {
        cout << a[i] << endl;
    }
    cout << endl;
}

vector<complex<double>> bit_reverse_copy(vector<complex<double>> a) {
    size_t n = a.size();
    vector<complex<double>> A(n);
    for (int k = 0; k < n; ++k) {

        //r = rev(k);
        size_t k0 = k;
        int r = 0;
        for (int s = 0; s < int(log2(n)); ++s, k0 >>= 1) {
            r = (r << 1) | (k0 & 0x01);
            //;
        }
        //cout << k << " " << r << endl;
        A[r] = a[k];
    }
    return A;
}

vector<complex<double>> fft_iter(vector<complex<double>> a) {
    vector<complex<double>> A = bit_reverse_copy(a);
    size_t n = a.size();
    //cout << "n="<<n<<endl;

    clock_t tic, toc;
    double elapsed_secs;

    for (int s = 1; s <= int(log2(n)); ++s) {
        //cout << "s=" << s << endl;
        int m = 1 << s;
        //cout << "m="<<m<<endl;
        complex<double> omega_m = polar(1.0, 2 * M_PI / m);

        //tic = clock();

        for (int k = 0; k < n; k += m) { // тут можно распараллелить!!! общая память!!!
            complex<double> omega(1, 0);
            for (int j = 0; j < m / 2; ++j) {
                complex<double> t = omega * A[k + j + m / 2];
                complex<double> u = A[k + j];
                A[k + j] = u + t;
                A[k + j + m / 2] = u - t;
                omega *= omega_m;
            }
        }

        //toc = clock();
        //elapsed_secs = double(toc - tic) / CLOCKS_PER_SEC;
        //cout << elapsed_secs << " sec" << endl;

    }
    return A;
}


vector<complex<double>> fft_iter_parallel(vector<complex<double>> a) {
    vector<complex<double>> A = bit_reverse_copy(a);
    size_t n = a.size();
    //cout << "n="<<n<<endl;

    clock_t tic, toc;
    double elapsed_secs;

    auto threads = omp_get_num_threads();

    for (int s = 1; s <= int(log2(n)); ++s) {
        //cout << "s=" << s << endl;
        int m = 1 << s;
        int m2 = m >> 1;
        //cout << "m="<<m<<endl;
        complex<double> omega_m = polar(1.0, 2 * M_PI / m);

        //tic = clock();

        //int k;
//#pragma omp parallel shared(A) private(k)
        int k;
#pragma omp parallel for default(none) shared(A,n,m,m2,omega_m) private(k)
            for (k = 0; k < n; k += m) { // тут можно распараллелить!!! общая память!!!
                complex<double> omega(1, 0);
                for (int j = 0; j < m2; ++j) {
                    complex<double> t = omega * A[k + j + m2];
                    complex<double> u = A[k + j];
                    A[k + j] = u + t;
                    A[k + j + m2] = u - t;
                    omega *= omega_m;
                }
            }
//        }

        //toc = clock();
        //elapsed_secs = double(toc - tic) / CLOCKS_PER_SEC;
        //cout << elapsed_secs << " sec" << endl;

    }
    return A;
}


vector<complex<double>> fft(vector<complex<double>> &a) {
    size_t n = a.size();
    if (n == 1)
        return a;
    complex<double> omega_n = polar(1.0, 2 * M_PI / n);
    complex<double> omega(1, 0);
    size_t n2 = n / 2;
    vector<complex<double>> a0(n2), a1(n2);
    for (int i = 0; i < n2; ++i) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }
    vector<complex<double>> y0 = fft(a0);
    vector<complex<double>> y1 = fft(a1);
    vector<complex<double>> y(n);
    for (int k = 0; k < n2; ++k) {
        complex<double> t = omega * y1[k];
        y[k] = y0[k] + t;
        y[k + n2] = y0[k] - t;
        omega *= omega_n;
    }
    return y;
}

vector<complex<double>> fft_parallel(vector<complex<double>> &a, int count) {
    size_t n = a.size();
    if (n == 1)
        return a;
    complex<double> omega_n = polar(1.0, 2 * M_PI / n);
    complex<double> omega(1, 0);
    size_t n2 = n / 2;
    vector<complex<double>> a0(n2), a1(n2);
    vector<complex<double>> y0, y1, y(n);

//#pragma omp parallel for shared(a,a0,a1)
    for (int i = 0; i < n2; ++i) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }

    if (count > 0) {
#pragma omp parallel shared(y0,y1)
        {
            y0 = fft_parallel(a0, count-1);
            y1 = fft_parallel(a1, count-1);
        }
    } else {
        y0 = fft_parallel(a0, count);
        y1 = fft_parallel(a1, count);
    }


//#pragma omp parallel for shared(y,y0,y1)
    for (int k = 0; k < n2; ++k) {
        complex<double> t = omega * y1[k];
        y[k] = y0[k] + t;
        y[k + n2] = y0[k] - t;
        omega *= omega_n;
    }

    return y;
}

vector<complex<double>> dft(vector<complex<double>> a) {
    size_t n = a.size();
    double alpha = 2 * M_PI / n;
    vector<complex<double>> y(n);
    for (int k = 0; k < n; ++k) {
        y[k] = {0, 0};
        for (int j = 0; j < n; ++j) {
            double arg = alpha * ((k * j) % n); // sign corresponds to the sign of the phase
            y[k] += polar(1.0, arg) * a[j];
        }
    }
    return y;
}

int main() {

    cout << fixed << setprecision(5);

    random_device rdev;
    mt19937 rgen(rdev());
    uniform_real_distribution<> dis(0, 1);

    const int bitcount = 20;
    const int n = 1 << bitcount; // 2^20

    vector<complex<double>> a(n);

    for (int i = 0; i < n; ++i) {
        a[i] = {dis(rgen), dis(rgen)};
    }

    clock_t tic, toc;
    double elapsed_secs;

    tic = clock();
    vector<complex<double>> y1 = fft(a);
    toc = clock();
    elapsed_secs = double(toc - tic) / CLOCKS_PER_SEC;
    //print_vec(y);
    cout << elapsed_secs << " sec" << endl;

    tic = clock();
    vector<complex<double>> y2 = fft_parallel(a,2);
    toc = clock();
    elapsed_secs = double(toc - tic) / CLOCKS_PER_SEC;
    //print_vec(y);
    cout << elapsed_secs << " sec" << endl;

    double err{0};
    for (int i = 0; i < n; ++i) {
        double delta = fabs(real(y1[i]) - real(y2[i]));
        if (err < delta)
            err = delta;
    }

    cout << "Error: " << err << endl;

//    auto nthreads = omp_get_num_threads();
//    cout << nthreads << endl;
//
//    int tid;
//
///* Fork a team of threads giving them their own copies of variables */
//#pragma omp parallel private(nthreads, tid)
//    {
//
//        /* Obtain thread number */
//        tid = omp_get_thread_num();
//        printf("Hello World from thread = %d\n", tid);
//
//        /* Only master thread does this */
//        if (tid == 0)
//        {
//            nthreads = omp_get_num_threads();
//            printf("Number of threads = %d\n", nthreads);
//        }
//
//    }  /* All threads join master thread and disband */

    return 0;
}