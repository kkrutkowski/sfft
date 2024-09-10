#ifndef FFTGRID_HPP
#define FFTGRID_HPP

//#include "/opt/intel/oneapi/2024.1/include/ipp/ipps.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <fftw3.h>


struct FFT {

    int order;
    int size;
    float* arr;
    const float cosh_tab[6];

    fftwf_plan p;

    FFT() : cosh_tab {0.30250242, 0.32175471, 0.17278529, 1.96204223, 0.00298915, 0.99996111} {}

    // Member function to initialize the FFT struct
    void init(int order_loc) {

        order = order_loc;
        size = int(ldexp(double(1.0), order));
        p = fftwf_plan_dft_r2c_1d(2 * size, arr, reinterpret_cast<fftwf_complex*>(arr), FFTW_ESTIMATE);
    }
        float cosh(float x){  // computes approximate cosh(2x) for x between 0 and 1
            float result = cosh_tab[5];
            float temp = 1;
            for (int i=1; i <6 ; i++){
                temp *= x;
                result += temp * cosh_tab[6-i];
            }
        return result;}

        // Member function to calculate nonuniform to uniform FFT
        float* NDFT_I(float *x, float *y, int n, double df, int &threadID, float FreqFactor = 1, float TrimFactor = 1.0){

            float* fftgrid = (float*) calloc(int(2 * size) + 2, sizeof(float));
            //p = fftwf_plan_dft_r2c_1d(int(2 * size), fftgrid, reinterpret_cast<fftwf_complex*>(fftgrid), FFTW_ESTIMATE); //, fftwf_FORWARD

            float freqmax = 2 * size * df * FreqFactor;


            float* xnorm = (float*) calloc(n, sizeof(float));
            for (int i = 0; i < n; i++) {xnorm[i] = x[i] - x[0];}
            for (int i = 0; i < n; i++){
                float tmp = xnorm[i] * freqmax;
                fftgrid[int(tmp) % (2 * size)] += y[i] * (1 - (tmp - float(int(tmp))));
                fftgrid[(int(tmp) + 1) % (2 * size)] += y[i] * (tmp - float(int(tmp)));
            }



            if (threadID == -1) {
                #pragma omp critical
                {
                    p = fftwf_plan_dft_r2c_1d(2 * size, fftgrid, reinterpret_cast<fftwf_complex*>(fftgrid), FFTW_ESTIMATE);
                    threadID = 0;
                    // FFTW_ESTIMATE --> FFTW_MEASURE --> FFTW_PATIENT --> FFTW_EXHAUSTIVE
                    // FFTW_ESTIMATE ~ 5.0 s execution
                    // FFTW_MEASURE ~ 4.0 s + 4s evalutaion
                    // FFTW_MEASURE ~ 4.0 s + 50s evalutaion (?)
                }
            }

            ::free(xnorm);
            fftwf_execute_dft_r2c(p, fftgrid, reinterpret_cast<fftwf_complex*>(fftgrid));
            //fftwf_execute(p);
            //::free(fftgrid);

            //for (int i = 0; i < int(2 * size) + 2; i++){fftgrid[i] *= std::sqrt(cosh(float(i) / (4 * ((float(size)) + 1))));}

        return fftgrid;}

        void free(){
            fftwf_print_plan(p);
            fftwf_destroy_plan(p);
        }

};



struct FFTGrid{
    std::vector<double> freq;
    int order;
    int size;
    std::vector<FFT> fft;
    double fstep;

    void generate(double f_max, int order_input){ //fftgrid always starts at 0 to simplify calculations
        order = order_input;
        size = int(ldexp(double(1.0), order));
        freq.resize(size);
        for (int i = 0; i < size; i++){freq[i] = double(i) * double(f_max) / double(size);}
        fstep = freq[1];
    }


};


#endif
