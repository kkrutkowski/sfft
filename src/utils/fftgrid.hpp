#ifndef FFTGRID_HPP
#define FFTGRID_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include <fftw3.h>

#include "extirpolation.hpp"


struct FFT {

    int order;
    int size;
    float* arr;
    //const float cosh_tab[6];

    fftwf_plan p;

    // Member function to initialize the FFT struct

    void init(int order_loc) {

        order = order_loc;
        size = int(ldexp(double(1.0), order));
        p = fftwf_plan_dft_r2c_1d(2 * size, arr, reinterpret_cast<fftwf_complex*>(arr), FFTW_ESTIMATE);
    }


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

            /*
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
            */

            ::free(xnorm);
            fftwf_execute_dft_r2c(p, fftgrid, reinterpret_cast<fftwf_complex*>(fftgrid));
            //fftwf_execute(p);
            //::free(fftgrid);

            //for (int i = 0; i < int(2 * size) + 2; i++){fftgrid[i] *= std::sqrt(cosh(float(i) / (4 * ((float(size)) + 1))));}

        return fftgrid;}


        // Member function to calculate nonuniform to uniform FFT
        float* NDFT_I_fasper(float *x, float *y, int n, double df, int &threadID, double FreqFactor = 1){

            df *= 2.0 * FreqFactor;
            double xmin=0;
            double* xnorm = new double[n];

            for (int i = 0; i < n; i++) {xnorm[i] = (x[i] - x[0]) * df;}
            for (int i = 0; i < n; i++) {xnorm[i] = (xnorm[i] - double(int(xnorm[i]))) * double(size);}
            //std::cout << n << " | " << size << std::endl;
            float* result = extirpolate(xnorm, y, n, size);
            //for (int i = 0; i < 2 * size; i+=2) {std::cout << result[i] << ", " << std::flush;}

            //ippsFFTInv_CToC_32fc_I(, pSpec, pMemBuffer);
            fftwf_execute_dft(p, reinterpret_cast<fftwf_complex*>(result), reinterpret_cast<fftwf_complex*>(result));
            //for (int i = 0; i < 2 * size; i+=2) {std::cout << result[i] * result[i] + result[i+1] * result[i+1] << " " << std::flush;}
            delete[] xnorm;

        return result;}


        void free(){
            //fftwf_print_plan(p);
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
