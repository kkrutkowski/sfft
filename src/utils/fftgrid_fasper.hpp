#ifndef FFTGRID_HPP
#define FFTGRID_HPP

//#include "/opt/intel/oneapi/2024.1/include/ipp/ipps.h"
#include <fftw3.h>
#include <vector>
#include <iostream>
#include <cmath>
#include "extirpolation.hpp"

/*

struct FFT {

    int order;
    int size;

    IppsFFTSpec_C_32fc *pSpec;

    Ipp8u* pMemSpec;
    Ipp8u* pMemInit;
    Ipp8u* pMemBuffer;

    int sizeSpec;
    int sizeInit;
    int sizeBuffer;

    // Member function to initialize the FFT struct
    void init(int order_loc) {

        order = 11;//order_loc;
        size = int(ldexp(double(1.0), order));

        pSpec = 0;

        pMemSpec = 0;
        pMemInit = 0;
        pMemBuffer = 0;

        sizeSpec   = 0;
        sizeInit   = 0;
        sizeBuffer = 0;

        // Query to get buffer sizes
        ippsFFTGetSize_C_32fc(order, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, &sizeSpec, &sizeInit, &sizeBuffer);

        // Alloc FFT buffer
        pMemSpec = (Ipp8u*) ippsMalloc_8u(sizeInit);
        pMemInit = (Ipp8u*) ippsMalloc_8u(sizeSpec);
        pMemBuffer = (Ipp8u*) ippsMalloc_8u(sizeBuffer);

        //Initialize DFT
        ippsFFTInit_C_32fc(&pSpec, order, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone, pMemSpec, pMemInit);
        if (sizeInit > 0) ippsFree(pMemInit);
    }

        // Member function to calculate nonuniform to uniform FFT
        double* NDFT_I(double *x, double *y, int n, double df, double FreqFactor = 1){

            df *= FreqFactor;
            double xmin=0;
            double* xnorm = new double[n];

            for (int i = 0; i < n; i++) {xnorm[i] = (x[i] - x[0]) * df;}
            for (int i = 0; i < n; i++) {xnorm[i] = (xnorm[i] - double(int(xnorm[i]))) * double(size);}
            double* result = extirpolate(xnorm, y, n, size);
            ippsFFTInv_CToC_32fc_I(reinterpret_cast<Ipp32fc*>(result), pSpec, pMemBuffer);
            delete[] xnorm;

        return result;}

        void free(){
            ippsFree(pMemSpec);
            ippsFree(pMemBuffer);
        }



};

*/


struct FFT {

    int order;
    int size;

    fftwf_plan p;

    // Member function to initialize the FFT struct
    void init(int order_loc) {

        order = order_loc;
        size = int(ldexp(double(1.0), order));
        double* arr;
        p = fftwf_plan_dft_1d(size, reinterpret_cast<fftwf_complex*>(arr) , reinterpret_cast<fftwf_complex*>(arr), FFTW_BACKWARD, FFTW_ESTIMATE);
    }

        // Member function to calculate nonuniform to uniform FFT
        float* NDFT_I(double *x, double *y, int n, double df, double FreqFactor = 1){

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
