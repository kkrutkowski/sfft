#include <ipps.h>
#include <cmath>
#include <iostream>
#include <iomanip>

struct FFT {

    int n;  //number of data points on which convolution is transformed

    IppsDFTSpec_R_32f *pDFTSpec;

    Ipp8u* pDFTInitBuf;
    Ipp8u* pDFTWorkBuf;

    Ipp32f* pSrc;
    Ipp32f* pDst;

    int sizeDFTSpec;
    int sizeDFTInitBuf;
    int sizeDFTWorkBuf;

    Ipp32f* kernel;

    // Member function to initialize the FFT struct
    void init(int size) {

        n = size;

            pDFTSpec=0;
        // Allocate buffers
            pSrc = ippsMalloc_32f(n);
            pDst = ippsMalloc_32f(n);

        // Query to get buffer sizes
            ippsDFTGetSize_C_32f(n, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &sizeDFTSpec, &sizeDFTInitBuf, &sizeDFTWorkBuf);

        // Alloc DFT buffers
            pDFTSpec    = (IppsDFTSpec_R_32f*)ippsMalloc_8u(sizeDFTSpec);
            pDFTInitBuf = ippsMalloc_8u(sizeDFTInitBuf);
            pDFTWorkBuf = ippsMalloc_8u(sizeDFTWorkBuf);

        // Initialize DFT
            ippsDFTInit_R_32f(n, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, pDFTSpec, pDFTInitBuf);
            if (pDFTInitBuf) ippsFree(pDFTInitBuf);
    }

    void generate_kernel(float a) {

        kernel = ippsMalloc_32f(n); // allocate memory for kernel
        double wsum;
        for (int i = 0; i <= n/2; i++) {kernel[i] = exp(float(-16) * float(i * i)/(a * a * float(n)));}
        for (int i = n-1; i > n/2; i--) {kernel[i] = kernel[n-i];}
        //for (int i = 0; i < n; i++) {std::cout << std::fixed << std::setprecision(3)  << kernel[i] << "\t";} //print the weights
        for (int i = 0; i < n; i++) {wsum += kernel[i];}
        for (int i = 0; i < n; i++) {kernel[i] /= wsum;}
        //for (int i = 0; i < n; i++) {std::cout << std::fixed << std::setprecision(3)  << kernel[i] << "\t";} //print normalized weights
        ippsDFTFwd_RToPack_32f(kernel, kernel, pDFTSpec, pDFTWorkBuf); //transform the kernel
        //ippsDFTInv_PackToR_32f(kernel, kernel, pDFTSpec, pDFTWorkBuf); //retrieve weights
        //for (int i = 0; i < n; i++) {std::cout << std::fixed << std::setprecision(3)  << kernel[i] << "\t";} //print the kernel
    }

    void generate_kernel(){generate_kernel(1);}

    void convolve(float* in, float* out) {

        ippsDFTFwd_RToPack_32f(reinterpret_cast<Ipp32f*>(in), reinterpret_cast<Ipp32f*>(out), pDFTSpec, pDFTWorkBuf);
        ippsMulPack_32f_I(kernel, reinterpret_cast<Ipp32f*>(out), n);
        ippsDFTInv_PackToR_32f(reinterpret_cast<Ipp32f*>(out), reinterpret_cast<Ipp32f*>(out), pDFTSpec, pDFTWorkBuf);
    }

        void convolve(float* y) {convolve(y, y);}

        // Destructor
    ~FFT() {
        if (pDFTSpec) ippsFree(pDFTSpec);
        if (pDFTWorkBuf) ippsFree(pDFTWorkBuf);
        if (pSrc) ippsFree(pSrc);
        if (pDst) ippsFree(pDst);
        if (kernel) ippsFree(kernel);
    }

};
