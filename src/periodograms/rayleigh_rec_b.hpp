#include <cmath>
#include <cstring>
#include <malloc.h>

#include "../utils/vertex.hpp"
#include "../utils/fftgrid.hpp"
#include "../utils/readout.hpp"
#include "../../include/mintrig.hpp"
//#include "utils/convolution.hpp"

output_data rayleigh_rec_b(const star &data, const FFTGrid &grid, FFT &fft) {
        output_data best_frequency;

   int n = data.x.size();
   uint n_simd = uint(n + (-n % 8) + 8);
   uint n_iter = n_simd / 8;

   //std::cout << n << "\t" << n_simd << std::endl; // ok


   float wsum=0, wsum_inv=0, Y=0, norm = 0;
   float *t = (float *) aligned_alloc(32, n_simd * sizeof(float)), //single precision float representation of time
          *w = (float *) aligned_alloc(32, n_simd * sizeof(float)),
          *wy = (float *) aligned_alloc(32, n_simd * sizeof(float)),
          *cosx = (float *) aligned_alloc(32, n_simd * sizeof(float)),
          *cosdx = (float *) aligned_alloc(32, n_simd * sizeof(float)),
          *sinx = (float *) aligned_alloc(32, n_simd * sizeof(float)),
          *sindx = (float *) aligned_alloc(32, n_simd * sizeof(float));

   float max_power = 0;
   uint i = 0; uint k = 0;

   for (i=0; i<n; ++i) {
      /* weights */
      w[i] = 1 / (data.dy[i] * data.dy[i]);
      wsum += w[i];
      t[i] = float(data.x[i]) - float(data.x[n/2]);
   }


   wsum_inv = 1 / wsum;

   for (i=0; i<n; ++i) {
      /* mean */
      w[i] *= wsum_inv;                 /* normalised weights */
      Y += w[i] * data.y[i];             /* Eq. (7) */
   }

   for (i=0; i<n; ++i) {
      /* variance */
      wy[i] = data.y[i] - Y;             /* Subtract weighted mean */
      wy[i] *= w[i];                /* attach weights */
      norm += wy[i] * wy[i];
   }

   for (i = 0; i < n; i++){
      cosdx[i] = cos(2 * M_PI * grid.fstep * t[i]);
      sindx[i] = sin(2 * M_PI * grid.fstep * t[i]);
   }

   for (i = n; i < n_simd; i++){
      cosdx[i] = 0;
      sindx[i] = 0;
   }

   for (i = n; i < n_simd; ++i){w[i] = 0; t[i] = 0; wy[i] = 0;}

   __m256* wy_simd = reinterpret_cast<__m256*>(wy);
   __m256* t_simd = reinterpret_cast<__m256*>(t);

   __m256* cosdx_simd = reinterpret_cast<__m256*>(cosdx);
   __m256* sindx_simd = reinterpret_cast<__m256*>(sindx);

   __m256* cosx_simd = reinterpret_cast<__m256*>(cosx);
   __m256* sinx_simd = reinterpret_cast<__m256*>(sinx);

   for (k=1; k<grid.freq.size(); ++k){
      float C = 0;
      float S = 0;
      __m256 tmp;
      float tmp_float;

      __m256 C_simd = _mm256_set1_ps(0.0f);
      __m256 S_simd = _mm256_set1_ps(0.0f);

      float power;

      // Two loops below are FMA equivalents of the loop above
      /*
      for (i=0; i<n; ++i) {
         C += wy[i] * cos(2 * M_PI * grid.freq[k] * t[i]);
         S += wy[i] * sin(2 * M_PI * grid.freq[k] * t[i]);
      }
      */


      if (k % 64 == 1){
      for (i=0; i<n_iter; i++){
         FTA::sincos_ps( _mm256_mul_ps(t_simd[i], _mm256_mul_ps(_mm256_set1_ps(grid.freq[k]), _mm256_set1_ps(twopi))), &sinx_simd[i], &cosx_simd[i]);
      }
      }

      for (i=0; i<n_iter; i++){
         C_simd = _mm256_fmadd_ps(wy_simd[i], cosx_simd[i], C_simd);
         S_simd = _mm256_fmadd_ps(wy_simd[i], sinx_simd[i], S_simd);
      }

      for (i=0; i<8; i++){
         C += reinterpret_cast<float*>(&C_simd)[i];
         S += reinterpret_cast<float*>(&S_simd)[i];
      }



      if (k % 64 != 0){
         for (i=0; i<n_iter; i++){
            //tmp = cosx_simd[i] * cosdx_simd[i] - sinx_simd[i] * sindx_simd[i];
            tmp = _mm256_fmsub_ps(cosx_simd[i], cosdx_simd[i], _mm256_mul_ps(sinx_simd[i], sindx_simd[i]));
            //sinx_simd[i] = cosx_simd[i] * sindx_simd[i] + sinx_simd[i] * cosdx_simd[i];
            sinx_simd[i] = _mm256_fmadd_ps(cosx_simd[i], sindx_simd[i], _mm256_mul_ps(sinx_simd[i], cosdx_simd[i]));
            cosx_simd[i] = tmp;
         }
      }


         power = ((C*C) + (S*S)) / norm;

         if (power > best_frequency.power && power < n)
               {
                  best_frequency.power = power;
                  best_frequency.frequency = grid.freq[k];
                     //std::cout << max_power << " " << best_freq << std::endl;
               }}

free(t); free(w); free(wy);
return best_frequency;}
