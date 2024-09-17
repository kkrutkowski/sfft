#include <cmath>
#include <cstring>
#include <malloc.h>

#include "../utils/vertex.hpp"
#include "../utils/fftgrid.hpp"
#include "../utils/readout.hpp"
#include "../../include/fma_trig.hpp"

output_data rayleigh_simd_b(const star &data, const FFTGrid &grid, FFT &fft) {
        output_data best_frequency;

   int n = data.x.size();
   uint n_simd = uint(n + (-n % 8) + 8);
   uint n_iter = n_simd / 8;

   //std::cout << n << "\t" << n_simd << std::endl; // ok


   float wsum=0, wsum_inv=0, Y=0, norm = 0;
   float *t = (float *) aligned_alloc(32, n_simd * sizeof(float)), //single precision float representation of time
          *w = (float *) aligned_alloc(32, n_simd * sizeof(float)),
          *wy = (float *) aligned_alloc(32, n_simd * sizeof(float));

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

   for (i = n; i < n_simd; ++i){w[i] = 0; t[i] = 0; wy[i] = 0;}

   __m256* wy_simd = reinterpret_cast<__m256*>(wy);
   __m256* t_simd = reinterpret_cast<__m256*>(t);

   for (k=1; k<grid.freq.size(); ++k){
      float C = 0;
      float S = 0;
      __m256 C_simd = _mm256_set1_ps(0.0f);
      __m256 S_simd = _mm256_set1_ps(0.0f);
      __m256 C_temp;
      __m256 S_temp;
      float power;

      // Two loops below are FMA equivalents of the loop above
      /*
      for (i=0; i<n; ++i) {
         C += wy[i] * cos(2 * M_PI * grid.freq[k] * t[i]);
         S += wy[i] * sin(2 * M_PI * grid.freq[k] * t[i]);
      }
      */


      for (i=0; i<n_iter; i++){
         FTA::sincos_2pi_ps( _mm256_mul_ps(t_simd[i], _mm256_set1_ps(grid.freq[k])), &S_temp, &C_temp);
         C_simd = _mm256_fmadd_ps(wy_simd[i], C_temp, C_simd);
         S_simd = _mm256_fmadd_ps(wy_simd[i], S_temp, S_simd);
      }

      for (i=0; i<8; i++){
         C += reinterpret_cast<float*>(&C_simd)[i];
         S += reinterpret_cast<float*>(&S_simd)[i];
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
