#include <cmath>
#include <cstring>

#include "../utils/vertex.hpp"
#include "../utils/fftgrid.hpp"
#include "../utils/readout.hpp"
//#include "utils/convolution.hpp"

output_data gls_slow_b(const star &data, const FFTGrid &grid, FFT &fft) {
        output_data best_frequency;

   int n = data.x.size();

   float wsum=0, wsum_inv=0, Y=0;
   float *t = (float *) malloc(n * sizeof(float)), //single precision float representation of time
          *w = (float *) malloc(n * sizeof(float)),
          *wy = (float *) malloc(n * sizeof(float));

   float max_power = 0;
   uint i = 0; uint k = 0;
   float YY = 0;

   for (i=0; i<n; ++i) {
      /* weights */
      w[i] = 1 / (data.dy[i] * data.dy[i]);
      wsum += w[i];
      t[i] = float(data.x[i]);
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
      YY += w[i] * wy[i] * wy[i];
      wy[i] *= w[i];                /* attach weights */
   }

   for (k=1; k<grid.freq.size(); ++k){
      float C = 0;
      float S = 0;
      float YC = 0;
      float YS = 0;
      float CC = 0;
      float CS = 0;
      float SS = 0;
      float D = 0;
      float power;
      for (i=0; i<n; ++i) {
         float cosx = cos(2 * M_PI * grid.freq[k] * t[i]);
         float sinx = sin(2 * M_PI * grid.freq[k] * t[i]);
         C += w[i] * cosx;
         S += w[i] * sinx;
         YC += wy[i] * cosx;
         YS += wy[i] * sinx;
         CC += w[i] * cosx * cosx;   /* Eq. (13) */
         CS += w[i] * cosx * sinx;   /* Eq. (15) */
      }

         SS = 1. - CC;
         CC -= C * C;
         SS -= S * S;
         CS -= C * S;
         D = CC*SS - CS*CS;
         power = (SS*YC*YC + CC*YS*YS - 2.*CS*YC*YS) / (YY*D);

         if (power > best_frequency.power)
               {
                  best_frequency.power = power;
                  best_frequency.frequency = grid.freq[k];
                     //std::cout << max_power << " " << best_freq << std::endl;
               }}

free(t); free(w); free(wy);
return best_frequency;}
