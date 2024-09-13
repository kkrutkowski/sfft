#include <cmath>
#include <cstring>
#include <iostream>

#include "../utils/vertex.hpp"
#include "../utils/fftgrid.hpp"
#include "../utils/readout.hpp"
//#include "utils/convolution.hpp"

output_data gls_fft_b(const star &data, const FFTGrid &grid, FFT &fft) {
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


   //std::cout << n << std::endl;
   float* YCYS= fft.NDFT_I(t, wy, n, grid.fstep, float(1));
   float* C2S2= fft.NDFT_I(t, w, n, grid.fstep, float(2));
   float* CS = fft.NDFT_I(t, w, n, grid.fstep, float(1));



   float power, T2, S2w, C2w, Cw, Sw, YC, YS, CC, SS;
   for (k=1; k< int(grid.freq.size() - 1); ++k) {

   T2 = (C2S2[2*k + 1] - (2 * CS[2*k + 1] * CS[2*k])) / (C2S2[2*k] - ((CS[2*k] * CS[2*k]) - (CS[2*k + 1] * CS[2*k + 1])));
   S2w = T2 / std::sqrt(1 + (T2 * T2));
   C2w = 1 / std::sqrt(1 + (T2 * T2));
   Cw = M_SQRT1_2 * std::sqrt(1 + C2w);
   Sw = std::copysign(M_SQRT1_2, S2w) * std::sqrt(1 - C2w);

   YC = (YCYS[2*k] * Cw) + (YCYS[2*k + 1] * Sw);
   YS = (YCYS[2*k+1] * Cw) - (YCYS[2*k] * Sw);
   CC = 0.5 * (1 + (C2S2[2*k] * C2w) + (C2S2[2*k+1] * S2w));
   SS = 0.5 * (1 - (C2S2[2*k] * C2w) - (C2S2[2*k+1] * S2w));

   CC -= ((CS[2*k] * Cw) + (CS[2*k+1] * Sw)) * ((CS[2*k] * Cw) + (CS[2*k+1] * Sw));
   SS -= ((CS[2*k+1] * Cw) - (CS[2*k] * Sw)) * ((CS[2*k+1] * Cw) - (CS[2*k] * Sw));

   power = (YC * YC / CC) + (YS * YS / SS);
   power /= YY;

   if (power > best_frequency.power)
         {
                  best_frequency.power = power;
                  best_frequency.frequency = grid.freq[k];
         }
   }

free(t); free(w); free(wy); free(CS); free(YCYS); free(C2S2);
return best_frequency;}
