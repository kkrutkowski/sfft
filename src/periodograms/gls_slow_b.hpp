#include <cmath>
#include <cstring>

#include "../utils/vertex.hpp"
#include "../utils/fftgrid.hpp"
#include "../utils/readout.hpp"
//#include "utils/convolution.hpp"

output_data gls_slow_b(const star &data, const FFTGrid &grid, FFT &fft) {
        output_data best_frequency;

   int n = data.x.size();

   float wsum=0, wsum_inv=0, Y=0, norm = 0;
   float *t = (float *) malloc(n * sizeof(float)), //single precision float representation of time
          *w = (float *) malloc(n * sizeof(float)),
          *wy = (float *) malloc(n * sizeof(float));

   float max_power = 0;
   int i = 0; int k = 0;

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
      wy[i] *= w[i];                /* attach weights */
      norm += wy[i] * wy[i];
   }

   for (k=1; k<grid.freq.size(); ++k){
      float C = 0;
      float S = 0;
      float power;
      for (i=0; i<n; ++i) {
         C += wy[i] * cos(2 * M_PI * grid.freq[k] * t[i]);
         S += wy[i] * sin(2 * M_PI * grid.freq[k] * t[i]);
      }
         power = ((C*C) + (S*S)) / norm;

         if (power > best_frequency.power)
               {
                  best_frequency.power = power;
                  best_frequency.frequency = grid.freq[k];
                     //std::cout << max_power << " " << best_freq << std::endl;
               }}

free(t); free(w); free(wy);
return best_frequency;}
