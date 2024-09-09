#include <cmath>
#include <cstring>

#include "../utils/vertex.hpp"
#include "../utils/periodograms.hpp"
#include "../utils/fftgrid.hpp"
#include "../utils/readout.hpp"
//#include "utils/convolution.hpp"

output_data rayleigh_fft_b(const star &data, const FFTGrid &grid, int terms, FFT &fft, int &threadID) {
        output_data best_frequency;

   int n = data.x.size();

   float wsum=0, wsum_inv=0, Y=0, norm = 0;
   float *t = (float *) malloc(n * sizeof(float)), //single precision float representation of time
          *w = (float *) malloc(n * sizeof(float)),
          *wy = (float *) malloc(n * sizeof(float));

   float *powers = (float *) calloc(grid.freq.size(), sizeof(float));

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

   //norm /= float(n);


   for(i=1; i <= terms; i++) {

   //std::cout << n << std::endl;
   float* CS= fft.NDFT_I(t, wy, n, grid.fstep, threadID, float(i));
   for (k=0; k< int(grid.freq.size() - 1); ++k) {powers[k] += float((CS[2*k] * CS[2*k]) + (CS[(2*k) + 1] * CS[(2*k) + 1])) / norm;} // Calibrate 'Z' value in Rayleigh's Z-test'
   free(CS);}

   for (k=1; k<grid.freq.size(); ++k){if (powers[k] > best_frequency.power)
         {
                  best_frequency.power = powers[k];
                  best_frequency.frequency = grid.freq[k];
                     //std::cout << max_power << " " << best_freq << std::endl;
               }}

free(t); free(w); free(wy); free(powers);
return best_frequency;}
