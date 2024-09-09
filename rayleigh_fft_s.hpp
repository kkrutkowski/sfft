#include <cmath>
#include <cstring>

#include "utils/vertex.hpp"
#include "utils/periodograms.hpp"
#include "utils/fftgrid.hpp"
#include "utils/readout.hpp"
//#include "utils/convolution.hpp"

void rayleigh_fft_s(const star &data, const FFTGrid &grid, double* powers, int terms) {

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

   //norm /= float(n);

   FFT fft;
   fft.init(grid.order);

   int threadID = 0;

   for(i=1; i <= terms; i++) {
      float* CS = fft.NDFT_I(t, wy, n, grid.fstep, threadID, float(i));
      for (k=0; k<grid.freq.size() - 1; ++k) {powers[k] += double((CS[2*k] * CS[2*k]) + (CS[(2*k) + 1] * CS[(2*k) + 1])) / norm;} // 'Z' w teście Rayleigha
      fft.free();
      //std::cout << grid.freq[k] << std::endl;
      //std::cout << powers[k] << std::endl;
      free(CS);
   }

free(t); free(w); free(wy);
return;}
