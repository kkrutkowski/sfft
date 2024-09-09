#include <cmath>
#include <iostream>
#include <cstring>
#include <x86simdsort.h>

#include "utils/vertex.hpp"
#include "utils/periodograms.hpp"
#include "utils/fftgrid.hpp"
#include "utils/readout.hpp"
//#include "utils/convolution.hpp"

output_data gls_b(const star &data, const Grid &grid, int terms) {
        output_data best_frequency;

   int n = data.x.size();

   double wsum=0, wsum_inv=0, Y=0, YY=0, C, S, YC, YS, CC, SS, CS, D, cos_amp, sin_amp, tmp;
   double *t = (double *) malloc(n * sizeof(double)), //single precision double representation of time
          *w = (double *) malloc(n * sizeof(double)),
          *wy = (double *) malloc(n * sizeof(double)),
          *wy_na = (double *) malloc(n * sizeof(double)),
          *cosx = (double *) malloc(n * sizeof(double)),
          *sinx = (double *) malloc(n * sizeof(double)),
          *cosdx = (double *) malloc(n * sizeof(double)),
          *sindx = (double *) malloc(n * sizeof(double)),

          *wy_temp = (double *) malloc(n * sizeof(double)),
          *y_temp = (double *) malloc(n * sizeof(double)),
          *sinx_temp = (double *) malloc(n * sizeof(double)),
          *cosx_temp = (double *) malloc(n * sizeof(double)),

          *powers = (double *) malloc(grid.freq.size() * sizeof(double));
   unsigned int i, k;

   double max_power = 0;

   for (i=0; i<n; ++i) {
      /* weights */
      w[i] = 1 / (data.dy[i] * data.dy[i]);
      wsum += w[i];
      t[i] = double(data.x[i]);
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
      wy_na[i] = wy[i];
      wy[i] *= w[i];                /* attach weights */

      /* Prepare trigonometric recurrences cos(dx)+i sin(dx) */
      cosdx[i] = cos(2 * M_PI * grid.fstep * t[i]);
      sindx[i] = sin(2 * M_PI * grid.fstep * t[i]);
   }

   for (k=0; k<grid.freq.size(); ++k) {
        powers[k] = 1;
        double amplitude = 0; double phase = 0;


         double C, D, S, YC, YS, CC, CS, tmp, cos_amp, sin_amp;
         for (uint i = 0; i < n; ++i) {
            if (k % int(exp2f(double(floor(log2f(double(1024 / terms))))) + 0.25) == 0) {   //roundabout way to call C++20 std::bit_floor without C++20
                  /* init/refresh recurrences to stop error propagation */
                  cosx[i] = cosf(2 * M_PI * grid.freq[k] * t[i]);
                  sinx[i] = sinf(2 * M_PI * grid.freq[k] * t[i]);
                  }}

            for (uint i = 0; i < n; ++i) {
                  wy_temp[i] = wy[i];
                  y_temp[i] = wy_na[i];
                  cosx_temp[i] = cosx[i];
                  sinx_temp[i] = sinx[i];
            }

            for (uint j = 0; j < terms; ++j){
               C = D = S = YC = YS = CC = CS = tmp = cos_amp = sin_amp = YY = 0;
               for (uint i = 0; i < n; ++i){
                  C += w[i] * cosx_temp[i];
                  S += w[i] * sinx_temp[i];

                  YY += w[i] * y_temp[i] * y_temp[i];   /* Eq. (10) */
                  YC += wy_temp[i] * cosx_temp[i];            /* Eq. (11) */
                  YS += wy_temp[i] * sinx_temp[i];

                  CC += w[i] * cosx_temp[i] * cosx_temp[i];
                  CS += w[i] * cosx_temp[i] * sinx_temp[i]; /* Eq. (12) */
                  }

                  /* power */
                  SS = 1. - CC;
                  CC -= C * C;
                  SS -= S * S;
                  CS -= C * S;
                  D = CC*SS - CS*CS;

                  cos_amp = (YC*SS-YS*CS) / D;
                  sin_amp = (YS*CC-YC*CS) / D;

                  //phase = std::atan2(sin_amp, cos_amp);
                  //amplitude = std::sqrt((cos_amp * cos_amp) + (sin_amp * sin_amp));
                  //std::cout << amplitude << std::endl;

                  powers[k] *= (1 - (SS*YC*YC + CC*YS*YS - 2.*CS*YC*YS) / (YY*D));  /* Eq. (5) in ZK09 */

                  for (uint i = 0; i < n; ++i){

                     //Fix the wy values for the normalization - slow implementation)
                        y_temp[i] -= cos_amp * cosx_temp[i];
                        y_temp[i] -= sin_amp * sinx_temp[i];
                        wy_temp[i] -= w[i] * YC * cosx_temp[i];
                        wy_temp[i] -= w[i] * YS * sinx_temp[i];
                     //std::cout << wy[i] << wy_temp[i] << std::endl;

                        /* increase freq for next harmonic */
                        if (j != terms - 1) {
                           tmp = cosx_temp[i] * cosx[i] - sinx_temp[i] * sinx[i];
                           sinx_temp[i] = cosx_temp[i] * sinx[i] + sinx_temp[i] * cosx[i];
                           cosx_temp[i] = tmp;
                        }
                     }
                  }

                  for (uint i = 0; i < n; ++i) {
                     /* increase freq for next loop */
                     tmp = cosx[i] * cosdx[i] - sinx[i] * sindx[i];
                     sinx[i] = cosx[i] * sindx[i] + sinx[i] * cosdx[i];
                     cosx[i] = tmp;
                  }



               powers[k] = 1 - powers[k];

               if (powers[k] > best_frequency.power && powers[k] < 5.0) {
                     best_frequency.power = powers[k];
                     best_frequency.frequency = grid.freq[k];
                     //std::cout << max_power << " " << best_freq << std::endl;
               }

     }

free(t); free(w); free(wy); free(cosdx); free(sindx); free(cosx); free(sinx); free(cosx_temp); free(sinx_temp); free(wy_temp); free(wy_na); free(y_temp); free(powers);
return best_frequency;}
