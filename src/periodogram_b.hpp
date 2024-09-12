#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm> //for sorting
#include <tuple>
#include <cmath>

#include "periodograms/rayleigh_slow_b.hpp"
#include "periodograms/rayleigh_simd_b.hpp"
#include "periodograms/rayleigh_rec_b.hpp"
#include "periodograms/rayleigh_fft_b.hpp"
#include "periodograms/rayleigh_fasper_b.hpp"

#include "periodograms/gls_slow_b.hpp"
#include "periodograms/gls_simd_b.hpp"

#include "utils/readout.hpp"

std::tuple<float, float, float> periodogram(FFTGrid &grid, std::filesystem::path in_file, FFT &fft,uint algorithm, uint method){

	star data; data.read(in_file);
	output_data best_frequency;

	if(algorithm == 0){
		if (method == 0){best_frequency = rayleigh_slow_b(data, grid, fft);}
		else if (method == 1){best_frequency = rayleigh_simd_b(data, grid, fft);}
		else if (method == 2){best_frequency = rayleigh_rec_b(data, grid, fft);}
		else if (method == 3){best_frequency = rayleigh_fft_b(data, grid, fft);}
		else if (method == 4){best_frequency = rayleigh_fasper_b(data, grid, fft);}
		else {std::cout << "Invalid method selected" << std::endl;}
	}
	else {
		if (method == 0){best_frequency = gls_slow_b(data, grid, fft);}
		else if (method == 1){best_frequency = gls_simd_b(data, grid, fft);}
	}

	float powers_average = best_frequency.sum_of_powers / double(grid.freq.size()); //calculates average power for the input data

	std::tuple<double, float, float> output_tuple = std::make_tuple(best_frequency.frequency, best_frequency.amplitude, best_frequency.power);
return output_tuple;
}
