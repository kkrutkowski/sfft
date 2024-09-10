#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm> //for sorting
#include <tuple>
#include <cmath>

#include "periodograms/rayleigh_fft_b.hpp"

#include "utils/grid.hpp"
#include "utils/readout.hpp"

std::tuple<float, float, float> periodogram(FFTGrid &grid, std::filesystem::path in_file, FFT &fft, int &threadID, uint method){

	star data; data.read(in_file);
	output_data best_frequency;

	if (method == 4){best_frequency = rayleigh_fft_b(data, grid, fft, threadID);} //binned conditional nonuniformity, requires fixing weighting function

	float powers_average = best_frequency.sum_of_powers / double(grid.freq.size()); //calculates average power for the input data

	std::tuple<double, float, float> output_tuple = std::make_tuple(best_frequency.frequency, best_frequency.amplitude, best_frequency.power);
return output_tuple;
}
