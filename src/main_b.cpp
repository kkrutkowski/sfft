#include <iostream>
#include <iomanip>
#include <filesystem>
#include <vector>
#include <fstream>
#include <tuple>
#include <thread>
#include <atomic>
#include <chrono>
#include <locale>
#include <cmath>

#include "utils/fftgrid.hpp"
#include "periodogram_b.hpp"

std::atomic<float> progressValue(0);
std::atomic<int> timeLeft(0), filesComputed(0);

void printProgress(int max_progress, const std::chrono::steady_clock::time_point startTime) {
    while (true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        if (filesComputed.load() == max_progress) {
            return;
        }
        int currentFilesComputed = filesComputed.load();
        float progress = 100. * static_cast<float>(currentFilesComputed) / static_cast<float>(max_progress);
        std::chrono::duration<double> elapsedTime = std::chrono::steady_clock::now() - startTime;
        float localTimeLeft = static_cast<float>(elapsedTime.count()) * (float(max_progress) - float(currentFilesComputed)) / float(currentFilesComputed);
        progressValue.store(10 * progress);
        timeLeft.store(int(localTimeLeft));
        std::cout << std::fixed << std::setprecision(2) << "\r" << progress << "% complete, " << localTimeLeft << " seconds left." << std::flush;
    }
}

int main(int argc, char *argv[]){

if (argc < 6){return 1;}
std::locale::global(std::locale("C"));

//defines and calculates constants used for calculations
const std::string files_location = argv[3];

std::string argv1 = argv[1];
std::string argv2 = argv[2];

uint method = 5;
uint algorithm = 2;

if (argv1 == "0" || argv1 == "naive" || argv1 == "slow"){method = 0; argv1 = "naive";}
else if (argv1 == "1" || argv1 == "direct" || argv1 == "simd" || argv1 == "fma"){method = 1; argv1 = "fma";}
else if (argv1 == "2" || argv1 == "rec" || argv1 == "recursive"){method = 2; argv1 = "recursive";}
else if (argv1 == "3" || argv1 == "fft"|| argv1 == "sfft" || argv1 == "sparse"){method = 3; argv1 = "sfft";}
else if (argv1 == "4" || argv1 == "fasper" || argv1 == "fast"){method = 4; argv1 = "fasper";}

if (argv2 == "0" || argv2 == "rz" || argv2 == "rzt" || argv1 == "rayleigh"){algorithm = 0; argv2 = "Rayleigh's Z-test";}
else if (argv2 == "1" || argv2 == "ls" || argv2 == "gls" || argv1 == "generalized"){algorithm = 1; argv2 = "Generalized Lomb-Scargle Periodogram";}

float max_frequency_temp = 10.0;
if (argc > 4 && ((argv[4][0]) != '\0')){max_frequency_temp = std::stof(argv[4]);}
const float max_frequency = max_frequency_temp; //0

double resolution;
if (argc > 5 && ((argv[5][0]) != '\0')){resolution = std::stod(argv[5]);}

FFTGrid grid; grid.generate(max_frequency, uint(std::exp2(resolution)), method);

std::cout << "\n" "Directory location: " << files_location << "\n";
std::cout << "Algorithm: " << argv2 << "\n";
std::cout << "Method: " << argv1 << "\n";
std::cout << "Min frequency: " << 0.0f << "\n";
std::cout << "Max frequency: " << max_frequency << "\n";
std::cout << "Step size: " << grid.fstep << "\n";

const unsigned int no_steps = grid.freq.size();


std::cout << "Number of steps: " << no_steps << "\n";

//const int max_thread_number = std::thread::hardware_concurrency();

//creates files array
std::vector<std::string> files;
auto directory_iterator = std::filesystem::directory_iterator(files_location);
unsigned int file_count = 0;
for (auto& entry : directory_iterator)
{files.push_back(entry.path());++file_count;}

std::cout <<"Number of files in directory: " << file_count << "\n" << std::endl;
// for(unsigned int i=0; i < files.size(); i++) std::cout << files.at(i) << ','; //prints list of files

const std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

std::thread printThread;
    printThread = std::thread([&]() {
        printProgress(file_count, startTime);
    });

std::filesystem::path path = std::filesystem::path(files_location).parent_path();

std::string output_path = path.string() + "/";
switch (method) {
    case 0:
        output_path += "naive"; break;
    case 1:
        output_path += "fma"; break;
    case 2:
        output_path += "recursive"; break;
    case 3:
        output_path += "sfft"; break;
    case 4:
        output_path += "fasper"; break;

    default: std::cout << "Invalid method selected" << std::endl; break;
}

if (algorithm == 0) {output_path += "_rz_output.tsv";}
else if (algorithm == 1) {output_path += "_gls_output.tsv";}
else {std::cout << "Invalid algorithm selected" << std::endl;}

std::cout << "Output file path: " << output_path << std::endl;

std::ofstream output_file(output_path);
 output_file.close(); //creates and closes output file

std::ofstream out(output_path);
out << "file	frequency	period	power" << std::endl;

FFT fft;
fft.init(grid.size);

#pragma omp parallel for
for (unsigned int i = 0; i < file_count; i++) {
    auto [frequency, amplitude, max_power] = periodogram(grid, files[i], fft, algorithm, method);

        #pragma omp critical
        {// Enter critical section to write to the file
            out << files[i] << "\t" << std::fixed << std::setprecision(6) << frequency << "\t" << std::fixed << std::setprecision(4) << 1/frequency << "\t" << std::fixed << std::setprecision(3) << max_power << "\n";
        };

#pragma omp critical
{filesComputed += 1;}
}

std::cout << "\r" << "Complete" << std::endl;
std::cout << "Elapsed time: "  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()<< " ms" << std::endl;
printThread.join();

fft.free();

return 0;}
