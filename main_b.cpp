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

std::locale::global(std::locale("C"));

//defines and calculates constants used for calculations
const std::string files_location = argv[2];

std::string argv1 = argv[1];

if (argc < 4){return 1;} // || argv[3][0] == '\0' - breaks gui

uint method = 2;

if (argv1 == "0" || argv1 == "slow"){method = 0; if (argv1 == "0"){argv1 = "slow";}}
else if (argv1 == "1" || argv1 == "direct" || argv1 == "simd" || argv1 == "fma"){method = 1; if (argv1 == "1" || argv1 == "direct" || argv1 == "simd"){argv1 = "fma";}}
else if (argv1 == "2" || argv1 == "rec" || argv1 == "recursive"){method = 2; if (argv1 == "2" || argv1 == "rec"){argv1 = "recursive";}}
else if (argv1 == "3" || argv1 == "sfft" || argv1 == "sparse"){method = 4; if (argv1 == "4" || argv1 == "sparse"){argv1 = "sfft";}}
else if (argv1 == "4" || argv1 == "fasper" || argv1 == "fast"){method = 5; if (argv1 == "5" || argv1 == "fast"){argv1 = "fasper";}}

float max_frequency_temp = 10.0;
if (argc > 3 && ((argv[3][0]) != '\0')){max_frequency_temp = std::stof(argv[3]);}
const float max_frequency = max_frequency_temp; //0

int terms = 1;
if (argc > 4 && ((argv[4][0]) != '\0')){terms = std::stoi(argv[4]);}

FFTGrid grid; grid.generate(max_frequency, int(12 + std::ceil(std::log2(terms * max_frequency))));

std::cout << "\n" "Directory location: " << files_location << "\n";
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

std::filesystem::path path = std::filesystem::path(files_location).parent_path(); std::string output_path = path.string() + "/rayleigh_output.tsv";
std::ofstream output_file(output_path);
 output_file.close(); //creates and closes output file

std::ofstream out(output_path);
out << "file	frequency	period	amplitude	power" << std::endl;

FFT fft;
fft.init(grid.order);

int threadID = -1;

#pragma omp parallel for private(threadID)
for (unsigned int i = 0; i < file_count; i++) {
    auto [frequency, amplitude, max_power] = periodogram(grid, files[i], fft, threadID, method);

        #pragma omp critical
        {// Enter critical section to write to the file
            out << files[i] << "\t" << std::fixed << std::setprecision(6) << frequency << "\t" << std::fixed << std::setprecision(4) << 1/frequency << "\t" << std::fixed << std::setprecision(3) << amplitude << "\t" << std::fixed << std::setprecision(3) << max_power << "\n";
        };

#pragma omp critical
{filesComputed += 1;}
}

std::cout << "\r" << "Complete" << std::endl;
std::cout << "Elapsed time: "  << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count()<< " ms" << std::endl;
printThread.join();

fft.free();

return 0;}
