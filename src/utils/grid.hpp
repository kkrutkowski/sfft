#ifndef GRID_HPP
#define GRID_HPP

#include <vector>

struct Grid{
    std::vector<double> freq;
    double fstep;

    void generate(double f_min, double f_max, int res){
        fstep = 1;
        for (int i = 0; i < res; i++){fstep *= 0.5;}
        for (double freq_temp = 0; freq_temp <= f_max; freq_temp += fstep){if (freq_temp >= f_min){freq.push_back(freq_temp);}}
        if (freq.size() % 2 != 0){freq.push_back(freq[freq.size() - 1] +fstep);}
    }
};


#endif
