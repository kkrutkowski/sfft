#ifndef EXTIRPOLATION_HPP
#define EXTIRPOLATION_HPP

#include <iostream>
//#include <vector>
//#include <numeric>
//#include <complex>
#include <cmath>
#include <algorithm>

double factorial(int n) {
    double result = 1.0;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

double* extirpolate(double *x, double *y, int size, int &N) {  // O(M) complexity
    const int M = 3;
    double max_x = 0;



    if (N <= 0) {
        for (int i = 0; i < size; i++){if (x[i] > max_x){max_x = x[i];}}
        N = int(max_x + 1.0 + (double(M) * 0.5));
    }

    //std::cout << N << std::endl;

    double*  result  =  (double*) calloc(2 * N, sizeof(double));
    int integers = 0;
    bool* IsInteger = (bool*)calloc(size, sizeof(bool));
    for (int i = 0; i < size; i++){if (x[i] == double(int(x[i]))){integers += 1; result[2*i] = y[i]; IsInteger[i] = 1;}} //else { if (IsInteger[i] == 1){integers -= 1;  IsInteger[i] = 0;}}}
    // np.add.at(result, x[integers].astype(int), y[integers]) - astropy
    //for (int i = 0; i < N; ++i) {std::cout << IsInteger[i] << " " << std::flush;}
    double* x_temp = (double*) calloc(size - integers, sizeof(double));
    double* y_temp = (double*) calloc(size - integers, sizeof(double));
    //std::cout << std::endl << std::endl;

    //std::cout << std::endl << std::endl << size - integers << std::endl << std::endl;


    int size_temp = 0;
    for (int i = 0; i < size; i++) {if (IsInteger[i] == 0) {x_temp[size_temp] = x[i]; y_temp[size_temp] = y[i]; size_temp++;};} //x, y = x[~integers], y[~integers]
    //std::cout << size_temp << std::endl;

    // for (int i = 0; i < size_temp; i++) {std::cout << x_temp[i] << " " << std::flush;} std::cout << std::endl << std::endl; OK
    // for (int i = 0; i < size_temp; i++) {std::cout << y_temp[i] << " " << std::flush;} std::cout << std::endl << std::endl; OK

    int* ilo =  (int*) calloc(size_temp, sizeof(int));
    for (int i = 0; i < size_temp; i++) {ilo[i] = std::clamp(int(x_temp[i] - (double(M) * 0.5)), 0, (2 * N) -M);}    // ilo = np.clip((x - M // 2).astype(int), 0, N - M)
    //for (int i = 0; i < size_temp; i++) {std::cout << ilo[i] << " " << std::flush;} std::cout << std::endl << std::endl; OK

    double denominator = factorial(M - 1);
    //std::cout << denominator << std::endl; // OK

    double* numerator = (double*) calloc(size_temp, sizeof(double));

    for (int i = 0; i < size_temp; i++) {    //numerator = y * np.prod(x - ilo - np.arange(M)[:, np.newaxis], 0)
        numerator[i] = y_temp[i];
        for (int j = 0; j < M; j++) {numerator[i] *= (x_temp[i] - ilo[i] - double(j));}
    }

    //for (int i = 0; i < size_temp; i++) {std::cout << numerator[i] << " " << std::flush;} std::cout << std::endl << std::endl; // OK

    for (int i = 0; i < M; i++){
        if (i > 0) {denominator *=  double(i) / double(i - M);}    //denominator *= j / (j - M)
        //for (int j = 0; j < size_temp; j++) {result[2* (ilo[j] + M - 1 - i)] += numerator[j] / (denominator * (x_temp[j] - double(ilo[j] + M - 1 - i)));} //np.add.at(result, ind, numerator / (denominator * (x - ind)))
        for (int j = 0; j < size_temp; j++) {result[(ilo[j] + M - 1 - i)] += numerator[j] / (denominator * (x_temp[j] - double(ilo[j] + M - 1 - i)));} //np.add.at(result, ind, numerator / (denominator * (x - ind)))
    }

    //for (int j = 0; j < N; j++) {std::cout << result[2*j] << " "<< std::flush;} // OK

    free(IsInteger); free(x_temp); free(y_temp); free(ilo); free(numerator);
return result;}

#endif
