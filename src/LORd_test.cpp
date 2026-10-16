#include <chrono>
#include <iostream>
#include "LORd.hpp"

int main(){


    LORd_model_params_t params{
        205e9, 
        200e9, 
        7.85e-9, 
        7.85e-9, 
        1.00e-3, 
        3.23e-3, 
        1.96e-3, 
        1.90e-3, 
        3.09e-5, 
        3.80e-3, 
        8.0e6, 
        8.0e6, 
        8.0e6, 
        8.0e6, 
        0.060, 
        0.250, 
        0.050, 
        0.90, 
        0.62, 
        1.80, 
        7330.4
    };
    LORd_model model;

    constexpr int N = 1000000000;
    using clock = std::chrono::steady_clock;

    // Warm-up so caches and branch predictors are primed
    model.make_H0(params);

    auto t0 = clock::now();
   for (int i = 0; i < N; ++i) {
       asm volatile("" : "+m"(params) : : "memory");
       model.make_H0(params);
       asm volatile("" : : "g"(&model.H0) : "memory");
   }
    auto t1 = clock::now();

    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / N;

    std::cout << model.H0 << "\n";
    std::cout << "make_H0: " << ns << " ns per call (avg over " << N << " calls)\n";
    return 0;
}