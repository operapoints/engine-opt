#include <chrono>
#include <iostream>
#include "LORd.hpp"

int main(){


    LORd_model_params_t params{
        // ---------------- Young's moduli [Pa] ----------------
        2.05e11,   // EC    Compressor shaft Young's modulus (steel)
        2.05e11,   // ES    Center shaft Young's modulus (unused)
        2.05e11,   // ET    Turbine shaft Young's modulus (steel)

        // ---------------- Geometry [m] ----------------
        0.010,     // rC    Compressor shaft radius
        0.012,     // rS    Center shaft radius
        0.010,     // rT    Turbine shaft radius
        0.050,     // LC    Compressor overhang length
        0.150,     // LS    Center shaft length (bearing to bearing)
        0.050,     // LT    Turbine overhang length

        // ---------------- Mass properties ----------------
        7800.0,    // rhoS  Center shaft density [kg/m^3] (steel)
        0.45,      // mC    Compressor mass [kg] (~7 cm radius aluminum/titanium wheel)
        0.60,      // mT    Turbine mass [kg] (~7 cm radius Inconel wheel)
        5.0e-4,    // JdC   Compressor diametral moment [kg*m^2]
        7.5e-4,    // JdT   Turbine diametral moment [kg*m^2]
        8.8e-4,    // JpC   Compressor polar moment [kg*m^2]
        1.3e-3,    // JpT   Turbine polar moment [kg*m^2]

        // ---------------- Bearings [N/m] ----------------
        1.0e7,     // KbyC  Compressor end bearing y axis stiffness
        1.2e7,     // KbyT  Turbine end bearing y axis stiffness
        1.0e7,     // KbzC  Compressor end bearing z axis stiffness
        1.2e7,     // KbzT  Turbine end bearing z axis stiffness

        // ---------------- Operating condition ----------------
        7330.4     // Omega Shaft speed [rad/s] (70,000 rpm = 70000 * 2*pi / 60)
    };
    LORd_model model;

    constexpr int N = 100000;
    using clock = std::chrono::steady_clock;

    // Warm-up so caches and branch predictors are primed
    auto lambda = model.mode_eigvals(params);

    auto t0 = clock::now();
   for (int i = 0; i < N; ++i) {
       asm volatile("" : "+m"(params) : : "memory");
       model.mode_eigvals(params);
       asm volatile("" : : "g"(&model.H0) : "memory");
   }
    auto t1 = clock::now();

    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / N;

    std::cout << model.H0 << "\n";
    std::cout << "mode_eigvals: " << ns << " ns per call (avg over " << N << " calls)\n";
    std::cout << "Eigenvalues: \n";
    std::cout << lambda;
    return 0;


    // // Usage
    // // Instantiate rd model (done only once)
    // LORd_model model;
    // // Optimization loop{
    //     // Package relevant parameters in struct
    //     LORd_model_params_t params{...};
    //     // Solve system eigenvalues
    //     auto lambda = model.mode_eigvals(params);
    // // }

}