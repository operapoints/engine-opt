// campbell.cpp
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "LORd.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(){
    LORd_model_params_t params{
        2.05e11, 2.05e11, 2.05e11,              // EC, ES (unused), ET   [Pa]
        0.010, 0.012, 0.010,                    // rC, rS, rT            [m]
        0.050, 0.150, 0.050,                    // LC, LS, LT            [m]
        7800.0,                                 // rhoS                  [kg/m^3]
        0.45, 0.60,                             // mC, mT                [kg]
        5.0e-4, 7.5e-4,                         // JdC, JdT              [kg*m^2]
        8.8e-4, 1.3e-3,                         // JpC, JpT              [kg*m^2]
        1.0e7, 1.2e7, 1.0e7, 1.2e7,             // KbyC, KbyT, KbzC, KbzT [N/m]
        0.0                                     // Omega (set in the sweep)
    };
    LORd_model model;

    // ---------------- Sweep configuration ----------------
    constexpr double rpm_min  = 0.0;
    constexpr double rpm_max  = 200000.0;
    constexpr int    n_points = 401;                  // 500 rpm steps
    constexpr double pi = 3.14159265358979323846;
    constexpr double rpm_to_rad_s  = 2.0 * pi / 60.0;
    constexpr double rad_s_to_krpm = 60.0 / (2.0 * pi) / 1000.0;
    constexpr int    n_eig = 16;                      // mode_eigvals returns a fixed 16x1

    using clock = std::chrono::steady_clock;

    std::vector<double> rpm_grid(n_points);
    for (int i = 0; i < n_points; ++i) {
        rpm_grid[i] = rpm_min + (rpm_max - rpm_min) * i / (n_points - 1);
    }

    // ---------------- Sweep helper ----------------
    // Runs the full speed sweep for one model variant, times it, and returns
    // scatter data (x = shaft speed [krpm], y = |Im(lambda)| [krpm]).
    struct SweepResult { std::vector<double> x, y; };

    auto run_sweep = [&](bool stiff_shaft, const std::string& name) {
        // Warm-up so caches and branch predictors are primed
        params.Omega = rpm_min * rpm_to_rad_s;
        model.mode_eigvals(params, stiff_shaft);

        // Preallocate so no allocation happens inside the timed region
        std::vector<double> imag_store(static_cast<size_t>(n_points) * n_eig);

        auto t0 = clock::now();
        for (int i = 0; i < n_points; ++i) {
            params.Omega = rpm_grid[i] * rpm_to_rad_s;
            auto lambda = model.mode_eigvals(params, stiff_shaft);
            for (int k = 0; k < n_eig; ++k) {
                imag_store[static_cast<size_t>(i) * n_eig + k] = lambda(k).imag();
            }
        }
        auto t1 = clock::now();

        const double total_us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        std::cout << name << ": evaluated " << n_points << " omega values ("
                  << n_eig << " eigenvalues each) in " << total_us / 1e3 << " ms\n"
                  << "  -> " << total_us / n_points << " us per omega value\n";

        SweepResult r;
        r.x.reserve(imag_store.size());
        r.y.reserve(imag_store.size());
        for (int i = 0; i < n_points; ++i) {
            for (int k = 0; k < n_eig; ++k) {
                r.x.push_back(rpm_grid[i] / 1000.0);
                r.y.push_back(std::abs(imag_store[static_cast<size_t>(i) * n_eig + k]) * rad_s_to_krpm);
            }
        }
        return r;
    };

    // ---------------- Run both variants ----------------
    SweepResult stiff = run_sweep(true,  "Stiff shaft    ");
    SweepResult flex  = run_sweep(false, "Non-stiff shaft");

    // 1X synchronous excitation line: frequency = shaft speed (y = x in krpm)
    std::vector<double> sync_x{rpm_min / 1000.0, rpm_max / 1000.0};
    std::vector<double> sync_y{rpm_min / 1000.0, rpm_max / 1000.0};

    // ---------------- Plot ----------------
    plt::figure_size(1000, 700);
    // Non-stiff drawn first (red, slightly larger) so the black stiff points stay visible on top
    plt::scatter(flex.x,  flex.y,  6.0, {{"c", "red"},   {"label", "Non-stiff shaft"}});
    plt::scatter(stiff.x, stiff.y, 3.0, {{"c", "black"}, {"label", "Stiff shaft"}});
    // 1X line changed to gray so it doesn't clash with the red data
    plt::plot(sync_x, sync_y, {{"color", "gray"}, {"linestyle", "--"}, {"label", "1X"}});
    plt::xlabel("Shaft speed [krpm]");
    plt::ylabel("Mode frequency [krpm]");
    plt::title("Campbell diagram: stiff vs. non-stiff shaft");
    plt::grid(true);
    plt::legend();
    plt::save("campbell.png");
    plt::show();
    return 0;
}