// Test the simple row physics independently of pagmo and stack propagation.
// Build with your implementation (do not #include its .cpp file):
// g++ -std=c++20 -Wall -Wextra -Wpedantic test_physics.cpp physics.cpp -o test_physics
// ./test_physics
//
// Assumptions: constant Cp/gamma, adiabatic rows, velocities normalized by
// sqrt(Cp * row inlet absolute total temperature). No dimensional inputs.
#include "rostack.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// Link this program with the source file defining YOUR function.
StagePerformance rostack_physics_simple(const RotorVariables&, double gamma);

namespace {
int checks = 0;
int failures = 0;

void near(const std::string& name, double actual, double expected)
{
    ++checks;
    const double tolerance = 1e-12 * std::max(1.0, std::abs(expected));
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        ++failures;
        std::cerr << std::setprecision(17) << "FAIL " << name
                  << ": got " << actual << ", expected " << expected << '\n';
    }
}

void check_case(const std::string& name, const RotorVariables& vars,
                const StagePerformance& expected, bool check_diffusion = true)
{
    try {
        const auto p = rostack_physics_simple(vars, 1.4);
        near(name + ".tau", p.tau, expected.tau);
        near(name + ".mach_rel_in", p.mach_rel_in, expected.mach_rel_in);
        near(name + ".mach_rel_out", p.mach_rel_out, expected.mach_rel_out);
        near(name + ".mach_abs_out", p.mach_abs_out, expected.mach_abs_out);
        if (check_diffusion) near(name + ".diffusion", p.diffusion, expected.diffusion);
    } catch (const std::exception& e) {
        ++checks;
        ++failures;
        std::cerr << "FAIL " << name << ": unexpected exception: " << e.what() << '\n';
    }
}

void expect_domain_error(const std::string& name, const RotorVariables& vars,
                         double gamma = 1.4)
{
    ++checks;
    try {
        (void)rostack_physics_simple(vars, gamma);
        ++failures;
        std::cerr << "FAIL " << name << ": expected std::domain_error\n";
    } catch (const std::domain_error&) {
        // Expected failure policy from the corrected physics implementation.
    } catch (const std::exception& e) {
        ++failures;
        std::cerr << "FAIL " << name << ": wrong exception: " << e.what() << '\n';
    }
}
} // namespace

int main()
{
    // RotorVariables: cth2, u1, cm1, k, lambda, sigma, nu, cth1.
    // Fixed reference values avoid reproducing the implementation in the test.
    // For example, compressor: tau=1.04, theta1=.95, theta2=.975,
    // W1^2=.18 and W2^2=C2^2=.13.
    const RotorVariables compressor{.2, .4, .3, 1, 1, 1, .5, .1};
    check_case("compressor", compressor,
               {1.04, .6882472016116854, .5773502691896257,
                .5773502691896257, .2680145445989607});

    // Fixed-zero blade speed: no shaft work, relative = absolute velocity.
    check_case("stator", {.1, 0, .3, 1, 1, 1, .5, .2},
               {1, .5895706767536022, .5129891760425771,
                .5129891760425771, .2616170297492780});

    // Unchanged velocity triangle: no work, equal inlet/outlet relative Mach,
    // and zero diffusion factor.
    check_case("unchanged", {.1, .4, .3, 1, 1, 1, .5, .1},
               {1, .6882472016116854, .6882472016116854,
                .5129891760425771, 0});

    // Check negative Euler work. Compressor DF is not validated for turbines.
    check_case("turbine", {.1, .4, .3, 1, 1, 1, .5, .3},
               {.92, .5241424183609591, .7191949522280762,
                .5360562674188973, 0}, false);

    // Exercise both u2=lambda*u1 and cm2=k*cm1. Do not assert that the
    // simple cascade DF correlation applies at changing radius.
    check_case("changing_radius", {.3, .4, .2, 1.5, 1.5, 1.2, .5, .1},
               {1.14, .5773502691896258, .6546536707079773,
                .6546536707079773, 0}, false);

    // Check the documented invalid-state policy, rather than accepting NaNs.
    expect_domain_error("gamma <= 1", compressor, 1.0);
    auto invalid = compressor;
    invalid.sigma = 0;
    expect_domain_error("zero solidity", invalid);
    invalid = compressor;
    invalid.cm1 = 2;
    expect_domain_error("nonpositive inlet static T", invalid);
    invalid = compressor;
    invalid.k = 10;
    expect_domain_error("nonpositive outlet static T", invalid);
    invalid = compressor;
    invalid.cm1 = 0;
    invalid.cth1 = invalid.u1;
    expect_domain_error("zero inlet relative speed", invalid);

    std::cout << checks - failures << '/' << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
