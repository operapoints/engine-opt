#pragma once
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>

struct OptimizationResult {
    pagmo::vector_double x;
    pagmo::vector_double fitness;
    bool feasible = false;
    double elapsed_seconds = 0.0;
};

// Single-objective optimization with optional equality/inequality constraints.
// Uses the existing GACO archipelago + COBYLA settings. Progress goes to stdout.
// Accepts a pagmo::problem directly; its constraint tolerance is set to 1e-7.
OptimizationResult optimize_udp(pagmo::problem problem);

// Accept any pagmo-compatible UDP, including Rostack_problem<YourObjective>.
// The non-template overload above handles the optimization implementation.
template<class UDP>
OptimizationResult optimize_udp(const UDP& udp)
{
    return optimize_udp(pagmo::problem{udp});
}