#include "rostack_optimizer.hpp"

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/population.hpp>
#include <pagmo/topologies/fully_connected.hpp>
#include <pagmo/topology.hpp>
#include <pagmo/utils/constrained.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {
struct ConstraintDiagnostics {
    std::size_t violated = 0;
    double max_violation = 0.0; // Excess beyond the configured tolerance.
    double max_inequality = -std::numeric_limits<double>::infinity();
};

ConstraintDiagnostics constraint_diagnostics(const pagmo::vector_double& f,
                                             const pagmo::problem& p) {
    ConstraintDiagnostics d;
    const auto& tol = p.get_c_tol();
    for (std::size_t j = 0; j < tol.size(); ++j) {
        const double c = f[p.get_nobj() + j];
        if (j >= p.get_nec()) d.max_inequality = std::max(d.max_inequality, c);
        const double residual = j < p.get_nec() ? std::abs(c) : c;
        const double violation = std::isfinite(c)
            ? std::max(0.0, residual - tol[j])
            : std::numeric_limits<double>::infinity();
        if (violation > 0.0) ++d.violated;
        d.max_violation = std::max(d.max_violation, violation);
    }
    return d;
}

void print_fitness_diagnostics(const std::string& label,
                               const pagmo::vector_double& f,
                               const pagmo::problem& p) {
    const auto d = constraint_diagnostics(f, p);
    std::cout << label << ": objective=" << f[0]
              << ", feasible=" << (p.feasibility_f(f) ? "yes" : "no")
              << ", violated=" << d.violated
              << ", max violation beyond tolerance=" << d.max_violation
              << ", largest inequality residual=" << d.max_inequality << '\n';
}

// Agreement is measured in objective space among feasible island champions only.
// It is a convergence proxy, not a certificate of global optimality; migration
// means these are not independent statistical samples.
struct ObjectiveAgreement {
    std::size_t count = 0;
    std::size_t near_best = 0;
    std::size_t top_count = 0;
    double median = 0.0;
    double iqr = 0.0;
    double top_spread_percent = 0.0;
};

ObjectiveAgreement objective_agreement(std::vector<double> values) {
    ObjectiveAgreement result;
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double v) { return !std::isfinite(v); }), values.end());
    if (values.empty()) return result;
    std::sort(values.begin(), values.end());
    result.count = values.size();
    const auto quantile = [&](double fraction) {
        const double position = fraction * static_cast<double>(values.size() - 1);
        const auto lo = static_cast<std::size_t>(position);
        const auto hi = std::min(lo + 1, values.size() - 1);
        return values[lo] + (position - lo) * (values[hi] - values[lo]);
    };
    result.median = quantile(0.5);
    result.iqr = quantile(0.75) - quantile(0.25);
    // Relative gaps also work for negative objectives. Near zero, use a small
    // absolute scale floor to avoid division by zero.
    const double scale = std::max(std::abs(values.front()), 1e-12);
    const double near_best_gap = 0.001 * scale; // Within 0.1% of the best.
    for (double value : values) {
        if (value - values.front() <= near_best_gap) ++result.near_best;
    }
    if (values.size() >= 4) {
        result.top_count = (values.size() + 3) / 4; // Best quarter, rounded up.
        result.top_spread_percent = 100.0
            * (values[result.top_count - 1] - values.front()) / scale;
    }
    return result;
}

void print_objective_agreement(const std::vector<pagmo::vector_double>& fitness,
                               const pagmo::problem& p, std::size_t island_count) {
    std::vector<double> feasible_objectives;
    for (std::size_t i = 0; i < std::min(island_count, fitness.size()); ++i) {
        if (p.feasibility_f(fitness[i])) feasible_objectives.push_back(fitness[i][0]);
    }
    const auto a = objective_agreement(std::move(feasible_objectives));
    std::cout << "  Objective agreement (feasible islands): ";
    if (a.count == 0) {
        std::cout << "n/a (none feasible)\n";
        return;
    }
    std::cout << "median=" << a.median << ", IQR=" << a.iqr;
    if (a.top_count > 0) {
        std::cout << ", best-quarter spread=" << a.top_spread_percent
                  << "% (n=" << a.top_count << ')';
    } else {
        std::cout << ", best-quarter spread=n/a (fewer than 4 feasible islands)";
    }
    std::cout << ", within 0.1% of best=" << a.near_best << '/' << a.count << '\n';
}

void print_progress_summary(const std::string& label,
                            const std::vector<pagmo::vector_double>& fitness,
                            const pagmo::problem& p, std::size_t island_count,
                            const std::string& extra = "") {
    if (fitness.empty()) return;
    const auto order = pagmo::sort_population_con(fitness, p.get_nec(), p.get_c_tol());
    const auto best = order[0];
    std::size_t feasible_count = 0;
    double best_feasible = std::numeric_limits<double>::infinity();
    double largest_objective = -std::numeric_limits<double>::infinity();
    for (const auto& f : fitness) {
        largest_objective = std::max(largest_objective, f[0]);
        if (p.feasibility_f(f)) {
            ++feasible_count;
            best_feasible = std::min(best_feasible, f[0]);
        }
    }
    std::cout << label << " — feasible champions: " << feasible_count << '/' << fitness.size()
              << ", best feasible objective: ";
    if (feasible_count > 0) std::cout << best_feasible;
    else std::cout << "n/a";
    std::cout << ", largest champion objective: " << largest_objective << extra << '\n';
    const std::string best_label = "  Best champion by constrained ranking (island "
        + std::to_string(best + 1) + ")";
    print_fitness_diagnostics(best_label, fitness[best], p);
    print_objective_agreement(fitness, p, island_count);
}
} // namespace

OptimizationResult optimize_udp(pagmo::problem problem)
{
    const auto start = std::chrono::steady_clock::now();
    problem.set_c_tol(1e-7);
    std::cout << std::setprecision(17) << problem;

    pagmo::algorithm global_algorithm{pagmo::gaco(
        100,     // generations per evolution
        31,      // archive size
        1.0,     // q
        100.0,   // oracle
        0.01,    // accuracy
        100,     // convergence-speed threshold
        7,       // generation mark
        100,     // improvement stop
        100000,  // evaluation stop
        0.0,     // focus
        true     // memory
    )};
    pagmo::topology topology{pagmo::fully_connected{}};
    pagmo::archipelago archipelago{topology, 32u, global_algorithm, problem, 6000u};

    // Populations are initialized by pagmo; no supplied incumbent or seed.
    for (int round = 0; round < 8; ++round) {
        archipelago.evolve(1);
        archipelago.wait_check();
        const auto fitness = archipelago.get_champions_f();
        print_progress_summary("Round " + std::to_string(round + 1) + " / 8",
                               fitness, problem, fitness.size());
    }

    pagmo::algorithm local_algorithm{pagmo::nlopt{"cobyla"}};
    auto* nlopt = local_algorithm.extract<pagmo::nlopt>();
    nlopt->set_maxeval(2000);
    nlopt->set_xtol_rel(1e-10);
    nlopt->set_ftol_rel(1e-10);
    nlopt->set_xtol_abs(1e-12);
    local_algorithm.set_verbosity(0);

    auto refined_x = archipelago.get_champions_x();
    std::vector<pagmo::vector_double> refined_f;
    refined_f.reserve(refined_x.size());
    std::vector<pagmo::population> local_populations;
    local_populations.reserve(refined_x.size());
    for (const auto& x : refined_x) {
        local_populations.emplace_back(problem, 0u);
        local_populations.back().push_back(x);
        refined_f.push_back(local_populations.back().champion_f());
    }

    std::vector<bool> active(refined_x.size(), true);
    for (int pass = 0; pass < 20; ++pass) {
        std::size_t attempted = 0;
        std::size_t errors = 0;
        std::string first_error;
        for (std::size_t i = 0; i < local_populations.size(); ++i) {
            if (!active[i]) continue;
            ++attempted;
            auto& population = local_populations[i];
            const double previous_objective = population.champion_f()[0];
            try {
                population = local_algorithm.evolve(population);
            } catch (const std::exception& e) {
                ++errors;
                active[i] = false;
                if (first_error.empty()) {
                    first_error = "island " + std::to_string(i + 1) + ": " + e.what();
                    std::replace(first_error.begin(), first_error.end(), '\n', ' ');
                    std::replace(first_error.begin(), first_error.end(), '\r', ' ');
                }
                continue; // Preserve this island's previous result.
            }
            refined_x[i] = population.champion_x();
            refined_f[i] = population.champion_f();
            if (std::abs(previous_objective - refined_f[i][0])
                < 1e-8 * std::abs(previous_objective)) {
                active[i] = false;
            }
        }
        if (attempted == 0) break;
        const auto remaining = std::count(active.begin(), active.end(), true);
        std::string extra = ", polished=" + std::to_string(attempted)
            + ", remaining active=" + std::to_string(remaining)
            + ", errors=" + std::to_string(errors);
        if (errors > 0) extra += ", first error=" + first_error;
        print_progress_summary("Local polish pass " + std::to_string(pass + 1) + " / 20",
                               refined_f, problem, refined_x.size(), extra);
        if (remaining == 0) break;
    }

    const auto ranking = pagmo::sort_population_con(
        refined_f, problem.get_nec(), problem.get_c_tol());
    const auto best = ranking.front();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    OptimizationResult result{refined_x[best], refined_f[best],
                              problem.feasibility_f(refined_f[best]), elapsed};

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
              << "No. refined champions: " << refined_f.size() << '\n'
              << "Rank 1: island " << best + 1 << '\n';
    print_fitness_diagnostics("Final candidate", result.fitness, problem);
    if (!result.feasible) {
        std::cout << "No feasible optimum found; returning the best constrained-ranked candidate.\n";
    }
    std::ostringstream duration;
    duration << std::fixed << std::setprecision(3) << elapsed;
    std::cout << "Optimization time: " << duration.str() << " s\n";
    return result;
}