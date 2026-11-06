// Two-stage compressor with an interstage stator/diffuser on one spool.
// Build with the inspection-capable rostack.hpp and an NLopt-enabled pagmo:
//   g++ -std=c++20 -O3 -pthread rostack_corotating.cpp -lpagmo -o rostack_corotating
//   ./rostack_corotating
// No separate physics .cpp is required; the corrected simple model is below.
// Check the exact study configuration without pagmo or running optimization:
//   g++ -std=c++20 -O2 -DROSTACK_CHECK_ONLY rostack_corotating.cpp -o rostack_check
//   ./rostack_check
#include "rostack.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef ROSTACK_CHECK_ONLY
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/topologies/fully_connected.hpp>
#include <pagmo/topology.hpp>
#include <pagmo/utils/constrained.hpp>
#endif

namespace {
constexpr double constraint_tolerance = 1e-7; // Optimizer feasibility tolerance.

using StudyProblem = Rostack_problem<MaximizeTemperatureRatio>;
// Corrected simple adiabatic row physics, copied from rostack_physics.cpp.
// Embedded here so this study needs only rostack.hpp, not another .cpp file.
StagePerformance study_physics(const RotorVariables& rv, double gamma,
                               Rostack_diagnostics performance_params){
    if (!std::isfinite(gamma) || gamma <= 1.0 ||
        !std::isfinite(rv.sigma) || rv.sigma <= 0.0) {
        throw std::domain_error("Invalid gamma or solidity");
    }
    // All vars in rv are in the local normalization
    auto cth2 = rv.cth2;
    auto u1 = rv.u1;
    auto cm1 = rv.cm1;
    auto k = rv.k;
    auto lambda = rv.lambda;
    auto sigma = rv.sigma;
    auto cth1 = rv.cth1;

    double u2 = lambda * u1;
    double cm2 = k * cm1;

    double tau = 1. + u2 * cth2 - u1 * cth1;


    double cm1sq = cm1 * cm1;
    double cth1sq = cth1 * cth1;
    double c1sq = cm1sq + cth1sq;
    double w1sq = cm1sq + (cth1 - u1)*(cth1 - u1);

    double cm2sq = cm2 * cm2;
    double cth2sq = cth2 * cth2;
    double c2sq = cm2sq + cth2sq;
    double w2sq = cm2sq + (cth2 - u2)*(cth2 - u2);

    // Both static-temperature ratios use the row inlet TOTAL temperature
    // as their reference, and subtract ABSOLUTE flow kinetic energy.
    const double theta1 = 1.0 - c1sq / 2.0;
    const double theta2 = tau - c2sq / 2.0;
    if (!std::isfinite(tau) || tau <= 0.0 ||
        !std::isfinite(theta1) || theta1 <= 0.0 ||
        !std::isfinite(theta2) || theta2 <= 0.0 ||
        !std::isfinite(w1sq) || w1sq <= 0.0 ||
        !std::isfinite(w2sq)) {
        throw std::domain_error("Invalid stage thermodynamic state");
    }

    double Mw1sq = w1sq / ((gamma - 1.) * theta1);
    double Mw2sq = w2sq / ((gamma - 1.) * theta2);

    double w1 = std::sqrt(w1sq);
    double wth1 = cth1 - u1;
    double w2 = std::sqrt(w2sq);
    double wth2 = cth2 - u2;

    double Mw1 = std::sqrt(Mw1sq);
    double Mw2 = std::sqrt(Mw2sq);

    double M2sq = c2sq / ((gamma - 1.) * theta2);
    double M2 = std::sqrt(M2sq);

    double DF = 1. - w2/w1 + std::abs(wth1 - wth2) / (2*sigma*w1);

    // These diagnostics are populated only when inspect output is requested.
    performance_params["M1 absolute"] = std::sqrt(c1sq / ((gamma - 1.) * theta1));
    performance_params["inlet static T / row inlet total T"] = theta1;
    performance_params["outlet static T / row inlet total T"] = theta2;
    performance_params["inlet u1 (local)"] = u1;
    performance_params["outlet u2 (local)"] = u2;

    StagePerformance performance{
        tau,
        Mw1,
        Mw2,
        M2,
        DF
    };

    return performance;
}


StudyProblem make_study()
{
    // Bounds and live-flag order: cth2, u1, cm1, k, lambda, sigma, nu, cth1.
    // Each add_rotor line is independent: edit its numbers and flags directly.
    // Downstream cm1/cth1 are inherited. The second rotor's u1 is shaft-inherited.
    Rostack_model stack{study_physics, 1.4};
    Rostack_spool spool;
    spool.add_rotor({{{0.0, 1.5}, {0.0, 1.5}, {0.4, 0.4}, {1.0, 1.0}, {1.0, 1.0}, {1.5, 1.5}, {0.5, 0.5}, {0.0, 0.0}}}, {true, true, false, false, false, false, false, false}, "Stage 1 compressor");
    spool.add_rotor({{{0.0, 1.5}, {0.0, 0.0}, {0.4, 0.4}, {1.0, 1.0}, {1.0, 1.0}, {1.5, 1.5}, {0.5, 0.5}, {0.0, 0.0}}}, {true, false, false, false, false, false, false, false}, "Interstage diffuser");
    spool.add_rotor({{{0.0, 1.5}, {0.0, 1.5}, {0.4, 0.4}, {1.0, 1.0}, {1.0, 1.0}, {1.5, 1.5}, {0.5, 0.5}, {0.0, 0.0}}}, {true, true, false, false, false, false, false, false}, "Stage 2 compressor");
    stack.add_spool(spool);
    if (stack.get_n_opt_vars() != 4 || stack.stage_count() != 3) {
        throw std::logic_error("Unexpected study variable/component layout");
    }

    // Constraint row indices include the diffuser: 0 = R1, 1 = stator, 2 = R2.
    // A stationary row has identical relative and absolute Mach numbers.
    std::vector<StageConstraint> constraints;
    constraints.push_back(upper(0, &StagePerformance::diffusion, 0.55));
    constraints.push_back(upper(0, &StagePerformance::mach_rel_in, 1.2));
    constraints.push_back(upper(0, &StagePerformance::mach_rel_out, 1.2));
    constraints.push_back(upper(0, &StagePerformance::mach_abs_out, 0.8));
    constraints.push_back(upper(1, &StagePerformance::diffusion, 0.55));
    constraints.push_back(upper(1, &StagePerformance::mach_rel_in, 0.8));
    constraints.push_back(upper(1, &StagePerformance::mach_rel_out, 0.8));
    constraints.push_back(upper(1, &StagePerformance::mach_abs_out, 0.8));
    constraints.push_back(upper(2, &StagePerformance::diffusion, 0.55));
    constraints.push_back(upper(2, &StagePerformance::mach_rel_in, 1.2));
    constraints.push_back(upper(2, &StagePerformance::mach_rel_out, 1.2));
    constraints.push_back(upper(2, &StagePerformance::mach_abs_out, 0.8));

    // Absolute inlet Mach is checked at the first row (fixed inputs).
    // At every following row it equals the preceding absolute outlet Mach:
    // inlet velocities and total temperature are inherited with no intervening
    // gap heat/work. Those outlet Mach numbers are already constrained above.
    // Keep this fixed-inlet check in sync if you edit the first cm1 or gamma above.
    const double first_absolute_mach = 0.4 / std::sqrt((1.4 - 1.0) * (1.0 - 0.5 * 0.4 * 0.4));
    if (!std::isfinite(first_absolute_mach) || first_absolute_mach >= 0.8) {
        throw std::invalid_argument("Fixed first-stage absolute inlet Mach exceeds limit");
    }
    return StudyProblem{std::move(stack), std::move(constraints)};
}

void validate_incumbent(const StudyProblem& udp, const std::vector<double>& x,
                        const Rostack_inspection& report)
{
    const auto [lower, upper] = udp.get_bounds();
    if (x.size() != lower.size()) {
        throw std::invalid_argument("Incumbent dimension does not match header bounds");
    }
    for (std::size_t j = 0; j < x.size(); ++j) {
        if (!std::isfinite(x[j]) || x[j] < lower[j] || x[j] > upper[j]) {
            throw std::invalid_argument("Incumbent lies outside study bounds");
        }
    }
    if (!report.valid) throw std::domain_error("Invalid incumbent: " + report.error);
    for (std::size_t j = 1; j < report.fitness.size(); ++j) {
        if (!std::isfinite(report.fitness[j]) || report.fitness[j] > constraint_tolerance) {
            throw std::domain_error("The supplied study incumbent must be feasible");
        }
    }
}

#ifndef ROSTACK_CHECK_ONLY
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
// and shared seeding mean these are not independent statistical samples.
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
    // The extra direct-incumbent polish participates in best-result selection,
    // but is excluded from the island agreement statistics.
    std::string best_label = "  Best champion by constrained ranking";
    if (!extra.empty()) {
        best_label += best < island_count
            ? " (island " + std::to_string(best + 1) + ")" : " (supplied incumbent)";
    }
    print_fitness_diagnostics(best_label, fitness[best], p);
    print_objective_agreement(fitness, p, island_count);
}
#endif
} // namespace

int main() try
{
    const auto start = std::chrono::steady_clock::now();
    const StudyProblem udp = make_study();

    // Exact decision-vector order after spool/inlet overrides:
    // [first rotor cth2, first rotor u1, diffuser cth2, second rotor cth2].
    // These are already global-normalized velocities, not a [0,1] encoding.
    const std::vector<double> incumbent_x{0.2, 0.4, 0.1, 0.3};
    const auto incumbent_report = udp.evaluate(incumbent_x, true);
    validate_incumbent(udp, incumbent_x, incumbent_report);

    std::cout << std::setprecision(17)
              << "Two-stage compressor / interstage diffuser / one spool\n"
              << "Independent variables: " << incumbent_x.size()
              << ", equality constraints: " << udp.get_nec()
              << ", inequality constraints: " << udp.get_nic() << '\n'
              << "Objective: minimize -tau_total (maximize tau_total)\n"
              << "Component bounds, fixed inputs and limits: see explicit lines in make_study().\n"
              << "Constraint residuals and resolved component inputs are printed below."
              << " (constraint tolerance " << constraint_tolerance << ")\n"
              << "Supplied incumbent tau: " << incumbent_report.result.tau_total << '\n';

#ifdef ROSTACK_CHECK_ONLY
    // Exercise the actual study setup and physics without pagmo dependencies.
    if (udp.get_nec() != 0 || udp.get_nic() != 12 ||
        std::abs(incumbent_report.result.tau_total - 1.16) > 1e-12 ||
        std::abs(incumbent_report.design.at("Stage 2 compressor, inlet u1") - 0.4) > 1e-12 ||
        incumbent_report.design.at("Interstage diffuser, inlet u1") != 0.0) {
        throw std::logic_error("Study configuration check failed");
    }
    std::cout << "CHECK ONLY: configuration and feasible seed verified; no optimization run.\n";
    print_inspection(std::cout, incumbent_report);
    std::cout << "Execution time: "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
              << " s\n";
    return 0;
#else
    pagmo::problem problem{udp};
    problem.set_c_tol(constraint_tolerance);
    std::cout << problem;
    const auto incumbent_f = problem.fitness(incumbent_x);
    print_fitness_diagnostics("Supplied incumbent", incumbent_f, problem);

    // GACO, island population sizes, migration topology, number of rounds and
    // COBYLA polishing settings below match the supplied engine_opt driver.
  pagmo::algorithm algo{pagmo::gaco(
        100,//gen
        31,//ker
        1,//q
        100,//oracle
        0.01,//acc
        100,//convergence speed threshold
        7,//n gen mark
        100,//impstop
        100000,//evalstop
        0,//focus
        true //memory
    )};
  pagmo::fully_connected fc_udt{};
  pagmo::topology topo{fc_udt};
  pagmo::archipelago archi(topo, 32u, algo, problem, 6000u);
  // Replace one individual on every island; keep the other random candidates.
  for (auto& island : archi) {
      auto pop = island.get_population();
      pop.set_xf(pop.worst_idx(), incumbent_x, incumbent_f);
      island.set_population(pop);
  }
  const int n_evolves = 8; // Same number of rounds as the supplied driver.
  for (int evolve = 0; evolve < n_evolves; evolve++) {
    archi.evolve(1);
    archi.wait_check();
    const auto champions_f = archi.get_champions_f();
    print_progress_summary("Round " + std::to_string(evolve + 1) + " / "
                           + std::to_string(n_evolves), champions_f, problem, champions_f.size());
  }


    // ----- Local refinement of island champions -----
    pagmo::algorithm local_algo{pagmo::nlopt("cobyla")};
    auto* nl = local_algo.extract<pagmo::nlopt>();
    nl->set_maxeval(2000);
    nl->set_xtol_rel(1e-10);
    nl->set_ftol_rel(1e-10);
    nl->set_xtol_abs(1e-12);
    local_algo.set_verbosity(0);

    auto champions_x = archi.get_champions_x();
    // Also polish the supplied incumbent directly, even if global search lost it.
    const std::size_t incumbent_candidate = champions_x.size();
    champions_x.push_back(incumbent_x);

    // Advance every candidate by one local pass before reporting the best.
    // This preserves each candidate's original stopping rule and solver budget.
    std::vector<pagmo::population> local_populations;
    local_populations.reserve(champions_x.size());
    std::vector<pagmo::vector_double> refined_x = champions_x;
    std::vector<pagmo::vector_double> refined_f;
    refined_f.reserve(champions_x.size());
    for (const auto& x : champions_x) {
        local_populations.emplace_back(problem, 0);
        local_populations.back().push_back(x);
        refined_f.push_back(local_populations.back().champion_f());
    }
    std::vector<bool> active(champions_x.size(), true);
    constexpr int max_passes = 20;
    for (int pass = 0; pass < max_passes; ++pass) {
        std::size_t attempted = 0;
        std::size_t errors = 0;
        std::string first_error;
        for (std::size_t i = 0; i < local_populations.size(); ++i) {
            if (!active[i]) continue;
            ++attempted;
            auto& pop = local_populations[i];
            const double previous_objective = pop.champion_f()[0];
            try {
                pop = local_algo.evolve(pop);
            } catch (const std::exception& e) {
                ++errors;
                active[i] = false;
                if (first_error.empty()) {
                    const std::string source = i < incumbent_candidate
                        ? "island " + std::to_string(i + 1) : "supplied incumbent";
                    first_error = source + ": " + e.what();
                    std::replace(first_error.begin(), first_error.end(), '\n', ' ');
                    std::replace(first_error.begin(), first_error.end(), '\r', ' ');
                }
                continue; // Retain the previous candidate if polishing fails.
            }
            refined_x[i] = pop.champion_x();
            refined_f[i] = pop.champion_f();
            if (std::abs(previous_objective - refined_f[i][0])
                < 1e-8 * std::abs(previous_objective)) {
                active[i] = false;
            }
        }
        if (attempted == 0) break;
        const std::size_t remaining = std::count(active.begin(), active.end(), true);
        std::string extra = ", polished=" + std::to_string(attempted)
                          + ", remaining active=" + std::to_string(remaining)
                          + ", errors=" + std::to_string(errors);
        if (errors > 0) extra += ", first error=" + first_error;
        print_progress_summary("Local polish pass " + std::to_string(pass + 1) + " / "
                               + std::to_string(max_passes), refined_f, problem,
                               incumbent_candidate, extra);
        if (remaining == 0) break;
    }

    // ----- Rank refined champions and report the best -----
    const auto order = pagmo::sort_population_con(
        refined_f, problem.get_nec(), problem.get_c_tol());
    const std::size_t i = order.front();
    const auto& x = refined_x[i];
    const auto& f = refined_f[i];
    const auto report = udp.evaluate(x, true);

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
              << "No. refined champions: " << refined_f.size() << '\n'
              << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
              << "Rank 1:\n"
              << "Source: " << (i == incumbent_candidate ? "Supplied incumbent" : "Island champion") << '\n'
              << "Score (-tau): " << f[0] << '\n';
    if (report.valid) std::cout << "Total temperature ratio: " << report.result.tau_total << '\n';
    print_fitness_diagnostics("Final candidate", f, problem);
    if (!problem.feasibility_f(f)) {
        std::cout << "No feasible optimum found; showing the best constrained-ranked candidate.\n";
    }
    if (report.fitness != f) {
        throw std::logic_error("Inspection fitness differs from the stored candidate fitness");
    }

    constexpr std::array<const char*, 4> variable_names{
        "Stage 1 compressor, outlet cth2", "Stage 1 compressor, inlet u1",
        "Interstage diffuser, outlet cth2", "Stage 2 compressor, outlet cth2"
    };
    std::cout << "Decision vector (global-normalized velocities):\n{\n";
    for (std::size_t j = 0; j < x.size(); ++j) {
        std::cout << "    " << std::quoted(variable_names[j]) << ": " << x[j]
                  << (j + 1 == x.size() ? "\n" : ",\n");
    }
    std::cout << "}\n";
    print_inspection(std::cout, report);

    const auto end = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "Execution time: " << std::fixed << std::setprecision(3)
              << elapsed << " s\n";
    return 0;
#endif
}
catch (const std::exception& e)
{
    std::cerr << "rostack_corotating: " << e.what() << '\n';
    return 1;
}
