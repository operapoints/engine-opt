#include <chrono>
#include <iomanip>
#include <iostream>
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/algorithms/sade.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/topologies/ring.hpp>
#include <pagmo/topologies/fully_connected.hpp>
#include <pagmo/topology.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <pagmo/utils/constrained.hpp>

#include "jet_calc.h"

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
} // namespace

int main() {
  // auto omega = x[0]; // rad/s - shaft speed
  // auto u_i = x[1];   // m/s   - compressor inlet velocity
  // auto T_4 = x[2];   // K     - combustor exit total temp
  // auto R_Cih = x[3]; // m     - compressor inlet hub radius
  // auto R_Cit = x[4]; // m     - compressor inlet tip radius
  // auto A_Co = x[5];  // m^2   - compressor outlet area
  // auto R_Com = x[6]; // m     - compressor outlet meanline radius
  // auto D_T_C = x[7]; // K     - compressor total temperature change
  // auto R_Tih = x[8]; // m     - turbine inlet hub radius
  // auto R_Tit = x[9]; // m     - Turbine inlet tip radius
  // auto A_To = x[10]; // m^2   - Turbine outlet area
  // auto R_Tom = x[11];// m     - Turbine exit meanline velocity
  problem_jet_calc pjc_obj;
  auto start = std::chrono::steady_clock::now();
#ifdef EVAL_JET_CALC
  // 10cm
  pagmo::vector_double x1 = {
      4648.1,    136.75,  1099.75,   0.0150023, 0.0507001,  0.00750535,
      0.0799953, 73.5045, 0.0511641, 0.0774573, 0.00759833, 0.0576068,
  };
  // 8cm
  pagmo::vector_double x2 = {
      4648.1,    136.75,  1099.75,   0.0150023, 0.0507001,  0.00750535,
      0.0799953, 72.5045, 0.0511641, 0.0774573, 0.00759833, 0.0576068,
  };
  // int nsteps = 10;
  // double min_OPR = 1.5;
  // double max_OPR = 2.2;
  // for (int i = 0; i < nsteps; i++){
  //     auto xi = x2;
  //     double OPRi = min_OPR + (max_OPR - min_OPR)*(float(i)/float(nsteps));
  //     //xi[7] = 298*(std::pow(OPRi,(0.4/1.4))-1)/0.75;
  //     auto res = pjc_obj.fitness(xi);
  //     std::cout << "D T C : " << OPRi << " || Objective : " << -res[0] <<
  //     "\n";
  // }

  // Max thrust
  pagmo::vector_double x3 = {
      14684.1,   125.722, 991.887,   0.00504372, 0.0164868,  0.000539293,
      0.0299262, 104.827, 0.0155055, 0.0245546,  0.00169747, 0.0209158,
  };
  pjc_obj.fitness(x1);
  pjc_obj.fitness(x2);
  pjc_obj.fitness(x3);
#endif
  pagmo::problem pjc{pjc_obj};
  pjc.set_c_tol(1e-7);
  std::cout << std::setprecision(17) << pjc;
  const pagmo::vector_double incumbent_physical = {
      7363.24, 253.786, 1097.32, 0.0230949, 0.0418834, 0.0114587,
      0.0799527, 183.448, 0.049076, 0.0613163, 0.00424502, 0.055194
  };
  if (incumbent_physical.size() != problem_jet_calc::design_size) {
      throw std::invalid_argument("Incumbent dimension does not match header bounds");
  }
  auto incumbent_x = incumbent_physical;
  for (std::size_t j = 0; j < incumbent_x.size(); ++j) {
      const double span = problem_jet_calc::physical_upper[j]
                        - problem_jet_calc::physical_lower[j];
      if (!(span > 0.0)) throw std::invalid_argument("Invalid physical bounds");
      incumbent_x[j] = (incumbent_x[j] - problem_jet_calc::physical_lower[j]) / span;
      if (!std::isfinite(incumbent_x[j]) || incumbent_x[j] < 0.0 || incumbent_x[j] > 1.0) {
          throw std::invalid_argument("Incumbent lies outside header bounds");
      }
  }
  const auto incumbent_f = pjc.fitness(incumbent_x);
  algorithm algo{gaco(
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
  archipelago archi(topo, 32u, algo, pjc, 6000u);
  // Replace one individual on every island; keep the other random candidates.
  for (auto& island : archi) {
      auto pop = island.get_population();
      pop.set_xf(pop.worst_idx(), incumbent_x, incumbent_f);
      island.set_population(pop);
  }
  int n_evolves = 8; // Increase this to 10 to get the last 0.1%
  for (int evolve = 0; evolve < n_evolves; evolve++) {
    archi.evolve(1);
    archi.wait_check();
    const auto champions_f = archi.get_champions_f();
    print_progress_summary("Round " + std::to_string(evolve + 1) + " / "
                           + std::to_string(n_evolves), champions_f, pjc, champions_f.size());
  }


    // ----- Local refinement of island champions -----
    algorithm local_algo{nlopt("cobyla")};
    auto* nl = local_algo.extract<nlopt>();
    nl->set_maxeval(2000);       // much bigger budget than 200
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
    std::vector<population> local_populations;
    local_populations.reserve(champions_x.size());
    std::vector<pagmo::vector_double> refined_x = champions_x;
    std::vector<pagmo::vector_double> refined_f;
    refined_f.reserve(champions_x.size());
    for (const auto& x : champions_x) {
        local_populations.emplace_back(pjc, 0);
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
                               + std::to_string(max_passes), refined_f, pjc,
                               incumbent_candidate, extra);
        if (remaining == 0) break;
    }

    // Keep refined_x normalized for further fitness() evaluations.
    // Store a separate copy in SI units for printing and physical calculations.
    auto physical_x = refined_x;
    for (auto& x : physical_x) {
        for (std::size_t j = 0; j < problem_jet_calc::design_size; ++j) {
            x[j] = problem_jet_calc::physical_lower[j]
                + x[j] * (problem_jet_calc::physical_upper[j]
                          - problem_jet_calc::physical_lower[j]);
        }
    }

    // ----- Rank refined champions and print the top k -----
    int k = 1; // Number of champions to report
    k = std::min<int>(k, static_cast<int>(refined_f.size()));

    // Feasible candidates first; among feasible candidates, smaller objective wins.
    const auto idx = pagmo::sort_population_con(refined_f, pjc.get_nec(), pjc.get_c_tol());
    auto nobj = pjc.get_nobj();
    auto nec  = pjc.get_nec();
    auto nic  = pjc.get_nic();
    auto ncon = nec + nic;

    for (int rank = 0; rank < k; ++rank) {
    std::size_t i = idx[rank];
    const auto& x = physical_x[i];
    const auto& f = refined_f[i];
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "No. refined champions: " << refined_f.size() <<'\n';
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "Rank " << rank + 1 << ":\n";
    std::cout << "Source: " << (i == incumbent_candidate ? "Supplied incumbent" : "Island champion") << '\n';
    std::cout << "Score: " << f[0] << "\n";

    std::cout << "Design vector (SI units): \n{\n";
    for (std::size_t j = 0; j < x.size(); ++j) {
        std::cout << "  "<< j+1<< ": " << x[j];
        if (j + 1 < x.size()) std::cout << ",";
        std::cout << "\n";
    }
    std::cout << "}\n";

    if (ncon > 0) {
        std::cout << "Constraints:\n{\n";
        for (std::size_t j = 0; j < ncon; ++j) {
        std::cout <<  "  "<< j+1<< ": " << f[nobj + j];
        if (j + 1 < ncon) std::cout << ",";
        std::cout << "\n";
        }
        std::cout << "}\n";
    }
    }

    auto end = std::chrono::steady_clock::now();

    double elapsed = std::chrono::duration<double>(end - start).count();

    std::cout << "Execution time: " << std::fixed << std::setprecision(3)
            << elapsed << " s\n";
  // pagmo::vector_double x0 =
  // {10183,127,1100,0.003,0.023,0.001,0.0250,65,0.02,0.03,0.002,0.025};
  // pagmo::vector_double ret = pjc.fitness(x0);

  // std::cout << "ret: " << '\n';
  // for (int i = 0 ; i < ret.size(); i++){
  //     std::cout << ret[i] << '\n';
  // }
  return 0;
}
