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

const char* nlopt_status_name(int status) {
    switch (status) {
        case -1: return "FAILURE";
        case -2: return "INVALID_ARGS";
        case -3: return "OUT_OF_MEMORY";
        case -4: return "ROUNDOFF_LIMITED";
        case -5: return "FORCED_STOP";
        case 1: return "SUCCESS";
        case 2: return "STOPVAL_REACHED";
        case 3: return "FTOL_REACHED";
        case 4: return "XTOL_REACHED";
        case 5: return "MAXEVAL_REACHED";
        case 6: return "MAXTIME_REACHED";
        default: return "UNKNOWN";
    }
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
  print_fitness_diagnostics("Supplied incumbent (reevaluated)", incumbent_f, pjc);
  std::cout << "Incumbent normalized vector:\n";
  for (std::size_t j = 0; j < incumbent_x.size(); ++j) {
      std::cout << "  " << j + 1 << ": " << incumbent_x[j] << '\n';
  }
  std::cout << "Incumbent constraint residuals and configured tolerances:\n";
  for (std::size_t j = 0; j < pjc.get_nc(); ++j) {
      std::cout << "  " << j + 1 << ": " << incumbent_f[pjc.get_nobj() + j]
                << ", tolerance=" << pjc.get_c_tol()[j] << '\n';
  }
  if (!pjc.feasibility_f(incumbent_f)) {
      std::cout << "NOTE: This rounded saved vector is not feasible when reevaluated. "
                   "No constraint tolerances have been relaxed.\n";
  }
  const double incumbent_r_cot = incumbent_physical[6]
      + incumbent_physical[5] / (4.0 * std::acos(-1.0) * incumbent_physical[6]);
  std::cout << "Incumbent geometric R_Cot [m]=" << incumbent_r_cot << '\n';

  // Small feasible descent checks distinguish a true boundary from solver stalls.
  const auto probe_descent = [&](const std::string& label, const pagmo::vector_double& x) {
      const auto base_f = pjc.fitness(x);
      for (const std::size_t j : {std::size_t(5), std::size_t(6)}) {
          for (const double requested_h : {1e-4, 1e-6, 1e-8}) {
              if (x[j] <= 0.0) continue;
              auto trial_x = x;
              const double h = std::min(requested_h, x[j]);
              trial_x[j] -= h;
              const auto trial_f = pjc.fitness(trial_x);
              std::cout << label << " probe: decrease " << (j == 5 ? "A_Co" : "R_Com")
                        << ", normalized step=" << h
                        << ", objective change=" << trial_f[0] - base_f[0] << '\n';
              print_fitness_diagnostics("  trial", trial_f, pjc);
          }
      }
  };
  probe_descent("Incumbent", incumbent_x);
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
  std::cout << "Seeded all 32 islands with the normalized incumbent.\n";
  int n_evolves = 8; // Increase this to 10 to get the last 0.1%
  for (int evolve = 0; evolve < n_evolves; evolve++) {
    archi.evolve(1);
    archi.wait_check();
    const auto champions_f = archi.get_champions_f();
    const auto order = pagmo::sort_population_con(champions_f, pjc.get_nec(), pjc.get_c_tol());
    std::size_t feasible_count = 0;
    double best_feasible = std::numeric_limits<double>::infinity();
    double largest_objective = -std::numeric_limits<double>::infinity();
    for (const auto& f : champions_f) {
        largest_objective = std::max(largest_objective, f[0]);
        if (pjc.feasibility_f(f)) {
            ++feasible_count;
            best_feasible = std::min(best_feasible, f[0]);
        }
    }
    std::cout << "Round " << evolve + 1 << " / " << n_evolves
              << " — feasible champions: " << feasible_count << '/' << champions_f.size()
              << ", best feasible objective: " << best_feasible
              << ", largest champion objective: " << largest_objective << '\n';
    print_fitness_diagnostics("  Best champion by constrained ranking", champions_f[order[0]], pjc);
  }


    // ----- Local refinement of island champions -----
    algorithm local_algo{nlopt("cobyla")};
    auto* nl = local_algo.extract<nlopt>();
    nl->set_maxeval(2000);       // much bigger budget than 200
    nl->set_xtol_rel(1e-10);
    nl->set_ftol_rel(1e-10);
    nl->set_xtol_abs(1e-12);
    local_algo.set_verbosity(100);

    auto champions_x = archi.get_champions_x();
    // Also polish the supplied incumbent directly, even if global search lost it.
    const std::size_t incumbent_candidate = champions_x.size();
    champions_x.push_back(incumbent_x);

    std::vector<pagmo::vector_double> refined_x;
    std::vector<pagmo::vector_double> refined_f;

    std::cout << "Local settings: maxeval=2000 per pass, max passes=20, "
                 "xtol_rel=1e-10, ftol_rel=1e-10, xtol_abs=1e-12\n";
    for (std::size_t i = 0; i < champions_x.size(); ++i) {
        const std::string label = i == incumbent_candidate
            ? "Supplied incumbent" : "Island " + std::to_string(i + 1);
        population pop(pjc, 0);
        pop.push_back(champions_x[i]);
        print_fitness_diagnostics(label + " before polish", pop.champion_f(), pjc);
        for (int pass = 0; pass < 20; ++pass) {
            const auto previous_x = pop.champion_x();
            const auto previous_f = pop.champion_f();
            const auto fevals_before = pop.get_problem().get_fevals();
            try {
                pop = local_algo.evolve(pop);
            } catch (const std::exception& e) {
                std::cout << label << " pass " << pass + 1
                          << " threw: " << e.what() << "; retaining previous candidate\n";
                break;
            }
            const auto current_f = pop.champion_f();
            double dx_inf = 0.0;
            for (std::size_t j = 0; j < previous_x.size(); ++j) {
                dx_inf = std::max(dx_inf, std::abs(pop.champion_x()[j] - previous_x[j]));
            }
            // Re-extract from the algorithm that actually ran.
            const int status = static_cast<int>(local_algo.extract<nlopt>()->get_last_opt_result());
            const auto fevals = pop.get_problem().get_fevals() - fevals_before;
            std::cout << label << " pass " << pass + 1
                      << ": " << nlopt_status_name(status) << " (" << status << ')'
                      << ", fitness evaluations=" << fevals
                      << ", objective before=" << previous_f[0]
                      << ", after=" << current_f[0]
                      << ", delta=" << current_f[0] - previous_f[0]
                      << ", max normalized coordinate change=" << dx_inf << '\n';
            print_fitness_diagnostics("  Retained candidate", current_f, pjc);
            // Keep the existing outer stopping rule, but make its activation explicit.
            if (std::abs(previous_f[0] - current_f[0]) < 1e-8 * std::abs(previous_f[0])) {
                std::cout << "  Outer objective-change stop triggered; "
                             "this does not certify local optimality.\n";
                break;
            }
        }
        refined_x.push_back(pop.champion_x());
        refined_f.push_back(pop.champion_f());
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
    print_fitness_diagnostics("Best retained result", refined_f[idx[0]], pjc);
    std::cout << "Best retained objective minus supplied incumbent objective: "
              << refined_f[idx[0]][0] - incumbent_f[0] << '\n';
    probe_descent("Best retained result", refined_x[idx[0]]);

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
    print_fitness_diagnostics("Result", f, pjc);
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
