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
#include <pagmo/topology.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <numeric>
#include <algorithm>

#include "jet_calc.h"

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
  std::cout << pjc;
  algorithm algo{gaco(100, 63, 1, 0, 0.01, 100, 7, 100)};
  pagmo::ring ring_udt{};
  pagmo::topology topo{ring_udt};
  archipelago archi(topo, 32u, algo, pjc, 6000u);
  // std::cout << archi;
  int n_evolves = 3; // Increase this to 10 to get the last 0.1%
  for (int evolve = 0; evolve < n_evolves; evolve++) {
    archi.evolve(1);
    archi.wait_check();
    const auto& champions_f = archi.get_champions_f();
    double max_fitness = champions_f[0][0];
    for (const auto& f : champions_f) {
      if (f[0] > max_fitness) {
        max_fitness = f[0];
      }
    }

    std::cout << "Round " << evolve + 1 << " / " << n_evolves
              << " — max champion fitness: " << max_fitness << '\n';
  }


    // ----- Local refinement of island champions -----
    algorithm local_algo{nlopt("cobyla")};
    auto* nl = local_algo.extract<nlopt>();
    nl->set_maxeval(2000);       // much bigger budget than 200
    nl->set_xtol_rel(1e-10);
    nl->set_ftol_rel(1e-10);
    nl->set_xtol_abs(1e-12);
    local_algo.set_verbosity(0);
    local_algo.extract<nlopt>()->set_maxeval(200); // tune budget as needed

    auto champions_x = archi.get_champions_x();
    auto champions_f = archi.get_champions_f();

    std::vector<pagmo::vector_double> refined_x;
    std::vector<pagmo::vector_double> refined_f;

    for (std::size_t i = 0; i < champions_x.size(); ++i) {
    // seed a 1-individual population at the champion and polish locally
    population pop(pjc, 0);
    pop.push_back(champions_x[i]);
    pop = local_algo.evolve(pop);
    refined_x.push_back(pop.champion_x());
    refined_f.push_back(pop.champion_f());
    }

    // ----- Rank refined champions and print the top k -----
    int k = 1; // Number of champions to report
    k = std::min<int>(k, static_cast<int>(refined_f.size()));

    std::vector<std::size_t> idx(refined_f.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
    return refined_f[a][0] < refined_f[b][0]; // descending, higher = better
    });

    auto nobj = pjc.get_nobj();
    auto nec  = pjc.get_nec();
    auto nic  = pjc.get_nic();
    auto ncon = nec + nic;

    for (int rank = 0; rank < k; ++rank) {
    std::size_t i = idx[rank];
    const auto& x = refined_x[i];
    const auto& f = refined_f[i];

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "Rank " << rank + 1 << ":\n";
    std::cout << "Score: " << f[0] << "\n";

    std::cout << "Design vector: \n{\n";
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