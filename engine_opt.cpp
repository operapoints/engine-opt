#include <iostream>
#include <chrono>

#include <pagmo/algorithm.hpp>

#include <pagmo/algorithms/sade.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/population.hpp>

#include "jet_calc.h"

int main(){
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
    #ifdef EVAL_JET_CALC
    // 10cm
    pagmo::vector_double x1 = {4648.1, 136.75, 1099.75, 0.0150023, 0.0507001, 0.00750535, 0.0799953, 73.5045, 0.0511641, 0.0774573, 0.00759833, 0.0576068, };
    // 8cm
    pagmo::vector_double x2 = {4648.1, 136.75, 1099.75, 0.0150023, 0.0507001, 0.00750535, 0.0799953, 72.5045, 0.0511641, 0.0774573, 0.00759833, 0.0576068, };
    // int nsteps = 10;
    // double min_OPR = 1.5;
    // double max_OPR = 2.2;
    // for (int i = 0; i < nsteps; i++){
    //     auto xi = x2;
    //     double OPRi = min_OPR + (max_OPR - min_OPR)*(float(i)/float(nsteps));
    //     //xi[7] = 298*(std::pow(OPRi,(0.4/1.4))-1)/0.75;
    //     auto res = pjc_obj.fitness(xi);
    //     std::cout << "D T C : " << OPRi << " || Objective : " << -res[0] << "\n";
    // }

    // Max thrust
    pagmo::vector_double x3 = {14684.1, 125.722, 991.887, 0.00504372, 0.0164868, 0.000539293, 0.0299262, 104.827, 0.0155055, 0.0245546, 0.00169747, 0.0209158, };
    pjc_obj.fitness(x1);
    pjc_obj.fitness(x2);
    pjc_obj.fitness(x3);
    #endif
    pagmo::problem pjc{pjc_obj};
    std::cout << pjc;
    algorithm algo{gaco(10000,63,1,0,0.01,100,7,10000)};    
    archipelago archi(32u, algo, pjc, 6000u);
    archi.evolve(1);
    archi.wait_check();

    for (const auto &isl : archi) {
        std::cout << isl.get_population().champion_f()[0] << '\n';
        pagmo::vector_double ch_x = isl.get_population().champion_x();
        std::cout << "{";
        for(int i =0; i<ch_x.size(); i++){
            std::cout << ch_x[i] << ", ";
        }
        std::cout << "}\n";
    }


    // pagmo::vector_double x0 = {10183,127,1100,0.003,0.023,0.001,0.0250,65,0.02,0.03,0.002,0.025};
    // pagmo::vector_double ret = pjc.fitness(x0);
    
    // std::cout << "ret: " << '\n';
    // for (int i = 0 ; i < ret.size(); i++){
    //     std::cout << ret[i] << '\n';
    // }
    return 0;
}