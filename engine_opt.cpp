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
    // Balance
    pagmo::vector_double x1 = {11475.8, 135.212, 1097.15, 0.00501128, 0.0206151, 0.00124931, 0.0299996, 62.201, 0.0154871, 0.0315287, 0.00210709, 0.0188079, };

    pagmo::vector_double x2 = {15044.9, 133.092, 1099.95, 0.0050458, 0.0157997, 0.000431372, 0.0299949, 109.67, 0.0154697, 0.0244883, 0.00123787, 0.0158204, };
    // Max thrust
    pagmo::vector_double x3 = {14684.1, 125.722, 991.887, 0.00504372, 0.0164868, 0.000539293, 0.0299262, 104.827, 0.0155055, 0.0245546, 0.00169747, 0.0209158, };
    pjc_obj.fitness(x1);
    pjc_obj.fitness(x2);
    pjc_obj.fitness(x3);
    #endif
    pagmo::problem pjc{pjc_obj};
    std::cout << pjc;
    algorithm algo{gaco(10000,63,1,0,0.01,1000,7,1000)};    
    archipelago archi(32u, algo, pjc, 10000u);
    archi.evolve(1);
    archi.wait_check();

    // 6 - Print the fitness of the best solution in each island.
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