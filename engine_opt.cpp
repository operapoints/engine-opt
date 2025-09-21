#include <iostream>
#include <chrono>

#include <limits>

#include <pagmo/algorithm.hpp>

#include <pagmo/algorithms/sade.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/population.hpp>

#include "jet_calc.h"


problem_jet_calc make_pjc_UDP(double k_F, double k_Isp, double k_mass){
    problem_jet_calc pjc_UDP;
    pjc_UDP.k_F = k_F;
    pjc_UDP.k_Isp = k_Isp;
    pjc_UDP.k_mass = k_mass;
    return pjc_UDP;
}

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
    problem_jet_calc pjc_UDP = make_pjc_UDP(0.,0.,1.);
    #ifdef EVAL_JET_CALC
    //centrif constraint
    pjc_UDP.fitness({5396.41, 92.9053, 1099.19, 0.0198369, 0.0469823, 0.00291403, 0.0394745, 42.5561, 0.0202302, 0.0603966, 0.00904774, 0.0367979, -31.3476, 2441.96, 0.0177454, 0.0556191, 81.6692, 0.0704051, });
    //centrif hub constraint
    pjc_UDP.fitness({4612.4, 94.2008, 1099.69, 0.0283437, 0.0554228, 0.00384227, 0.0418329, 36.2184, 0.0172472, 0.0698659, 0.00622032, 0.0343891, -25.2571, 2047.3, 0.0206522, 0.0710636, 76.4735, 0.0631419, });
    //Pure Isp
    pjc_UDP.fitness({3478.44, 55.5469, 988.474, 0.00683918, 0.0432945, 0.00213543, 0.105641, 94.2251, 0.0767441, 0.08577, 0.00581199, 0.070395, -58.198, 889.036, 0.0081343, 0.199999, 59.4054, 0.199999, });
    #endif
    pagmo::problem pjc{pjc_UDP};
    std::cout << pjc;
    algorithm algo{gaco(5000,63,1,0,0.01,100,7,5000)};    
    archipelago archi(32u, algo, pjc, 3000u);
    archi.evolve(3);
    archi.wait_check();

    // 6 - Print the fitness of the best solution in each island.
    double best_champion_score = std::numeric_limits<double>::infinity();
    int best_isl_idx = 0;
    int isl_idx = 0;
    for (const auto &isl : archi) {
        auto champion_score = isl.get_population().champion_f()[0];
        std::cout << isl.get_population().champion_f()[0] << '\n';
        pagmo::vector_double ch_x = isl.get_population().champion_x();
        std::cout << "{";
        for(int i =0; i<ch_x.size(); i++){
            std::cout << ch_x[i] << ", ";
        }
        std::cout << "}\n";
        if (champion_score<best_champion_score){
            best_champion_score = champion_score;
            best_isl_idx = isl_idx;
        }
        isl_idx++;
    }
    std::cout << '\n';
    auto champion_x = archi[best_isl_idx].get_population().champion_x();
    std::cout << "Best score: " << best_champion_score << '\n';
    std::cout << "{";
    for(int i =0; i<champion_x.size(); i++){
        std::cout << champion_x[i] << ", ";
    }
    std::cout << "}\n";



    // pagmo::vector_double x0 = {10183,127,1100,0.003,0.023,0.001,0.0250,65,0.02,0.03,0.002,0.025};
    // pagmo::vector_double ret = pjc.fitness(x0);
    
    // std::cout << "ret: " << '\n';
    // for (int i = 0 ; i < ret.size(); i++){
    //     std::cout << ret[i] << '\n';
    // }
    return 0;
}