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
    pjc_UDP.fitness({7106.95, 130.78, 1099.4, 0.00507257, 0.0335215, 0.00210258, 0.0411106, 60.5994, 0.022099, 0.0496305, 0.00484863, 0.0335061, -41.6881, 3323.66, 0.0113265, 0.0470107, 219.465, });
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