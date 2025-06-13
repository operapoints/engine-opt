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
    // Balance
    pagmo::vector_double x0 = {11666.3, 138.623, 1099.9, 0.00300108, 0.0200928, 0.000985933, 0.029994, 64.4693, 0.0146342, 0.0307445, 0.00159195, 0.0173278, };
    // Max Isp
    pagmo::vector_double x1 = {14906.3, 144.981, 981.448, 0.00301898, 0.0154312, 0.000510475, 0.029995, 108.346, 0.0155338, 0.0246004, 0.00111854, 0.0163724, };
    pjc_obj.fitness(x1);
    pagmo::problem pjc{pjc_obj};
    std::cout << pjc;
    algorithm algo{gaco(50000,63,1,0,0.01,10000,7,10000)};    
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