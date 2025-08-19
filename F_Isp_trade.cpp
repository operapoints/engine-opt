#include <iostream>
#include <chrono>

#include <string>
#include <fstream>

#include <limits>

#include <pagmo/algorithm.hpp>

#include <pagmo/algorithms/sade.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/population.hpp>

#include "jet_calc.h"


void save_to_csv(const std::vector<std::vector<double>> &data, const std::string &filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    for (const auto &row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
}

problem_jet_calc make_pjc_UDP(double k_F, double k_Isp, double k_mass){
    problem_jet_calc pjc_UDP;
    pjc_UDP.k_F = k_F;
    pjc_UDP.k_Isp = k_Isp;
    pjc_UDP.k_mass = k_mass;
    return pjc_UDP;
}

std::pair<double, pagmo::vector_double> get_best_objective(problem_jet_calc& pjc_UDP){
    pagmo::problem pjc{pjc_UDP};
    algorithm algo{gaco(5000,63,1,0,0.01,100,7,5000)};    
    archipelago archi(32u, algo, pjc, 1000u);
    archi.evolve(3);
    archi.wait_check();
    double best_champion_score = std::numeric_limits<double>::infinity();
    int best_isl_idx = 0;
    int isl_idx = 0;
    for (const auto &isl : archi) {
        double champion_score = isl.get_population().champion_f()[0];
        //std::cout << isl.get_population().champion_f()[0] << '\n';
        pagmo::vector_double ch_x = isl.get_population().champion_x();
        //std::cout << "{";
        //for(int i =0; i<ch_x.size(); i++){
        //    std::cout << ch_x[i] << ", ";
        //}
        //std::cout << "}\n";
        if (champion_score<best_champion_score){
            best_champion_score = champion_score;
            best_isl_idx = isl_idx;
        }
        isl_idx++;
    }
    //std::cout << '\n';
    pagmo::vector_double champion_x = archi[best_isl_idx].get_population().champion_x();
    return std::pair<double, pagmo::vector_double>{best_champion_score,champion_x};
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
    problem_jet_calc pjc_obj;
    #ifdef EVAL_JET_CALC
    // Balance
    pagmo::vector_double x1 = {11475.8, 135.212, 1097.15, 0.00501128, 0.0206151, 0.00124931, 0.0299996, 62.201, 0.0154871, 0.0315287, 0.00210709, 0.0188079, };

    pagmo::vector_double x2 = {11095.5, 140.063, 1099.58, 0.00500161, 0.0210375, 0.00149369, 0.0312871, 63.1404, 0.0118135, 0.03243, 0.0014575, 0.0180168, };
    // Max thrust
    pagmo::vector_double x3 = {14684.1, 125.722, 991.887, 0.00504372, 0.0164868, 0.000539293, 0.0299262, 104.827, 0.0155055, 0.0245546, 0.00169747, 0.0209158, };
    pjc_obj.fitness(x1);
    pjc_obj.fitness(x2);
    pjc_obj.fitness(x3);
    #endif

    // This code calculates n_steps separate objective functions, with the relative weightings
    // of mass and Isp ranging linearly. The first steps prioritize mass, the last steps prioritize Isp

    int n_steps = 20;
    double weighting_step = 1/(n_steps-1);
    std::vector<std::vector<double>> data = {};
    std::vector<std::vector<double>> ch_xs = {};

    auto pjc_UDP_start = make_pjc_UDP(0,0,1);
    auto start_res = get_best_objective(pjc_UDP_start);
    double best_mass = -1/start_res.first;
    data.push_back(pjc_obj.evaluate(start_res.second));
    ch_xs.push_back(start_res.second);
    printf("Wrote %u step out of %u\n", 1, n_steps);

    auto pjc_UDP_end = make_pjc_UDP(0,1,0);
    auto end_res = get_best_objective(pjc_UDP_end);
    double best_Isp = -1.*end_res.first;
    data.push_back(pjc_obj.evaluate(end_res.second));
    ch_xs.push_back(end_res.second);
    printf("Wrote %u step out of %u\n", 2, n_steps);

    for (int i=1;i<=n_steps-2;i++){
        auto pjc_UDP = make_pjc_UDP(0,(1./best_Isp)*i*weighting_step,(1./best_mass)*(1-(1*weighting_step)));
        auto res = get_best_objective(pjc_UDP);
        data.insert(data.end()-1,pjc_obj.evaluate(res.second));
        ch_xs.insert(ch_xs.end()-1,res.second);
    printf("Wrote %u step out of %u\n", i+2, n_steps);
    }



    save_to_csv(data, "100_lbf.csv");
    save_to_csv(ch_xs, "100_lbf_designs.csv");

    std::cout << "Done";

    // pagmo::vector_double x0 = {10183,127,1100,0.003,0.023,0.001,0.0250,65,0.02,0.03,0.002,0.025};
    // pagmo::vector_double ret = pjc.fitness(x0);
    
    // std::cout << "ret: " << '\n';
    // for (int i = 0 ; i < ret.size(); i++){
    //     std::cout << ret[i] << '\n';
    // }
    return 0;
}