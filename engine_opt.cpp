#include <iostream>
#include <chrono>

#include <pagmo/algorithm.hpp>

#include <pagmo/algorithms/sade.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/problem.hpp>

#include "jet_calc.h"

int main(){
    problem_jet_calc pjc;
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono:: high_resolution_clock::now();
    double computed_u_a;
    try{
        auto start = std::chrono::high_resolution_clock::now();
        computed_u_a = pjc.compute_u_a(0.25,0.4,1100,200,1.36,0.01);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "computed u_a: " << std::setprecision(10) << computed_u_a << '\n';
    }
    catch(boost::math::evaluation_error& e){
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Invalid input to u_a:\n";
        std::cout << e.what() << '\n';
    }
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken: " << elapsed.count() * 1e+6 << " microseconds\n";
    return 0;
}