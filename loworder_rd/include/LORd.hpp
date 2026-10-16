#pragma once
#include <vector>
#include <Eigen/Dense>

//Defines everything the lord model needs to calculate rotordynamic modes
//EC, ET, IC, IT, JdC, JdS, JdT, JpC, JpS, JpT, KbyC, KbyT, KbzC, KbzT, \
LC, LS, LT, mC, mS, mT, Omega
struct LORd_model_params_t{
    double EC;   // Compressor shaft Young's modulus
    double ET;   // Turbine shaft Young's modulus
    double IC;   // Compressor shaft MoI
    double IT;   // Turbine shaft MoI
    double JdC;  // Compressor diametral moment
    double JdS;  // Center shaft diametral moment
    double JdT;  // Turbine diametral moment
    double JpC;  // Compressor polar moment
    double JpS;  // Center shaft polar moment
    double JpT;  // Turbine polar moment
    double KbyC; // Compressor end bearing y axis stiffness
    double KbyT; // Turbine end bearing y axis stiffness
    double KbzC; // Compressor end bearing z axis stiffness
    double KbzT; // Turbine end bearing z axis stiffness
    double LC;   // Compressor overhang length
    double LS;   // Center shaft length
    double LT;   // Turbine overhang length
    double mC;   // Compressor mass
    double mS;   // Center shaft mass
    double mT;   // Turbine mass
    double Omega;// Shaft speed, positive means CW looking at the compressor
};

class LORd_model{
public:
    // Returns mode frequencies
    std::vector<double> mode_frequencies(LORd_model_params_t params);
//private:
    // Assembles the Hessian of the Lagrangian about equilibrium
    // Takes <40 ns
    Eigen::Matrix<double,16,16> make_H0(LORd_model_params_t params);
};