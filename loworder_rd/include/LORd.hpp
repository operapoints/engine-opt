#pragma once
#include <vector>
#include <Eigen/Dense>

//Defines everything the lord model needs to calculate rotordynamic modes

struct LORd_model_params_t{
    // ---------------- Young's moduli ----------------
    double EC;   // Compressor shaft Young's modulus
    double ES;   // Center shaft Young's modulus (unused)
    double ET;   // Turbine shaft Young's modulus

    // ---------------- Geometry ----------------
    double rC;   // Compressor shaft radius (used for bending MoI: IC = pi*rC^4/4)
    double rS;   // Center shaft radius (used for mS, JpS, JdS)
    double rT;   // Turbine shaft radius (used for bending MoI: IT = pi*rT^4/4)
    double LC;   // Compressor overhang length
    double LS;   // Center shaft length
    double LT;   // Turbine overhang length

    // ---------------- Mass properties ----------------
    double rhoS; // Center shaft density (used for mS = rhoS*pi*rS^2*LS)
    double mC;   // Compressor mass
    double mT;   // Turbine mass
    double JdC;  // Compressor diametral moment of inertia
    double JdT;  // Turbine diametral moment of inertia
    double JpC;  // Compressor polar moment of inertia
    double JpT;  // Turbine polar moment of inertia

    // ---------------- Bearings ----------------
    double KbyC; // Compressor end bearing y axis stiffness
    double KbyT; // Turbine end bearing y axis stiffness
    double KbzC; // Compressor end bearing z axis stiffness
    double KbzT; // Turbine end bearing z axis stiffness

    // ---------------- Operating condition ----------------
    double Omega; // Shaft speed, positive means CW looking at the compressor
};

struct ASTRA_LORd_model_params_t{
    double EC; 
    double ES; 
    double ET; 
    double rC; 
    double rS; 
    double rT; 
    double LC; 
    double LS; 
    double LT; 
    double rhoS;
    double mC; 
    double mT; 
    double JdC;
    double JdT;
    double JpC;
    double JpT;
    double KbyC;
    double KbyT;
    double KbzC;
    double KbzT;
    double Omega;
};

class LORd_model{
public:
    // Returns mode frequencies
    LORd_model();
    // Solves mode eigenvalues
    // Takes 16250ns
    Eigen::Matrix<std::complex<double>,16,1>  mode_eigvals(LORd_model_params_t params);
//private:
    Eigen::Matrix<double, 16, 16> H0;
    // Assembles the Hessian of the Lagrangian about equilibrium
    // Takes about 20ns
    void make_ASTRA_H0(const LORd_model_params_t params);
    void make_H0(const LORd_model_params_t params);
};