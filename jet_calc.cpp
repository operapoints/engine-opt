#include "jet_calc.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <initializer_list>
#include <utility>


using namespace pagmo;

vector_double::size_type problem_jet_calc::get_nec() const{
    return auto(0);
}
vector_double::size_type problem_jet_calc::get_nic() const{
    return auto(0);
}

// Calculates the fitness function and associated constraints,
// in the form {fitness, eq..., ineq...}
vector_double problem_jet_calc::fitness(vector_double &x) const{
    auto omega = x[0]; // rad/s - shaft speed
    auto u_i = x[1];   // m/s   - compressor inlet velocity
    auto T_4 = x[2];   // K     - combustor exit total temp
    auto R_Cih = x[3]; // m     - compressor inlet hub radius
    auto R_Cit = x[4]; // m     - compressor inlet tip radius
    auto A_Co = x[5];  // m^2   - compressor outlet area
    auto R_Com = x[6]; // m     - compressor outlet meanline radius
    auto D_T_C = x[7]; // K     - compressor total temperature change
    auto R_Tih = x[8]; // m     - turbine inlet hub radius
    auto R_Tit = x[9]; // m     - Turbine inlet tip radius
    auto A_To = x[10]; // m^2   - Turbine outlet area
    auto R_Tom = x[11];// m     - Turbine exit meanline velocity
    try{
        // Stagnation quantities
        double T_0 = Ts_0 + ((u_0*u_0)/(2*C_pc));
        double P_0 = Ps_0*std::pow((T_0/Ts_0),(gam_c/(gam_c-1)));
        // Duct
        double Ts_2 = T_0 - ((u_i*u_i)/(2*C_pc));
        double Ps_2 = std::pow((Ts_2/T_0),(gam_c/(gam_c-1)));
        double rhos_2 = Ps_2/(R*Ts_2);
        double A_Ci = M_PI*(R_Cit*R_Cit - R_Cih*R_Cih);
        double m_dot = rhos_2*u_i*A_Ci;
        // Compressor
        double P_spC = C_pc*D_T_C;
        double T_3 = T_0 + D_T_C;
        double P_3 = std::pow((T_0+eta_C*D_T_C)/T_0,(gam_c/(gam_c-1)));

        // Compressor constraints
        // Velocities are positive in the direction of rotation
        double phi = u_i / (omega*R_Com);
        double psi = (C_pc*D_T_C) / (omega*omega*R_Com*R_Com);

        // TODO: Set phi and psi constraints

        double M_Ci_tip = std::pow((u_i*u_i + omega*omega*R_Cit*R_Cit)/(gam_c*R*Ts_2),0.5);
        double con_M_Ci_tip = M_Ci_tip - 0.8;
        double beta_Ci_tip = (180/M_PI)*std::atan((omega*R_Cit)/u_i);
        double con_beta_Ci_tip = beta_Ci_tip - 70;
        double u_Coth_STAT = (C_pc*D_T_C)/(omega+R_Com);// Calculated from Euler work eqn
        double u_Coth_ROT = u_Coth_STAT - (omega*R_Com);
        double u_Coa = compute_u_a(m_dot,(P_3/(R*T_3)), T_3, u_Coth_STAT, gam_c, A_Co);
        double u_Co_ROT = std::pow((u_Coth_STAT-(omega*R_Com))*(u_Coth_STAT-(omega*R_Com)) + u_Coa*u_Coa,0.5);
        double u_Co_STAT = std::pow(u_Coth_STAT*u_Coth_STAT + u_Coa*u_Coa,0.5);
        double u_Ci_ROT = std::pow((u_i*u_i + omega*omega*R_Cit*R_Cit),0.5);
        double Diff_C = 1 - (u_Co_ROT/u_Ci_ROT)+u_Coth_STAT/(2*sigma_C);// Liebler diffusion factor
        double con_Diff_C = Diff_C - 0.55;
        double beta_Co = (180/M_PI)*std::atan(-u_Coth_ROT/u_Coa);// u_Coth_ROT is negative because the angle is defined as positive in the direction opposing rotation
        double con_beta_Co = beta_Co - 40; // Compressor blade exit angle below 40 degrees
        double Ts_3 = T_3 - (u_Co_STAT*u_Co_STAT)/(2*C_pc);
        double a_3 = std::pow(gam_c*R*Ts_3, 0.5);
        double M_Co_ROT = u_Co_ROT/a_3;
        double con_M_Co_ROT = M_Co_ROT = 0.8;
        double M_Co_STAT = std::pow((u_Coa*u_Coa + u_Coth_STAT*u_Coth_STAT)/(gam_c*R*Ts_3),0.5);
        double con_M_Co_STAT = M_Co_STAT = 0.8;


        // Combustor
        double f = (C_ph*T_4 - C_pc*T_3)/(h_ker - C_ph*T_4);
        // Turbine
        double D_T_T = -(P_spC/(C_ph*(1+f)));
        double T_5 = T_4+D_T_T;
        double P_5 = std::pow((T_4+(D_T_T/eta_T))/T_4,(gam_h/(gam_h-1)));
        // Turbine constraints

        // Nozzle
        constexpr const double& Ps_6 = Ps_0;
        double Ts_6 = T_5*std::pow(Ps_6/P_5,(gam_h-1)/gam_h);
        double a_6 = std::pow(gam_h*R*Ts_6,0.5);
        double u_6 = a_6*std::pow((2/(gam_h-1))*(std::pow((P_5/Ps_6),(gam_h-1)/gam_h)-1),0.5);
        double F = m_dot*((1+f)*u_6 - u_0);
        //Calculate objective
        double Isp = (F/(m_dot*f*9.8066));

        vector_double ret = {Isp};
        // A physically impossible engine will usually result in a bunch of NaNs,
        // and bad inputs might give Inf due to division by zero.
        // If this happens, return a big penalty
        if (invalid_ret(ret)){
            vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
        }
        return ret;
    }catch(boost::math::evaluation_error){
        // If the compute_u_a throws an error, return high penalty immediately
        vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
        return ret;
    }
}

inline bool problem_jet_calc::invalid_ret(vector_double &x) const{
    for (double xi : x){
        if(std::isinf(xi) || std::isnan(xi)){
            return true;
        }
    }
    return false;
}

std::pair<vector_double, vector_double> problem_jet_calc::get_bounds() const{}


// TODO: make this use newton method instead
// Utility
// Given a 6D double array describing adiabatic flow in a duct,
// computes the axial velocity
// Parameters: mass flow, stag density, stag temp,
// tangential velocity, gamma, area
// Cost is around 200ns a pop
double problem_jet_calc::compute_u_a(double m_dot,
                                    double rho_t,
                                    double T_t,
                                    double u_th,
                                    double gam,
                                    double A) const{
    // The bisect method might raise the boost::math::evaluation_error if the input is non physical.
    // If this happens, catch it in the optimization loop by returning some large penalty.
    double C_p = 287.0 * (gam/(gam-1));
    // Accuracy tolerance
    double eps = 1e-6;
    // a_max is the maximum possible speed of sound, corresponding to no axial velocity.
    // Since only supersonic solutions are necessary, this is used as the upper bound for numerical solution.
    double a_max = std::pow(gam*287.0*(T_t - (u_th*u_th)/(2*C_p)),0.5);
    auto compute_u_a_residual = [=](double u_a){
        double T_s = T_t - ((std::pow(u_a,2)+std::pow(u_th,2))/(2*C_p));
        double tau = T_t/T_s;
        return m_dot - u_a*A*rho_t*std::pow(tau,(1/(1-gam)));
    };
    // For some reason if a NaN comes up at any point bisect will hang indefinitely
    // Max iteration limit is set to 50 because of this
    uintmax_t max_iter = 50;
    auto results = boost::math::tools::bisect(compute_u_a_residual, 
        0.0, // Lower bound
        a_max, // Upper bound
        [=](auto lambda_arg_1, auto lambda_arg_2){return std::abs(lambda_arg_1-lambda_arg_2)<eps;},
        max_iter); // Termination lambda
    double u_a_guess = (results.first + results.second)/2;
    double T_s_guess = T_t - (u_th*u_th + u_a_guess*u_a_guess)/(2*C_p);
    double M_guess = std::pow((u_th*u_th + u_a_guess*u_a_guess)/(gam*287.0*T_s_guess),0.5);
    // If the found solution is supersonic, try again with the upper bound restricted
    // To the lower solution
    // The idea here is that if the flow has solutions, there is no reason to design for the
    // supersonic case due to resulting losses
    if(M_guess>1){
        auto results = boost::math::tools::bisect(compute_u_a_residual, 
            0.0, // Lower bound
            u_a_guess - eps, // Upper bound; subtract eps so that a sign change is guaranteed
            [=](auto lambda_arg_1, auto lambda_arg_2){return std::abs(lambda_arg_1-lambda_arg_2)<eps;},
            max_iter); // Termination lambda
        double u_a = (results.first+results.second) / 2;
        if (std::isnan(u_a) || std::isinf(u_a)){
            throw boost::math::evaluation_error("Computed u_a value was inf or nan");
        }
        return u_a;
    // If the original solution is good, return it
    }else{
        if (std::isnan(u_a_guess) || std::isinf(u_a_guess)){
            throw boost::math::evaluation_error("Computed u_a value was inf or nan");
        }
        return u_a_guess;
    }
}