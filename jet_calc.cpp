#include "jet_calc.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <initializer_list>
#include <utility>
#include <sstream>


using namespace pagmo;

vector_double::size_type problem_jet_calc::get_nec() const{
    return 0;
}
vector_double::size_type problem_jet_calc::get_nic() const{
    return 19;
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
        double P_3 = P_0*std::pow((T_0+eta_C*D_T_C)/T_0,(gam_c/(gam_c-1)));

        // Compressor constraints
        // Velocities are positive in the direction of rotation
        // phi and psi are calculated according to https://manual.cfturbo.com/en/index.html?md_parameters_axvent.html
        double phi_C = (u_i*(A_Ci/(M_PI*R_Com*R_Com))) / (omega*R_Com);
        double psi_C = (2*C_pc*D_T_C) / (omega*omega*R_Com*R_Com);
        double spec_speed_C = std::pow(phi_C,0.5)/std::pow(psi_C,0.75);
        double spec_dia_C = std::pow(psi_C,0.25)/std::pow(phi_C,0.5);
        // The expression for the cordier line is based on the compressor cordier line at https://manual.cfturbo.com/en/index.html?cordier.html
        double con_cordier_compressor = is_cordier(spec_speed_C, spec_dia_C)?-1:1;// The compressor must be in range of the Cordier line
        double M_Ci_tip = std::pow((u_i*u_i + omega*omega*R_Cit*R_Cit)/(gam_c*R*Ts_2),0.5);
        double con_M_Ci_tip = M_Ci_tip - 0.8;// Mach at compressor inlet less than 0.8
        double beta_Ci_tip = (180/M_PI)*std::atan((omega*R_Cit)/u_i);
        double con_beta_Ci_tip = beta_Ci_tip - 70;// Inlet blade angle less than 70 degrees
        double u_Coth_STAT = (C_pc*D_T_C)/(omega*R_Com);
        double u_Coth_ROT = u_Coth_STAT - (omega*R_Com);
        double u_Coa = compute_u_a(m_dot,(P_3/(R*T_3)), T_3, u_Coth_STAT, gam_c, A_Co);
        if(std::isnan(u_Coa)){
            vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
            return ret;
        }
        double u_Co_ROT = std::pow((u_Coth_STAT-(omega*R_Com))*(u_Coth_STAT-(omega*R_Com)) + u_Coa*u_Coa,0.5);
        double u_Co_STAT = std::pow(u_Coth_STAT*u_Coth_STAT + u_Coa*u_Coa,0.5);
        double u_Ci_ROT = std::pow((u_i*u_i + omega*omega*R_Cit*R_Cit),0.5);
        double Diff_C = 1 - (u_Co_ROT/u_Ci_ROT)+u_Coth_STAT/(2*sigma_C);
        double con_Diff_C = Diff_C - 0.55;// Lieblein diffusion factor less than 0.55
        // u_Coth_ROT is negative because the angle is defined as positive in the direction opposing rotation
        double beta_Co = (180/M_PI)*std::atan(-u_Coth_ROT/u_Coa);
        double con_beta_Co = beta_Co - 60; // Compressor blade exit angle below 40 degrees
        double Ts_3 = T_3 - (u_Co_STAT*u_Co_STAT)/(2*C_pc);
        double a_3 = std::pow(gam_c*R*Ts_3, 0.5);
        double M_Co_ROT = u_Co_ROT/a_3;
        double con_M_Co_ROT = M_Co_ROT - 0.8;// Compressor exit relative Mach less than 0.8
        double M_Co_STAT = std::pow((u_Coa*u_Coa + u_Coth_STAT*u_Coth_STAT)/(gam_c*R*Ts_3),0.5);
        double con_M_Co_STAT = M_Co_STAT - 0.8;// Compressor diffuser inlet Mach less than 0.8
        double sigma_max_Ci = 0.5*omega*omega*rho_C*(R_Cit*R_Cit - R_Cih*R_Cih);
        double com_sigma_max_Ci = FOS_C * sigma_max_Ci - sigma_max_C;// FOS at compressor inlet blade root
        double R_Coh = R_Com - (A_Co/(4*M_PI*R_Com));
        double R_Cot = R_Com + (A_Co/(4*M_PI*R_Com));
        double sigma_max_Co = 0.5*omega*omega*rho_C*(R_Cot*R_Cot - R_Coh*R_Coh);
        double com_sigma_max_Co = FOS_C * sigma_max_Co - sigma_max_C;// FOS at compressor outlet blade root



        // Combustor
        double f = (C_ph*T_4 - C_pc*T_3)/(h_ker - C_ph*T_4);
        // Turbine
        double D_T_T = -(P_spC/(C_ph*(1+f)));
        double T_5 = T_4+D_T_T;
        double P_5 = P_3*std::pow((T_4+(D_T_T/eta_T))/T_4,(gam_h/(gam_h-1)));
        // Turbine constraints
        double A_Ti = M_PI*(R_Tit*R_Tit - R_Tih*R_Tih);
        double R_Tim = 0.5*(R_Tit+R_Tih);
        double u_Tith_STAT = (C_ph)*D_T_T/(omega*R_Tim);
        double u_Tia = compute_u_a(m_dot*(1+f),P_3/(R*T_4), T_4, u_Tith_STAT, gam_h, A_Ti);
        if(std::isnan(u_Tia)){
            vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
            return ret;
        }
        double beta_NGV = (180/M_PI)*std::atan(u_Tith_STAT/u_Tia);
        double con_beta_NGV = beta_NGV - 70;
        // Phi and psi follow the definitions in Dixon and Hall
        // It is calculated with reference to the turbine tip inlet speed. This is apparently the 
        // standard but might apply only to axial turbines where blade speed doesn't vary much axially
        double phi_T = u_Tia / (omega*R_Tit);
        double psi_T = (C_ph*D_T_T) / (omega*omega*R_Tom*R_Tom);
        // Phi and psi constraints are for the smith chart for axial gas turbines in dixon and hall.
        // Smith believed that the losses were proportional to the average kinetic energy in the row,
        // and correlated this with empirically measured losses.
        // If this is true, these phi and psi constraints should hold for any topology.
        double con_phi_T = std::abs(phi_T - (0.5*(min_phi_T+max_phi_T))) - (0.5*(max_phi_T - min_phi_T));// Flow coefficient in appropriate range
        double con_psi_T = std::abs(psi_T - (0.5*(min_psi_T+max_psi_T))) - (0.5*(max_psi_T - min_psi_T));// Loading coefficient in appropriate range
        double u_Ti_STAT = std::pow(u_Tith_STAT*u_Tith_STAT+u_Tia*u_Tia,0.5);
        double D_Ts_NGV = (u_Ti_STAT*u_Ti_STAT) / (2*C_ph);
        double Ts_Ti = T_4 - D_Ts_NGV;
        double DoR = 1 - (D_Ts_NGV/D_T_T);
        double con_DoR = std::abs(0.4 - DoR) - 0.1; // Degree of reaction of turbine between 0.3 and 0.5
        double M_NGVo = std::pow((u_Tith_STAT*u_Tith_STAT+u_Tia*u_Tia)/(gam_h*R*Ts_Ti),0.5);
        double con_M_NGVo = M_NGVo - 0.8; // Mach in NGV exit less than 0.8
        double u_Tith_ROT = u_Tith_STAT - omega*R_Tim;
        double M_Ti = std::pow((u_Tith_ROT*u_Tith_ROT + u_Tia*u_Tia)/(gam_h*R*Ts_Ti),0.5);
        double con_M_Ti = M_Ti - 0.8;// Turbine inlet mach less than 0.8
        double beta_Tim = (180/M_PI)*std::atan(u_Tith_ROT/u_Tia);
        double con_beta_Tim = beta_Tim - 65;// Turbine inlet angle less than 65 degrees
        double u_Toa = compute_u_a(m_dot, P_5/(R*T_5), T_5, omega*R_Tom, gam_h, A_To);
        double beta_Tom = (180/M_PI)*std::atan((omega*R_Tom)/u_Toa);
        double con_beta_Tom = beta_Tom - 65; // Turbine outlet meridional blade angle less than 65 degrees - this shouldn't be active
        double con_T_width = 0.009 - (R_Tit - R_Tih);// Turbine inlet annulus width greater than 9mm
        double con_turbine_diffusion = u_Tia - u_Toa; // Flow must accelerate through turbine to avoid separation
        double Ts_To = T_5 - (u_Toa*u_Toa)/(2*C_ph);
        double con_turbine_outlet_width = 0.09 - A_To/(2*M_PI*R_Tom); // Turbine outlet width more than 9mm
        double R_Tot = R_Tom + A_To/(4*M_PI*R_Tom);
        double M_Tom = std::pow((u_Toa*u_Toa + omega*omega*R_Tot*R_Tot)/(gam_h*R*Ts_To),0.5);
        double con_M_Tom = M_Tom - 0.8; // Relative Mach at turbine exit less than 0.8
        double sigma_max_Ti = 0.5*omega*omega*rho_C*(R_Tit*R_Tit - R_Tih*R_Tih);
        double com_sigma_max_Ti = FOS_T * sigma_max_Ti - sigma_max_T;// FOS at compressor inlet blade root
        double R_Toh = R_Tom - A_To/(4*M_PI*R_Tom);
        double sigma_max_To = 0.5*omega*omega*rho_C*(R_Tot*R_Tot - R_Toh*R_Toh);
        double com_sigma_max_To = FOS_T * sigma_max_To - sigma_max_T;// FOS at compressor inlet blade root



        // Nozzle
        double Ts_6 = T_5*std::pow(Ps_6/P_5,(gam_h-1)/gam_h);
        double a_6 = std::pow(gam_h*R*Ts_6,0.5);
        double u_6 = a_6*std::pow((2/(gam_h-1))*(std::pow((P_5/Ps_6),(gam_h-1)/gam_h)-1),0.5);
        double F = m_dot*((1+f)*u_6 - u_0);
        //Calculate objective
        double Isp = (F/(m_dot*f*9.8066));
        vector_double ret = {Isp, 
            con_beta_Ci_tip, 
            con_beta_Co, 
            con_beta_NGV, 
            con_beta_Tim, 
            con_beta_Tom, 
            con_cordier_compressor, 
            con_Diff_C, 
            con_DoR, 
            con_M_Ci_tip, 
            con_M_Co_ROT, 
            con_M_Co_STAT, 
            con_M_NGVo, 
            con_M_Ti, 
            con_M_Tom, 
            con_phi_T, 
            con_psi_T, 
            con_T_width, 
            con_turbine_diffusion, 
            con_turbine_outlet_width};
        // A physically impossible engine will usually result in a bunch of NaNs,
        // and bad inputs might give Inf due to division by zero.
        // If this happens, return a big penalty
        if (invalid_ret(ret)){
            ret = vector_double(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
        }
        return ret;
    }catch(boost::math::evaluation_error){
        // If the compute_u_a throws an error, return high penalty immediately
        vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
        return ret;
    }
}

inline bool problem_jet_calc::is_cordier(double sigma, double delta) const{
    double ideal_delta{};
    double min_delta{};
    double max_delta{};

    if(sigma<0.05 || sigma>2.5){
        return false;
    }
    if(sigma>=0.05 && sigma < 0.4){
        ideal_delta = 1.1173*std::pow(sigma, -0.963);
    }
    if(sigma>=0.4 && sigma < 0.8){
        ideal_delta = 1.46477*std::pow(sigma, -0.66742);
    }
    if(sigma>=0.8 && sigma < 2.5){
        ideal_delta = 1.61299*std::pow(sigma, -0.23543);
    }
    min_delta = 0.9 * ideal_delta;
    max_delta = 1.1 * ideal_delta;
    return (delta <= max_delta && delta >= min_delta);
}

// Checks if anything in a vector_double is inf or nan
inline bool problem_jet_calc::invalid_ret(vector_double &x) const{
    for (double xi : x){
        if(std::isinf(xi) || std::isnan(xi)){
            return true;
        }
    }
    return false;
}

std::pair<vector_double, vector_double> problem_jet_calc::get_bounds() const{}

// Utility
// Given a 6D double array describing adiabatic flow in a duct,
// computes the axial velocity
// Parameters: mass flow, stag density, stag temp,
// tangential velocity, gamma, area
// Cost is around 200ns
double problem_jet_calc::compute_u_a(double m_dot,
                                    double rho_t,
                                    double T_t,
                                    double u_th,
                                    double gam,
                                    double A) const{
    // For now, if an evaluation_error ever comes up in optimizer runs, try to find a way to handle this without relying on catching
    // try{
        double C_p = 287.0 * (gam/(gam-1));
        // Positive point at which the derivative is zero, the subsonic solution
        // will always be less than this so it is used as the upper bound
        double u_a_max = std::pow(((gam-1)*rho_t*(2*C_p*T_t-(u_th*u_th)))/((1+gam)*rho_t),0.5);
        auto compute_u_a_residual = [=](double u_a)->std::tuple<double,double>{
            double T_s = T_t - ((std::pow(u_a,2)+std::pow(u_th,2))/(2*C_p));
            double tau = T_t/T_s;
            double u_a_residual = m_dot - u_a*A*rho_t*std::pow(tau,(1/(1-gam)));
            double D_u_a_residual = (A*rho_t*(2*C_p*(-1+gam)*T_t-(1+gam)*std::pow(u_a,2)-(-1+gam)*std::pow(u_th,2))*
                                    std::pow(T_t/(T_t-(std::pow(u_a,2)+std::pow(u_th,2))/(2.*C_p)),1/(1-gam)))/
                                    ((-1+gam)*(-2*C_p*T_t+std::pow(u_a,2)+std::pow(u_th,2)));
            return {u_a_residual,D_u_a_residual};
        };
        int digits = 8;
        std::uintmax_t max_iter = 50;
        double u_a = boost::math::tools::newton_raphson_iterate(compute_u_a_residual, 0., 0., u_a_max, digits, max_iter);
        // Error checking that 1. u_a is in fact a root and 2. u_a is not inf or nan
        if (std::abs(std::get<0>(compute_u_a_residual(u_a))-0)>1e-6 || std::isnan(u_a) || std::isinf(u_a)){
            return static_cast<double>(NAN);
        }
        return u_a;
    // }catch(boost::math::evaluation_error){
        // return static_cast<double>(NAN);
    // }
}