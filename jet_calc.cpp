#include "jet_calc.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <initializer_list>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <limits>


using namespace pagmo;

namespace {
struct TurbineSpanMach {
    double inlet_hub = 0., inlet_shroud = 0.;
    double ngv_hub = 0., ngv_shroud = 0.;
    double outlet_hub = 0., outlet_shroud = 0.;
};

// Same free-vortex/uniform-axial convention as the blading evaluation.
// Local static temperature is a physical scalar, independent of reference frame.
bool turbine_span_mach(TurbineSpanMach& mach, double omega,
                       double r_ih, double r_it, double r_oh, double r_ot,
                       double c_theta_mean, double c_axial_in, double c_axial_out,
                       double T_total_in, double T_total_out, double cp,
                       double gamma, double gas_R) {
    if (!(r_ih > 0. && r_it > r_ih && r_oh > 0. && r_ot > r_oh)) return false;
    const double vortex = 0.5*(r_ih+r_it)*c_theta_mean;
    const auto inlet_mach = [&](double r, double& absolute_mach) {
        const double swirl = vortex/r;
        const double ts = T_total_in-(c_axial_in*c_axial_in+swirl*swirl)/(2*cp);
        if (!(ts>0.)) return std::numeric_limits<double>::quiet_NaN();
        const double sound_speed = std::sqrt(gamma*gas_R*ts);
        absolute_mach = std::hypot(c_axial_in,swirl)/sound_speed;
        return std::hypot(c_axial_in,swirl-omega*r)/sound_speed;
    };
    // Both rotors have zero absolute exit swirl and uniform exit axial velocity.
    const double ts_out = T_total_out-c_axial_out*c_axial_out/(2*cp);
    if (!(ts_out>0.)) return false;
    const double a_out = std::sqrt(gamma*gas_R*ts_out);
    mach.inlet_hub = inlet_mach(r_ih,mach.ngv_hub);
    mach.inlet_shroud = inlet_mach(r_it,mach.ngv_shroud);
    mach.outlet_hub = std::hypot(c_axial_out,omega*r_oh)/a_out;
    mach.outlet_shroud = std::hypot(c_axial_out,omega*r_ot)/a_out;
    return std::isfinite(mach.inlet_hub) && std::isfinite(mach.inlet_shroud) &&
           std::isfinite(mach.ngv_hub) && std::isfinite(mach.ngv_shroud) &&
           std::isfinite(mach.outlet_hub) && std::isfinite(mach.outlet_shroud);
}

void record_turbine_mach(std::unordered_map<std::string,double>& out,
                         const std::string& stage, const TurbineSpanMach& mach,
                         double mean_in, double mean_out, double mean_ngv) {
    out["[] "+stage+" NGV exit Mach (hub)"] = mach.ngv_hub;
    out["[] "+stage+" NGV exit Mach (shroud)"] = mach.ngv_shroud;
    out["[] "+stage+" NGV exit Mach (mean)"] = mean_ngv;
    out["[] "+stage+" inlet relative Mach (hub)"] = mach.inlet_hub;
    out["[] "+stage+" inlet relative Mach (shroud)"] = mach.inlet_shroud;
    out["[] "+stage+" outlet relative Mach (hub)"] = mach.outlet_hub;
    out["[] "+stage+" outlet relative Mach (shroud)"] = mach.outlet_shroud;
    out["[] "+stage+" inlet relative Mach (mean)"] = mean_in;
    out["[] "+stage+" outlet relative Mach (mean)"] = mean_out;
}

// Ideal free vortex: r*C_theta is constant and axial velocity is constant across
// the span. Meanline continuity remains the existing quasi-1D approximation.
// Rotor-inlet beta follows rotation; outlet beta is the positive opposing angle,
// matching the original HPT evaluation's sign convention. All angles in degrees.
void evaluate_turbine_blading(std::unordered_map<std::string,double>& out,
                             const std::string& stage, double omega,
                             double r_ih, double r_it, double r_oh, double r_ot,
                             double c_theta_mean, double c_axial_in, double c_axial_out,
                             double T_total_in, double cp) {
    const double r_mean = 0.5*(r_ih+r_it);
    const double vortex = r_mean*c_theta_mean;
    bool valid = true;
    for (int i=0;i<3;++i) {
        const std::string span = i==0 ? "hub" : (i==1 ? "mean" : "tip");
        const double r_in = i==0 ? r_ih : (i==1 ? r_mean : r_it);
        const double r_out = i==0 ? r_oh : (i==1 ? 0.5*(r_oh+r_ot) : r_ot);
        if (!(r_in>0. && r_out>0.)) { valid=false; continue; }
        const double swirl = vortex/r_in;
        const double ts = T_total_in-(swirl*swirl+c_axial_in*c_axial_in)/(2*cp);
        if (!(ts>0.)) valid=false;
        const auto key = [&](const std::string& units, const std::string& name) {
            return units+" "+stage+" "+name+" ("+span+")";
        };
        out[key("[m/s]","NGV exit swirl")] = swirl;
        out[key("[m/s]","NGV exit axial velocity")] = c_axial_in;
        out[key("[K]","NGV exit static temperature")] = ts;
        out[key("[deg]","NGV exit angle")] = 180./M_PI*std::atan2(swirl,c_axial_in);
        out[key("[deg]","rotor inlet angle")] = 180./M_PI*std::atan2(swirl-omega*r_in,c_axial_in);
        out[key("[deg]","rotor outlet angle")] = 180./M_PI*std::atan2(omega*r_out,c_axial_out);
    }
    // Failed radial states remain visible; inspection must not change fitness.
    out["[0/1] "+stage+" free-vortex span valid"] = valid ? 1. : 0.;
}
}




vector_double::size_type problem_jet_calc::get_nec() const{
    return 0;
}
vector_double::size_type problem_jet_calc::get_nic() const{
    return aft_fan_ ? 42 : 24;
}

std::pair<vector_double, vector_double> problem_jet_calc::get_bounds() const{
    return {vector_double(get_design_size(), 0.0), vector_double(get_design_size(), 1.0)};
}

std::pair<vector_double, vector_double> problem_jet_calc::get_physical_bounds() const {
    vector_double lo(physical_lower, physical_lower + design_size);
    vector_double hi(physical_upper, physical_upper + design_size);
    if (aft_fan_) {
        lo.insert(lo.end(), aft_lower, aft_lower + 8);
        hi.insert(hi.end(), aft_upper, aft_upper + 8);
    }
    return {lo, hi};
}

vector_double problem_jet_calc::normalize(const vector_double &physical) const {
    const auto bounds = get_physical_bounds();
    if (physical.size() != get_design_size())
        throw std::invalid_argument("Physical design dimension does not match engine mode");
    auto x = physical;
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] = (x[i] - bounds.first[i]) / (bounds.second[i] - bounds.first[i]);
        if (!std::isfinite(x[i]) || x[i] < 0. || x[i] > 1.)
            throw std::invalid_argument("Physical design lies outside engine bounds");
    }
    return x;
}

vector_double problem_jet_calc::denormalize(const vector_double &x) const {
    const auto bounds = get_physical_bounds();
    if (x.size() != get_design_size())
        throw std::invalid_argument("Normalized design dimension does not match engine mode");
    auto physical = x;
    for (std::size_t i = 0; i < x.size(); ++i)
        physical[i] = bounds.first[i] + x[i] * (bounds.second[i] - bounds.first[i]);
    return physical;
}

// Calculates the fitness function and associated constraints,
// in the form {fitness, eq..., ineq...}
// x contains normalized design variables; all engine calculations use SI units.
vector_double problem_jet_calc::fitness(const vector_double &x, std::unordered_map<std::string,double> *results_ptr) const{
    if (results_ptr) results_ptr->clear();
    const auto penalty = [this]() {
        return vector_double(1 + get_nec() + get_nic(), 1.e6);
    };
    if (x.size() != get_design_size())
        throw std::invalid_argument("Design dimension does not match engine mode");
    for (double value : x)
        if (!std::isfinite(value) || value < 0. || value > 1.) return penalty();
    if (aft_fan_ && (!(eta_fan > 0. && eta_fan <= 1.) ||
                     !std::isfinite(M_max_fan) || !(M_max_fan > 0.)))
        throw std::invalid_argument("Invalid fan efficiency or Mach limit");
    // Decode directly: avoid allocating bounds/vectors in the fitness hot path.
    const auto physical = [&x](std::size_t i) {
        const double lo = i < design_size ? physical_lower[i] : aft_lower[i-design_size];
        const double hi = i < design_size ? physical_upper[i] : aft_upper[i-design_size];
        return lo + x[i]*(hi-lo);
    };
    const double omega = physical(0); // rad/s - shaft speed
    const double u_i = physical(1);   // m/s   - compressor inlet velocity
    const double T_4 = physical(2);   // K     - combustor exit total temp
    const double R_Cih = physical(3); // m     - compressor inlet hub radius
    const double R_Cit = physical(4); // m     - compressor inlet tip radius
    const double A_Co = physical(5);  // m^2   - compressor outlet area
    const double R_Com = physical(6); // m     - compressor outlet meanline radius
    const double D_T_C = physical(7); // K     - compressor total temperature change
    const double R_Tih = physical(8); // m     - turbine inlet hub radius
    const double R_Tit = physical(9); // m     - turbine inlet tip radius
    const double A_To = physical(10);// m^2   - turbine outlet area
    const double R_Tom = physical(11);// m     - turbine exit meanline radius
    const double D_T_LPT = aft_fan_ ? physical(12) : 0.; // K, negative extraction
    const double omega_LPT = aft_fan_ ? physical(13) : 0.; // rad/s
    const double A_LPTo = aft_fan_ ? physical(14) : 0.; // m^2
    const double R_LPTom = aft_fan_ ? physical(15) : 0.; // m
    const double u_fani = aft_fan_ ? physical(16) : 0.; // m/s
    const double fan_blade_width = aft_fan_ ? physical(17) : 0.; // radial span, m
    const double R_LPTih = aft_fan_ ? physical(18) : 0.; // m, LPT inlet hub
    const double R_LPTit = aft_fan_ ? physical(19) : 0.; // m, LPT inlet tip
    // Reject nonphysical annuli and zero-work/zero-flow singularities.
    if (!(omega > 0. && u_i > 0. && D_T_C > 0. && R_Com > 0. &&
          R_Cit > R_Cih && R_Tit > R_Tih)) return penalty();
    if (aft_fan_ && !(D_T_LPT < 0. && omega_LPT > 0. && A_LPTo > 0. &&
                     R_LPTom > 0. && u_fani > 0. && fan_blade_width > 0. && R_LPTit > R_LPTih))
        return penalty();
    // Variable subscript convention:
    // if u:
    // u_<Component : n><Component axial location : 1><Velocity direction : 1 - 2>_<component radial location : n>
    // otherwise:
    // <var symbol : n>_<Component : n><component axial location : 1><component radial location : n>
    try{
        // Stagnation quantities
        double T_0 = Ts_0 + ((u_0*u_0)/(2*C_pc));
        double P_0 = Ps_0*std::pow((T_0/Ts_0),(gam_c/(gam_c-1)));
        // Duct
        double Ts_2 = T_0 - ((u_i*u_i)/(2*C_pc));
        double Ps_2 = P_0 * std::pow((Ts_2/T_0),(gam_c/(gam_c-1)));
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
        double con_cordier_compressor = is_cordier(spec_speed_C, spec_dia_C);// The compressor must be in range of the Cordier line
        double M_Ci_tip = std::pow((u_i*u_i + omega*omega*R_Cit*R_Cit)/(gam_c*R*Ts_2),0.5);
        double con_M_Ci_tip = M_Ci_tip - 1.2;// Main-branch compressor inlet Mach limit
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
        double Diff_C = 1 - (u_Co_ROT/u_Ci_ROT)+u_Coth_STAT/(2*sigma_C*u_Ci_ROT);
        double con_Diff_C = Diff_C - 0.55;// Lieblein diffusion factor less than 0.55
        // u_Coth_ROT is negative because the angle is defined as positive in the direction opposing rotation
        double beta_Co = (180/M_PI)*std::atan(-u_Coth_ROT/u_Coa);
        double con_beta_Co = beta_Co - 60; // Compressor blade exit angle below 40 degrees
        double Ts_3 = T_3 - (u_Co_STAT*u_Co_STAT)/(2*C_pc);
        double Ps_3 = std::pow((Ts_3/T_3),gam_c/(gam_c-1)) * P_3;
        double a_3 = std::pow(gam_c*R*Ts_3, 0.5);
        double M_Co_ROT = u_Co_ROT/a_3;
        double con_M_Co_ROT = M_Co_ROT - 0.8;// Compressor exit relative Mach less than 0.8
        double M_Co_STAT = std::pow((u_Coa*u_Coa + u_Coth_STAT*u_Coth_STAT)/(gam_c*R*Ts_3),0.5);
        double con_M_Co_STAT = M_Co_STAT - 0.8;// Compressor diffuser inlet Mach less than 0.8
        double sigma_max_Ci = 0.5*omega*omega*rho_C*(R_Cit*R_Cit - R_Cih*R_Cih);
        double con_sigma_max_Ci = FOS_C * sigma_max_Ci - sigma_max_C;// FOS at compressor inlet blade root
        double R_Coh = R_Com - (A_Co/(4*M_PI*R_Com));
        double con_R_Coh = -1*R_Coh;
        double R_Cot = R_Com + (A_Co/(4*M_PI*R_Com));
        double sigma_max_Co = 0.5*omega*omega*rho_C*(R_Cot*R_Cot - R_Coh*R_Coh);
        double con_sigma_max_Co = FOS_C * sigma_max_Co - sigma_max_C;// FOS at compressor outlet blade root

        #ifdef EVAL_JET_CALC
        double beta_Di = (180/M_PI)*std::atan(u_Coth_STAT / u_Coa);

        double R_Co_x = R_Com + 0.003;
        double u_Coth_x_STAT = u_Coth_STAT * (R_Com/R_Co_x);
        double u_Coth_x_ROT = u_Coth_x_STAT - omega*R_Co_x;
        double u_Coa_x = compute_u_a(m_dot,(P_3/(R*T_3)), T_3, u_Coth_x_STAT, gam_c, A_Co);
        double beta_Co_x = (180/M_PI)*std::atan(-u_Coth_x_ROT/u_Coa_x);
        double beta_Di_x = (180/M_PI)*std::atan(u_Coth_x_STAT / u_Coa_x);
        #endif

        // Combustor
        double f = (C_ph*T_4 - C_pc*T_3)/(h_ker - C_ph*T_4);
        if (!(f > 0. && std::isfinite(f))) return penalty();
        double P_4 = 0.94*P_3;
        // Turbine
        double D_T_T = -(P_spC/(C_ph*(1+f)));
        double T_5 = T_4+D_T_T;
        double P_5 = P_4*std::pow((T_4+(D_T_T/eta_T))/T_4,(gam_h/(gam_h-1)));
        // Turbine constraints
        double A_Ti = M_PI*(R_Tit*R_Tit - R_Tih*R_Tih);
        double R_Tim = 0.5*(R_Tit+R_Tih);
        double u_Tith_STAT = (C_ph*(-1*D_T_T))/(omega*R_Tim);// Since D_T_T is negative, it has to be flipped to calculate the absolute NGV exit u_th
        double u_Tia = compute_u_a(m_dot*(1+f),P_4/(R*T_4), T_4, u_Tith_STAT, gam_h, A_Ti);
        if(std::isnan(u_Tia)){
            vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
            return ret;
        }
        double beta_NGV = (180/M_PI)*std::atan(u_Tith_STAT/u_Tia);
        double con_beta_NGV = beta_NGV - 70;
        // Phi and psi follow the definitions in Dixon and Hall
        // It is calculated with reference to the turbine tip inlet speed. This is apparently the 
        // standard but might apply only to axial turbines where blade speed doesn't vary much axially
        double R_T_mean = 0.5*(R_Tim+R_Tom); // NOTE: this is an arbitrarily selected "mean" radius to
        // match the turbine phi and psi estimation to existing empirical efficiency data assuming a constant
        // meanline radius. Variable radius may cause the map from phi and psi to efficiency to change.
        double phi_T = u_Tia / (omega*R_T_mean);
        double psi_T = (C_ph*(-D_T_T)) / (omega*omega*R_T_mean*R_T_mean);//D_T_T is flipped here also to follow work coefficient convention
        // Phi and psi constraints are for the smith chart for axial gas turbines in dixon and hall.
        // Smith believed that the losses were proportional to the average kinetic energy in the row,
        // and correlated this with empirically measured losses.
        // If this is true, these phi and psi constraints should hold for any topology.
        // In reality, higher loading will force longer passage length to prevent separation losses
        // But as long as psi remains on the low side, the topology will be axial and the assumptions of Smith will hold
        double con_phi_T = std::abs(phi_T - (0.5*(min_phi_T+max_phi_T))) - (0.5*(max_phi_T - min_phi_T));// Flow coefficient in appropriate range
        double con_psi_T = std::abs(psi_T - (0.5*(min_psi_T+max_psi_T))) - (0.5*(max_psi_T - min_psi_T));// Loading coefficient in appropriate range
        double u_Ti_STAT = std::pow(u_Tith_STAT*u_Tith_STAT+u_Tia*u_Tia,0.5);
        double D_Ts_NGV = -1*(u_Ti_STAT*u_Ti_STAT) / (2*C_ph);
        double Ts_Ti = T_4 + D_Ts_NGV;
        double Ps_Ti = P_4 * std::pow(Ts_Ti/T_4, (gam_h)/(gam_h-1));
        double M_NGVo = std::pow((u_Tith_STAT*u_Tith_STAT+u_Tia*u_Tia)/(gam_h*R*Ts_Ti),0.5);
        double u_Tith_ROT = u_Tith_STAT - omega*R_Tim;
        double M_Ti = std::pow((u_Tith_ROT*u_Tith_ROT + u_Tia*u_Tia)/(gam_h*R*Ts_Ti),0.5);
        double beta_Tim = (180/M_PI)*std::atan(u_Tith_ROT/u_Tia);
        double con_beta_Tim = std::abs(beta_Tim) - 65;// Turbine inlet angle less than 65 degrees
        double u_Toa = compute_u_a(m_dot * (1+f), P_5/(R*T_5), T_5, 0, gam_h, A_To);
        if(std::isnan(u_Toa)){
            vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
            return ret;
        }
        double beta_Tom = (180/M_PI)*std::atan((omega*R_Tom)/u_Toa);
        double con_beta_Tom = beta_Tom - 65; // Turbine outlet meridional blade angle less than 65 degrees - this shouldn't be active
        double con_T_width = 0.009 - (R_Tit - R_Tih);// Turbine inlet annulus width greater than 9mm
        // double con_turbine_diffusion = u_Tia - u_Toa; // Flow must accelerate through turbine to avoid separation
        double Ts_To = T_5 - (u_Toa*u_Toa)/(2*C_ph);
        double con_turbine_outlet_width = 0.009 - A_To/(2*M_PI*R_Tom); // Turbine outlet width more than 9mm
        double R_Tot = R_Tom + A_To/(4*M_PI*R_Tom);
        double M_Tom = std::pow((u_Toa*u_Toa + omega*omega*R_Tom*R_Tom)/(gam_h*R*Ts_To),0.5);
        double M_Toa = u_Toa / std::pow((gam_h*R*Ts_To),0.5);

        // Degree of reaction calculation
        double D_Ts_T = Ts_To - T_4;
        double DoR = 1 - (D_Ts_NGV/D_Ts_T);
        double con_DoR = std::abs(0.4 - DoR) - 0.1; // Degree of reaction of turbine between 0.3 and 0.5
        // TODO: Change this back when done
        double sigma_max_Ti = 0.5*omega*omega*rho_T*(R_Tit*R_Tit - R_Tih*R_Tih);
        double con_sigma_max_Ti = FOS_T * sigma_max_Ti - sigma_max_T;// FOS at compressor inlet blade root
        double R_Toh = R_Tom - A_To/(4*M_PI*R_Tom);
        double con_R_Toh = 0.004 - R_Toh;// Inner radius of turbine outlet must be at least 4mm
        double sigma_max_To = 0.5*omega*omega*rho_T*(R_Tot*R_Tot - R_Toh*R_Toh);
        double con_sigma_max_To = FOS_T * sigma_max_To - sigma_max_T;// FOS at compressor inlet blade root
        TurbineSpanMach hpt_mach;
        if (!turbine_span_mach(hpt_mach,omega,R_Tih,R_Tit,R_Toh,R_Tot,
                              u_Tith_STAT,u_Tia,u_Toa,T_4,T_5,C_ph,gam_h,R))
            return penalty();
        // Opposing-sign hub/shroud migration, allowing a stationary boundary.
        // Normalize the product (m^2) by a 10 mm reference length squared.
        const double con_HPT_migration = (R_Tih-R_Toh)*(R_Tit-R_Tot)/(0.01*0.01);
        double con_turbine_tip_geom = R_Tit - R_Tot;
        double con_turbine_hub_geom = R_Toh - R_Tih;

        double T_6 = T_5, P_6 = P_5;
        double R_LPTot = 0., R_LPToh = 0., R_fan = 0., m_dot_fan = 0.;
        double D_T_fan = 0., T_fano = 0., P_fano = 0., F_fan = 0., u_LPToa = 0.;
        double phi_LPT = 0., psi_LPT = 0., phi_fan = 0., psi_fan = 0.;
        double spec_speed_fan = 0., spec_dia_fan = 0.;
        double con_LPT_backpres = 0., con_sigma_max_LPTo = 0., con_Mach_fan_tip = 0.;
        double con_sigma_max_fano = 0., con_phi_LPT = 0., con_psi_LPT = 0.;
        double con_cordier_fan = 0., M_LPTo = 0.;
        TurbineSpanMach lpt_mach;
        double u_LPTia = 0., u_LPTith_STAT = 0., Ts_LPTi = 0., Ps_LPTi = 0.;
        double Ts_LPTo = 0., DoR_LPT = 0., R_LPT_mean = 0.;
        double M_LPTNGVo = 0., M_LPTi = 0., M_LPToa = 0.;
        double beta_LPTNGV = 0., beta_LPTim = 0., beta_LPTom = 0.;
        double con_DoR_LPT = 0., con_LPT_migration = 0.;
        double con_LPT_width = 0., con_LPT_outlet_width = 0., con_sigma_max_LPTi = 0.;
        if (aft_fan_) {
            // LPT: same inlet/rotor/outlet convention as HPT; zero exit swirl.

            T_6 = T_5 + D_T_LPT;
            if (!(T_6 > 0. && T_5 + D_T_LPT/eta_T > 0.)) return penalty();
            P_6 = P_5*std::pow((T_5+(D_T_LPT/eta_T))/T_5,(gam_h/(gam_h-1)));
            con_LPT_backpres = Ps_6 - P_6;
            double Ts_fani = T_0 - ((u_fani*u_fani)/(2*C_pc));
            double Ps_fani = P_0 * std::pow((Ts_fani/T_0),(gam_c/(gam_c-1)));
            double rhos_fani = Ps_fani/(R*Ts_fani);

            R_LPTot = R_LPTom + A_LPTo/(4*M_PI*R_LPTom);
            R_LPToh = R_LPTom - A_LPTo/(4*M_PI*R_LPTom);


            double sigma_max_LPTo = 0.5*omega_LPT*omega_LPT*rho_T*(R_LPTot*R_LPTot - R_LPToh*R_LPToh);
            con_sigma_max_LPTo = FOS_T * sigma_max_LPTo - sigma_max_T;// FOS at compressor inlet blade root

            R_fan = fan_blade_width + R_LPTot;
            // Direct relative Mach residual stays finite when axial speed alone exceeds the limit.
            con_Mach_fan_tip = std::sqrt((u_fani*u_fani + omega_LPT*omega_LPT*R_fan*R_fan)
                                                /(gam_c*R*Ts_fani)) - M_max_fan;
            double A_fan = M_PI*(R_fan*R_fan - R_LPTot*R_LPTot);
            m_dot_fan = A_fan * u_fani* rhos_fani;
            D_T_fan = -(C_ph*m_dot*(1+f)*D_T_LPT)/(C_pc*m_dot_fan);
            T_fano = D_T_fan+T_0;
            P_fano = P_0*std::pow((T_0+eta_fan*D_T_fan)/T_0,(gam_c/(gam_c-1)));
            // Architecture assumption: an independently supported outer ring carries
            // fan centrifugal loads; its LPT coupling transmits torque only.
            // Do not add fan radial load to the LPT blade-root stress estimate.
            // Ring/attachment strength and mass are not sized by this surrogate.
            double sigma_max_fano = 0.5*omega_LPT*omega_LPT*rho_T*(R_fan * R_fan - R_LPTot*R_LPTot);
            con_sigma_max_fano = FOS_T * sigma_max_fano - sigma_max_T;// FOS at compressor inlet blade root

            // Fan nozzle
            double Ts_fannoz = T_fano*std::pow(Ps_6/P_fano,(gam_c-1)/gam_c);
            double a_fannoz = std::pow(gam_c*R*Ts_fannoz,0.5);
            double u_fannoz = a_fannoz*std::pow((2/(gam_c-1))*(std::pow((P_fano/Ps_6),(gam_c-1)/gam_c)-1),0.5);
            F_fan = m_dot_fan*(u_fannoz - u_0);




            // LPT inlet (NGV exit), analogous to HPT inlet.
            const double A_LPTi = M_PI*(R_LPTit*R_LPTit-R_LPTih*R_LPTih);
            const double R_LPTim = 0.5*(R_LPTit+R_LPTih);
            u_LPTith_STAT = C_ph*(-D_T_LPT)/(omega_LPT*R_LPTim);
            u_LPTia = compute_u_a(m_dot*(1+f), P_5/(R*T_5), T_5,
                                  u_LPTith_STAT, gam_h, A_LPTi);
            if (!std::isfinite(u_LPTia)) return penalty();
            beta_LPTNGV = 180./M_PI*std::atan2(u_LPTith_STAT,u_LPTia);
            const double u_LPTi_STAT = std::hypot(u_LPTith_STAT,u_LPTia);
            const double D_Ts_LPTNGV = -u_LPTi_STAT*u_LPTi_STAT/(2*C_ph);
            Ts_LPTi = T_5+D_Ts_LPTNGV;
            Ps_LPTi = P_5*std::pow(Ts_LPTi/T_5,gam_h/(gam_h-1));
            M_LPTNGVo = u_LPTi_STAT/std::sqrt(gam_h*R*Ts_LPTi);
            const double u_LPTith_ROT = u_LPTith_STAT-omega_LPT*R_LPTim;
            beta_LPTim = 180./M_PI*std::atan2(u_LPTith_ROT,u_LPTia);
            M_LPTi = std::hypot(u_LPTith_ROT,u_LPTia)/std::sqrt(gam_h*R*Ts_LPTi);
            // Same mean-radius definitions, coefficient ranges and residuals as HPT.
            R_LPT_mean = 0.5*(R_LPTim+R_LPTom);
            phi_LPT = u_LPTia/(omega_LPT*R_LPT_mean);
            psi_LPT = C_ph*(-D_T_LPT)/(omega_LPT*omega_LPT*R_LPT_mean*R_LPT_mean);
            con_phi_LPT = std::abs(phi_LPT-0.5*(min_phi_T+max_phi_T))-0.5*(max_phi_T-min_phi_T);
            con_psi_LPT = std::abs(psi_LPT-0.5*(min_psi_T+max_psi_T))-0.5*(max_psi_T-min_psi_T);

            // LPT outlet, analogous to HPT outlet.
            u_LPToa = compute_u_a(m_dot*(1+f), P_6/(R*T_6), T_6, 0., gam_h, A_LPTo);
            if (!std::isfinite(u_LPToa)) return penalty();
            Ts_LPTo = T_6-u_LPToa*u_LPToa/(2*C_ph);
            M_LPTo = std::hypot(u_LPToa,omega_LPT*R_LPTom)/std::sqrt(gam_h*R*Ts_LPTo);
            if (!turbine_span_mach(lpt_mach,omega_LPT,R_LPTih,R_LPTit,R_LPToh,R_LPTot,
                                  u_LPTith_STAT,u_LPTia,u_LPToa,T_5,T_6,C_ph,gam_h,R))
                return penalty();
            con_LPT_migration = (R_LPTih-R_LPToh)*(R_LPTit-R_LPTot)/(0.01*0.01);
            M_LPToa = u_LPToa/std::sqrt(gam_h*R*Ts_LPTo);
            beta_LPTom = 180./M_PI*std::atan2(omega_LPT*R_LPTom,u_LPToa);
            // Angle limits remain diagnostic-only for both turbines, as in main.
            // Same total-to-static reference as HPT: negligible stator-entry KE.
            const double D_Ts_LPT = Ts_LPTo-T_5;
            DoR_LPT = 1.-D_Ts_LPTNGV/D_Ts_LPT;
            con_DoR_LPT = std::abs(0.4-DoR_LPT)-0.1;
            con_LPT_width = 0.009-(R_LPTit-R_LPTih);
            con_LPT_outlet_width = 0.009-A_LPTo/(2*M_PI*R_LPTom);
            const double sigma_max_LPTi = 0.5*omega_LPT*omega_LPT*rho_T*
                                         (R_LPTit*R_LPTit-R_LPTih*R_LPTih);
            con_sigma_max_LPTi = FOS_T*sigma_max_LPTi-sigma_max_T;
            // Use Q/(pi*omega*r^3), consistent with main's Cordier convention.
            phi_fan = (u_fani*A_fan/(M_PI*R_fan*R_fan)) / (omega_LPT*R_fan);
            psi_fan = (2*C_pc*D_T_fan) / (omega_LPT*omega_LPT*R_fan*R_fan);

            spec_speed_fan = std::pow(phi_fan,0.5)/std::pow(psi_fan,0.75);
            spec_dia_fan = std::pow(psi_fan,0.25)/std::pow(phi_fan,0.5);
            // The expression for the cordier line is based on the compressor cordier line at https://manual.cfturbo.com/en/index.html?cordier.html
            con_cordier_fan = is_cordier(spec_speed_fan, spec_dia_fan, 10.);// The compressor must be in range of the Cordier line


        }

        // Nozzle: station 6 for turbojet, station 7 after aft-fan LPT.
        double Ts_6 = T_6*std::pow(Ps_6/P_6,(gam_h-1)/gam_h);
        double a_6 = std::pow(gam_h*R*Ts_6,0.5);
        double u_6 = a_6*std::pow((2/(gam_h-1))*(std::pow((P_6/Ps_6),(gam_h-1)/gam_h)-1),0.5);
        double M_6 = u_6/a_6;
        double R_6 = std::pow((m_dot*(1+f))/(u_6*(Ps_0/(R*Ts_6)))/M_PI,0.5);
        double F_core = m_dot*((1+f)*u_6 - u_0);
        double F = F_core + F_fan;
        double con_F = 445.-F; // Thrust at least 100lbf
        // Mass model
        double m_Compressor_conical = std::abs((M_PI*(R_Cih*R_Cih*R_Cih - R_Coh*R_Coh*R_Coh)))*(rho_C/3);
        double m_Compressor_prismatic = (R_Cot - R_Coh)*M_PI*R_Coh*R_Coh
            * rho_C; // volume in m^3 -> mass in kg
        double m_Compressor = std::max(m_Compressor_conical,m_Compressor_prismatic);
        double m_Turbine_conical = std::abs((M_PI*(R_Tih*R_Tih*R_Tih - R_Toh*R_Toh*R_Toh)))*(rho_T/3);
        double m_Turbine_prismatic = (R_Tot - R_Toh)*M_PI*R_Tot*R_Tot
            * rho_T; // volume in m^3 -> mass in kg
        double m_Turbine = std::max(m_Turbine_conical, m_Turbine_prismatic);
        // Diffuser mass proportional to spec speed and component radius squared
        double m_Diffuser = (0.07146/(0.03*0.03))*(R_Com*R_Com)*(0.4/(spec_speed_C<0.4?spec_speed_C:0.4))*(0.4/(spec_speed_C<0.4?spec_speed_C:0.4));
        double m_NGV = (0.2885/(0.0315287*0.0315287))*(R_Tit*R_Tit);

        double m_Duct = 2*M_PI*R_Com*R_Com*0.001*2700;
        double m_Nozzle = 2*M_PI*R_Tit*R_Tot*0.001*8000;

        double len_Combustor = u_Tia * 0.0005 + 0.01;
        double m_Shaft = (len_Combustor+R_Com+0.02) * M_PI * 0.005 * 0.005 * 8000;
        double R_overall = (R_Com * (0.036/0.03) * (0.4/(spec_speed_C<0.4?spec_speed_C:0.4)));
        // Core casing must enclose compressor and HPT geometry, not just diffuser.
        R_overall = std::max({R_overall,R_Cit,R_Cot,R_Tit,R_Tot});
        double m_Wall = 2* M_PI * R_overall * (len_Combustor + 0.02) * 0.00065 * 8000;
        double m_Combustor = 2*M_PI*(0.5*R_overall + 0.8*R_overall)*len_Combustor*0.00065*8000;
        
        double m_LPT_conical = std::abs(M_PI*(std::pow(R_LPTih,3)-std::pow(R_LPToh,3)))*rho_T/3.;
        double m_LPT_prismatic = (R_LPTot-R_LPToh)*M_PI*R_LPTot*R_LPTot*rho_T;
        double m_LPT = std::max(m_LPT_conical,m_LPT_prismatic);
        double m_LPTNGV = (0.2885/(0.0315287*0.0315287))*R_LPTit*R_LPTit;
        double m_transition = M_PI*(std::abs(R_LPTit*R_LPTit - R_Tot*R_Tot)+std::abs(R_LPTih*R_LPTih - R_Toh*R_Toh))*1.414*0.0015*rho_T;
        double m_fan = M_PI*(R_fan*R_fan - R_LPTot*R_LPTot)*0.0015*rho_T;
        double m_nozzle = M_PI*(std::abs(R_6*R_6 - R_LPTot*R_LPTot))*1.414*0.0015*rho_T;

        double m_total = 1.12*(m_Compressor + m_Turbine+
                        m_Diffuser + m_NGV + m_Duct+
                        (aft_fan_ ? m_nozzle : m_Nozzle)+m_Shaft+m_Wall+m_Combustor
                        +(aft_fan_ ? m_LPT+m_LPTNGV+m_transition+m_fan : 0.));


        //Calculate objective and normalize
        double Isp = (F/(m_dot*f*g));
        double con_min_Isp = (2870.) - Isp;
        vector_double ret = {m_total,//-(35*(std::min(F/355.86,1.25))+25*std::min(1.50/(3600/Isp),1.25)+10*std::min(F/(g*m_total*10),1.25)), 
            // 0*con_beta_Ci_tip/70., 
            // 0*con_beta_Co/70., 
            // 0*con_beta_NGV/70., 
            // 0*con_beta_Tim/70., 
            // 0*con_beta_Tom/70., // All the beta constraints are zeroed, so that 
            // // designs are not rejected on the basis of high metal angle
            con_cordier_compressor, 
            con_Diff_C/0.55, 
            con_DoR/0.4, 
            con_M_Ci_tip, 
            con_M_Co_ROT, 
            con_M_Co_STAT, 
            con_R_Coh,
            hpt_mach.ngv_hub-0.9,
            hpt_mach.inlet_hub-0.8,
            hpt_mach.outlet_hub-0.9,
            con_phi_T, 
            con_psi_T, 
            con_T_width/0.01, // con_turbine_diffusion, 
            con_turbine_outlet_width/0.01,
            con_R_Toh/0.005,
            con_sigma_max_Ci/1.0e9,
            con_sigma_max_Co/1.0e9,
            con_sigma_max_Ti/1.0e9,
            con_sigma_max_To/1.0e9,
            con_F/500.,
            con_min_Isp/3000.
            };

        if (aft_fan_) {
            ret.insert(ret.end(), {
                con_LPT_backpres/101325., con_phi_LPT, con_psi_LPT,
                con_Mach_fan_tip, con_cordier_fan,
                con_sigma_max_fano/1.e9, con_sigma_max_LPTo/1.e9,
                (0.004-R_LPToh)/0.005,
                lpt_mach.outlet_hub-0.9,
                con_DoR_LPT/0.4, lpt_mach.ngv_hub-0.9, lpt_mach.inlet_hub-0.8,
                con_LPT_width/0.01, con_LPT_outlet_width/0.01,
                con_sigma_max_LPTi/1.e9
            });
        }

        // Preserve prior indices: former meanline slots now constrain the hub.
        // Append shroud constraints after the existing engine-specific entries.
        ret.insert(ret.end(), {hpt_mach.inlet_shroud-0.8,hpt_mach.outlet_shroud-0.9});
        if (aft_fan_)
            ret.insert(ret.end(), {lpt_mach.inlet_shroud-0.8,lpt_mach.outlet_shroud-0.9});

        // Append one geometry heuristic per turbine, preserving all prior slots.
        ret.push_back(con_HPT_migration);
        if (aft_fan_) ret.push_back(con_LPT_migration);

        if (invalid_ret(ret)) return penalty();

            if(results_ptr != nullptr){
                (*results_ptr)["[N] Total thrust"] = F;
                (*results_ptr)["[N] Core thrust"] = F_core;
                (*results_ptr)["[N] Fan thrust"] = F_fan;
                (*results_ptr)["[s] Specific impulse"] = Isp;
                (*results_ptr)["[kg] Total mass"] = m_total;
                (*results_ptr)["[kg/s] Core air mass flow"] = m_dot;
                (*results_ptr)["[kg/s] Fuel mass flow"] = m_dot*f;
                if (aft_fan_) {
                    (*results_ptr)["[m] LPT hub radius change (inlet-outlet)"] = R_LPTih-R_LPToh;
                    (*results_ptr)["[m] LPT shroud radius change (inlet-outlet)"] = R_LPTit-R_LPTot;
                    (*results_ptr)["[] LPT radius migration residual"] = con_LPT_migration;
                    (*results_ptr)["[m] LPT inlet hub radius"] = R_LPTih;
                    (*results_ptr)["[m] LPT inlet tip radius"] = R_LPTit;
                    (*results_ptr)["[m] LPT coefficient reference radius"] = R_LPT_mean;
                    (*results_ptr)["[m/s] LPT inlet axial velocity"] = u_LPTia;
                    (*results_ptr)["[m/s] LPT inlet swirl"] = u_LPTith_STAT;
                    (*results_ptr)["[m/s] LPT outlet axial velocity"] = u_LPToa;
                    (*results_ptr)["[K] LPT inlet static temperature"] = Ts_LPTi;
                    (*results_ptr)["[Pa] LPT inlet static pressure"] = Ps_LPTi;
                    (*results_ptr)["[K] LPT outlet static temperature"] = Ts_LPTo;
                    (*results_ptr)["[] LPT NGV exit Mach"] = M_LPTNGVo;
                    (*results_ptr)["[] LPT inlet relative Mach"] = M_LPTi;
                    (*results_ptr)["[] LPT outlet absolute Mach"] = M_LPToa;
                    (*results_ptr)["[deg] LPT NGV mean angle"] = beta_LPTNGV;
                    (*results_ptr)["[deg] LPT rotor inlet mean angle"] = beta_LPTim;
                    (*results_ptr)["[deg] LPT rotor outlet mean angle"] = beta_LPTom;
                    (*results_ptr)["[] LPT degree of reaction"] = DoR_LPT;
                    (*results_ptr)["[kg] LPT NGV mass estimate"] = m_LPTNGV;
                    evaluate_turbine_blading(*results_ptr,"LPT",omega_LPT,
                        R_LPTih,R_LPTit,R_LPToh,R_LPTot,u_LPTith_STAT,u_LPTia,u_LPToa,T_5,C_ph);
                    (*results_ptr)["[Pa] P_6 (LPT discharge)"] = P_6;
                    (*results_ptr)["[K] T_6 (LPT discharge)"] = T_6;
                    (*results_ptr)["[] M_7 (core nozzle outlet)"] = M_6;
                    (*results_ptr)["[kg/s] Fan air mass flow"] = m_dot_fan;
                    (*results_ptr)["[] Bypass ratio"] = m_dot_fan/m_dot;
                    (*results_ptr)["[K] Fan temperature rise"] = D_T_fan;
                    (*results_ptr)["[Pa] Fan discharge total pressure"] = P_fano;
                    (*results_ptr)["[K] Fan discharge total temperature"] = T_fano;
                    (*results_ptr)["[m] Fan tip radius"] = R_fan;
                    (*results_ptr)["[m] LPT outlet hub radius"] = R_LPToh;
                    (*results_ptr)["[] LPT flow coefficient"] = phi_LPT;
                    (*results_ptr)["[] LPT work coefficient"] = psi_LPT;
                    (*results_ptr)["[] Fan flow coefficient"] = phi_fan;
                    (*results_ptr)["[] Fan work coefficient"] = psi_fan;
                    (*results_ptr)["[] Fan specific speed"] = spec_speed_fan;
                    (*results_ptr)["[W] LPT extracted power"] = -m_dot*(1+f)*C_ph*D_T_LPT;
                    (*results_ptr)["[N m] LPT-to-fan transmitted torque"] = -m_dot*(1+f)*C_ph*D_T_LPT/omega_LPT;
                    (*results_ptr)["[W] Fan absorbed power"] = m_dot_fan*C_pc*D_T_fan;
                    (*results_ptr)["[] Fan relative tip Mach"] = con_Mach_fan_tip + M_max_fan;
                    (*results_ptr)["[] LPT exit relative Mach"] = M_LPTo; // Legacy meanline key
                    record_turbine_mach(*results_ptr,"LPT",lpt_mach,M_LPTi,M_LPTo,M_LPTNGVo);
                    (*results_ptr)["[] Fan efficiency (assumed)"] = eta_fan;
                    (*results_ptr)["[kg] LPT mass estimate"] = m_LPT;
                    (*results_ptr)["[kg] Fan mass estimate"] = m_fan;
                }
                (*results_ptr)["[m] HPT hub radius change (inlet-outlet)"] = R_Tih-R_Toh;
                (*results_ptr)["[m] HPT shroud radius change (inlet-outlet)"] = R_Tit-R_Tot;
                (*results_ptr)["[] HPT radius migration residual"] = con_HPT_migration;
                (*results_ptr)["[kg] Compressor prismatic mass"] = m_Compressor_prismatic;
                (*results_ptr)["[kg] HPT prismatic mass"] = m_Turbine_prismatic;
                (*results_ptr)["[m] Core casing radius"] = R_overall;
                (*results_ptr)["[] HPT degree of reaction"] = DoR;
                record_turbine_mach(*results_ptr,"HPT",hpt_mach,M_Ti,M_Tom,M_NGVo);
                evaluate_turbine_blading(*results_ptr,"HPT",omega,R_Tih,R_Tit,R_Toh,R_Tot,
                                         u_Tith_STAT,u_Tia,u_Toa,T_4,C_ph);
                (*results_ptr)["[] Compressor specific speed"] = spec_speed_C;
                (*results_ptr)["[Pa] P_0 (ambient)"] = P_0;
                (*results_ptr)["[K] T_0 (ambient)"] = T_0;
                (*results_ptr)["[Pa] P_3 (compressor discharge)"] = P_3;
                (*results_ptr)["[K] T_3 (compressor discharge)"] = T_3;
                (*results_ptr)["[Pa] P_4 (combustor discharge)"] = P_4;
                (*results_ptr)["[K] T_4 (combustor discharge)"] = T_4;
                (*results_ptr)["[Pa] P_5 (turbine discharge)"] = P_5;
                (*results_ptr)["[K] T_5 (turbine discharge)"] = T_5;
                if (!aft_fan_) (*results_ptr)["[] M_6 (nozzle outlet)"] = M_6;
                (*results_ptr)["[] Compressor flow coefficient"] = phi_C;
                (*results_ptr)["[] Turbine flow coefficient"] = phi_T;
                (*results_ptr)["[] Compressor work coefficient"] = psi_C;
                (*results_ptr)["[] Turbine work coefficient"] = psi_T;
            }
        return ret;
    }catch(const boost::math::evaluation_error&){
        // If the compute_u_a throws an error, return high penalty immediately
        vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),1e+6);
        return ret;
    }
}
// Calculates if a specific speed and specific work lie on the Cordier line
double problem_jet_calc::is_cordier(double sigma, double delta, double max_sigma) const{
    double ideal_delta{};
    double min{};
    double max{};

    if(!std::isfinite(sigma) || !std::isfinite(delta) || sigma<0.05 || sigma>max_sigma){
        return 1.;
    }
    if(sigma>=0.05 && sigma < 0.4){
        ideal_delta = 1.1173*std::pow(sigma, -0.963);
    }
    if(sigma>=0.4 && sigma < 0.8){
        ideal_delta = 1.46477*std::pow(sigma, -0.66742);
    }
    if(sigma>=0.8 && sigma <= max_sigma){
        ideal_delta = 1.61299*std::pow(sigma, -0.23543);
    }
    min = 0.9 * ideal_delta;
    max = 1.1 * ideal_delta;

    return 4*((delta-min)*(delta-max))/(std::pow(max-min,2));
}

double problem_jet_calc::cordier(double sigma) const{
    if(sigma<0.05 || sigma>2.5){
        return 0;
    }
    if(sigma>=0.05 && sigma < 0.4){
        return 1.1173*std::pow(sigma, -0.963);
    }
    if(sigma>=0.4 && sigma < 0.8){
        return 1.46477*std::pow(sigma, -0.66742);
    }
    if(sigma>=0.8 && sigma <= 2.5){
        return 1.61299*std::pow(sigma, -0.23543);
    }
    return 0;
}

// Checks if anything in a vector_double is inf or nan
bool problem_jet_calc::invalid_ret(vector_double &x) const{
    for (double xi : x){
        if(std::isinf(xi) || std::isnan(xi)){
            return true;
        }
    }
    return false;
}

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
    if (!(std::isfinite(m_dot) && std::isfinite(rho_t) && std::isfinite(T_t) &&
          std::isfinite(u_th) && std::isfinite(gam) && std::isfinite(A) &&
          m_dot > 0. && rho_t > 0. && T_t > 0. && gam > 1. && A > 0. &&
          2.*R*(gam/(gam-1.))*T_t > u_th*u_th))
        return std::numeric_limits<double>::quiet_NaN();
    // For now, if an evaluation_error ever comes up in optimizer runs, try to find a way to handle this without relying on catching
    // try{
        double C_p = this->R * (gam/(gam-1));
        // Positive point at which the derivative is zero, the subsonic solution
        // will always be less than this so it is used as the upper bound
        double u_a_max = std::pow(((gam-1)*rho_t*(2*C_p*T_t-(u_th*u_th)))/((1+gam)*rho_t),0.5);
        if (!(std::isfinite(u_a_max) && u_a_max > 0.0001))
            return std::numeric_limits<double>::quiet_NaN();
        const double Ts_limit = T_t-(u_a_max*u_a_max+u_th*u_th)/(2*C_p);
        const double max_flow = u_a_max*A*rho_t*std::pow(Ts_limit/T_t,1./(gam-1.));
        if (!std::isfinite(max_flow) || m_dot > max_flow)
            return std::numeric_limits<double>::quiet_NaN();
        auto compute_u_a_residual = [=](double u_a)->std::tuple<double,double>{
            double T_s = T_t - ((std::pow(u_a,2)+std::pow(u_th,2))/(2*C_p));
            double tau = T_t/T_s;
            double u_a_residual = m_dot - u_a*A*rho_t*std::pow(tau,(1/(1-gam)));
            double D_u_a_residual = (A*rho_t*(2*C_p*(-1+gam)*T_t-(1+gam)*std::pow(u_a,2)-(-1+gam)*std::pow(u_th,2))*
                                    std::pow(T_t/(T_t-(std::pow(u_a,2)+std::pow(u_th,2))/(2.*C_p)),1/(1-gam)))/
                                    ((-1+gam)*(-2*C_p*T_t+std::pow(u_a,2)+std::pow(u_th,2)));
            return {u_a_residual,D_u_a_residual};
        };
        int digits = 32;
        std::uintmax_t max_iter = 50;
        double u_a = boost::math::tools::newton_raphson_iterate(compute_u_a_residual, 0., 0.0001, u_a_max, digits, max_iter);
        // Error checking that 1. u_a is in fact a root and 2. u_a is not inf or nan
        const double residual = std::get<0>(compute_u_a_residual(u_a));

        if (!std::isfinite(u_a) || !std::isfinite(residual)
            || std::abs(residual) > 1e-6) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return u_a;
    // }catch(boost::math::evaluation_error){
        // return static_cast<double>(NAN);
    // }
}