#include <matplot/matplot.h>
#include <cmath>
#include "jet_calc.h"


std::vector<double> problem_jet_calc::compute_contour_val(double T_4, double OPR, double spec_speed_C) const{
    std::vector<double> x = {11475.8, 135.212, 1097.15, 0.00501128, 0.0206151, 0.00124931, 0.0299996, 62.201, 0.0154871, 0.0315287, 0.00210709, 0.0188079, };
    double des_psi_C = 1/std::pow(cordier(spec_speed_C)*spec_speed_C,2);
    double des_phi_C = spec_speed_C*spec_speed_C*std::pow(des_psi_C,3/2);
    // TODO: find x from x_opt and args
    auto R_Com = 0.03; // m     - compressor outlet meanline radius
    auto u_i = x[1];   // m/s   - compressor inlet velocity
    auto D_T_C = ((Ts_0 + ((u_i*u_i)/(2*C_pc)))/eta_C)*(std::pow(OPR,((gam_c-1)/gam_c))-1); // K     - compressor total temperature change
    auto omega = std::pow(((2*D_T_C*C_pc)/(R_Com*R_Com*des_psi_C)),0.5); // rad/s - shaft speed
    

    // auto T_4 = T_4;   // K     - combustor exit total temp
    auto R_Cih = 0.0; // m     - compressor inlet hub radius
    auto R_Cit = std::pow((des_phi_C*omega*R_Com*R_Com*R_Com)/(u_i),0.5); // m     - compressor inlet tip radius
    auto A_Co = x[5];  // m^2   - compressor outlet area

    auto R_Tih = x[8]; // m     - turbine inlet hub radius
    auto R_Tit = x[9]; // m     - Turbine inlet tip radius
    auto A_To = x[10]; // m^2   - Turbine outlet area
    auto R_Tom = x[11];// m     - Turbine exit meanline velocity
    // Variable subscript convention:
    // if u:
    // u_<Component : n><Component axial location : 1><Velocity direction : 1 - 2>_<component radial location : n>
    // otherwise:
    // <var symbol : n>_<Component : n><component axial location : 1><component radial location : n>
    pagmo::vector_double res;
    double F;
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

        // // Compressor constraints
        // // Velocities are positive in the direction of rotation
        // // phi and psi are calculated according to https://manual.cfturbo.com/en/index.html?md_parameters_axvent.html
        // double phi_C = (u_i*(A_Ci/(M_PI*R_Com*R_Com))) / (omega*R_Com);
        // double psi_C = (2*C_pc*D_T_C) / (omega*omega*R_Com*R_Com);
        // double spec_speed_C = std::pow(phi_C,0.5)/std::pow(psi_C,0.75);
        // double spec_dia_C = std::pow(psi_C,0.25)/std::pow(phi_C,0.5);
        // // The expression for the cordier line is based on the compressor cordier line at https://manual.cfturbo.com/en/index.html?cordier.html
        // double con_cordier_compressor = is_cordier(spec_speed_C, spec_dia_C)?-1:1;// The compressor must be in range of the Cordier line
        double M_Ci_tip = std::pow((u_i*u_i + omega*omega*R_Cit*R_Cit)/(gam_c*R*Ts_2),0.5);
        double con_M_Ci_tip = M_Ci_tip - 0.8;// Mach at compressor inlet less than 0.8
        // double beta_Ci_tip = (180/M_PI)*std::atan((omega*R_Cit)/u_i);
        // double con_beta_Ci_tip = beta_Ci_tip - 70;// Inlet blade angle less than 70 degrees
        // double u_Coth_STAT = (C_pc*D_T_C)/(omega*R_Com);
        // double u_Coth_ROT = u_Coth_STAT - (omega*R_Com);
        // double u_Coa = compute_u_a(m_dot,(P_3/(R*T_3)), T_3, u_Coth_STAT, gam_c, A_Co);
        // if(std::isnan(u_Coa)){
        //     vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),0.0);
        //     return ret;
        // }
        // double u_Co_ROT = std::pow((u_Coth_STAT-(omega*R_Com))*(u_Coth_STAT-(omega*R_Com)) + u_Coa*u_Coa,0.5);
        // double u_Co_STAT = std::pow(u_Coth_STAT*u_Coth_STAT + u_Coa*u_Coa,0.5);
        // double u_Ci_ROT = std::pow((u_i*u_i + omega*omega*R_Cit*R_Cit),0.5);
        // double Diff_C = 1 - (u_Co_ROT/u_Ci_ROT)+u_Coth_STAT/(2*sigma_C*u_Ci_ROT);
        // double con_Diff_C = Diff_C - 0.55;// Lieblein diffusion factor less than 0.55
        // // u_Coth_ROT is negative because the angle is defined as positive in the direction opposing rotation
        // double beta_Co = (180/M_PI)*std::atan(-u_Coth_ROT/u_Coa);
        // double con_beta_Co = beta_Co - 60; // Compressor blade exit angle below 40 degrees
        // double Ts_3 = T_3 - (u_Co_STAT*u_Co_STAT)/(2*C_pc);
        // double a_3 = std::pow(gam_c*R*Ts_3, 0.5);
        // double M_Co_ROT = u_Co_ROT/a_3;
        // double con_M_Co_ROT = M_Co_ROT - 0.8;// Compressor exit relative Mach less than 0.8
        // double M_Co_STAT = std::pow((u_Coa*u_Coa + u_Coth_STAT*u_Coth_STAT)/(gam_c*R*Ts_3),0.5);
        // double con_M_Co_STAT = M_Co_STAT - 0.8;// Compressor diffuser inlet Mach less than 0.8
        // double sigma_max_Ci = 0.5*omega*omega*rho_C*(R_Cit*R_Cit - R_Cih*R_Cih);
        // double con_sigma_max_Ci = FOS_C * sigma_max_Ci - sigma_max_C;// FOS at compressor inlet blade root
        // double R_Coh = R_Com - (A_Co/(4*M_PI*R_Com));
        // double R_Cot = R_Com + (A_Co/(4*M_PI*R_Com));
        // double sigma_max_Co = 0.5*omega*omega*rho_C*(R_Cot*R_Cot - R_Coh*R_Coh);
        // double con_sigma_max_Co = FOS_C * sigma_max_Co - sigma_max_C;// FOS at compressor outlet blade root

        // #ifdef EVAL_JET_CALC
        // double beta_Di = (180/M_PI)*std::atan(u_Coth_STAT / u_Coa);

        // double R_Co_x = R_Com + 0.003;
        // double u_Coth_x_STAT = u_Coth_STAT * (R_Com/R_Co_x);
        // double u_Coth_x_ROT = u_Coth_STAT - omega*R_Co_x;
        // double u_Coa_x = compute_u_a(m_dot,(P_3/(R*T_3)), T_3, u_Coth_x_STAT, gam_c, A_Co);
        // double beta_Co_x = (180/M_PI)*std::atan(-u_Coth_x_ROT/u_Coa_x);
        // double beta_Di_x = (180/M_PI)*std::atan(u_Coth_x_STAT / u_Coa_x);
        // #endif

        // Combustor
        double f = (C_ph*T_4 - C_pc*T_3)/(h_ker - C_ph*T_4);
        // Turbine
        double D_T_T = -(P_spC/(C_ph*(1+f)));
        double T_5 = T_4+D_T_T;
        double P_5 = P_3*std::pow((T_4+(D_T_T/eta_T))/T_4,(gam_h/(gam_h-1)));
        // // Turbine constraints
        // double A_Ti = M_PI*(R_Tit*R_Tit - R_Tih*R_Tih);
        // double R_Tim = 0.5*(R_Tit+R_Tih);
        // double u_Tith_STAT = (C_ph*(-1*D_T_T))/(omega*R_Tim);// Since D_T_T is negative, it has to be flipped to calculate the absolute NGV exit u_th
        // double u_Tia = compute_u_a(m_dot*(1+f),P_3/(R*T_4), T_4, u_Tith_STAT, gam_h, A_Ti);
        // if(std::isnan(u_Tia)){
        //     vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),0);
        //     return ret;
        // }
        // double beta_NGV = (180/M_PI)*std::atan(u_Tith_STAT/u_Tia);
        // double con_beta_NGV = beta_NGV - 70;
        // // Phi and psi follow the definitions in Dixon and Hall
        // // It is calculated with reference to the turbine tip inlet speed. This is apparently the 
        // // standard but might apply only to axial turbines where blade speed doesn't vary much axially
        // double phi_T = u_Tia / (omega*R_Tit);
        // double psi_T = (C_ph*(-D_T_T)) / (omega*omega*R_Tom*R_Tom);//D_T_T is flipped here also to follow work coefficient convention
        // // Phi and psi constraints are for the smith chart for axial gas turbines in dixon and hall.
        // // Smith believed that the losses were proportional to the average kinetic energy in the row,
        // // and correlated this with empirically measured losses.
        // // If this is true, these phi and psi constraints should hold for any topology.
        // double con_phi_T = std::abs(phi_T - (0.5*(min_phi_T+max_phi_T))) - (0.5*(max_phi_T - min_phi_T));// Flow coefficient in appropriate range
        // double con_psi_T = std::abs(psi_T - (0.5*(min_psi_T+max_psi_T))) - (0.5*(max_psi_T - min_psi_T));// Loading coefficient in appropriate range
        // double u_Ti_STAT = std::pow(u_Tith_STAT*u_Tith_STAT+u_Tia*u_Tia,0.5);
        // double D_Ts_NGV = -1*(u_Ti_STAT*u_Ti_STAT) / (2*C_ph);
        // double Ts_Ti = T_4 - D_Ts_NGV;
        // double DoR = 1 - (D_Ts_NGV/D_T_T);
        // double con_DoR = std::abs(0.4 - DoR) - 0.1; // Degree of reaction of turbine between 0.3 and 0.5
        // double M_NGVo = std::pow((u_Tith_STAT*u_Tith_STAT+u_Tia*u_Tia)/(gam_h*R*Ts_Ti),0.5);
        // double con_M_NGVo = M_NGVo - 0.8; // Mach in NGV exit less than 0.8
        // double u_Tith_ROT = u_Tith_STAT - omega*R_Tim;
        // double M_Ti = std::pow((u_Tith_ROT*u_Tith_ROT + u_Tia*u_Tia)/(gam_h*R*Ts_Ti),0.5);
        // double con_M_Ti = M_Ti - 0.8;// Turbine inlet mach less than 0.8
        // double beta_Tim = (180/M_PI)*std::atan(u_Tith_ROT/u_Tia);
        // double con_beta_Tim = std::abs(beta_Tim) - 65;// Turbine inlet angle less than 65 degrees
        // double u_Toa = compute_u_a(m_dot * (1+f), P_5/(R*T_5), T_5, 0, gam_h, A_To);
        // if(std::isnan(u_Toa)){
        //     vector_double ret(1+static_cast<int>(get_nec())+static_cast<int>(get_nic()),0);
        //     return ret;
        // }
        // double beta_Tom = (180/M_PI)*std::atan((omega*R_Tom)/u_Toa);
        // double con_beta_Tom = beta_Tom - 65; // Turbine outlet meridional blade angle less than 65 degrees - this shouldn't be active
        // double con_T_width = 0.009 - (R_Tit - R_Tih);// Turbine inlet annulus width greater than 9mm
        // // double con_turbine_diffusion = u_Tia - u_Toa; // Flow must accelerate through turbine to avoid separation
        // double Ts_To = T_5 - (u_Toa*u_Toa)/(2*C_ph);
        // double con_turbine_outlet_width = 0.009 - A_To/(2*M_PI*R_Tom); // Turbine outlet width more than 9mm
        // double R_Tot = R_Tom + A_To/(4*M_PI*R_Tom);
        // double M_Tom = std::pow((u_Toa*u_Toa + omega*omega*R_Tot*R_Tot)/(gam_h*R*Ts_To),0.5);
        // double con_M_Tom = M_Tom - 0.8; // Relative Mach at turbine exit less than 0.8
        // double sigma_max_Ti = 0.5*omega*omega*rho_C*(R_Tit*R_Tit - R_Tih*R_Tih);
        // double con_sigma_max_Ti = FOS_T * sigma_max_Ti - sigma_max_T;// FOS at compressor inlet blade root
        // double R_Toh = R_Tom - A_To/(4*M_PI*R_Tom);
        // double con_R_Toh = 0.004 - R_Toh;// Inner radius of turbine outlet must be at least 4mm
        // double sigma_max_To = 0.5*omega*omega*rho_C*(R_Tot*R_Tot - R_Toh*R_Toh);
        // double con_sigma_max_To = FOS_T * sigma_max_To - sigma_max_T;// FOS at compressor inlet blade root
        // double con_turbine_tip_geom = R_Tit - R_Tot;
        // double con_turbine_hub_geom = R_Toh - R_Tih;

        // #ifdef EVAL_JET_CALC
        // double u_NGVoth_tip = (R_Tim / R_Tit) * u_Tith_STAT;// NGV outlet theta tip
        // double u_NGVoa_tip = compute_u_a(m_dot*(1+f),P_3/(R*T_4), T_4, u_NGVoth_tip, gam_h, A_Ti);
        // double beta_NGVot = (180/M_PI)*std::atan(u_NGVoth_tip / u_NGVoa_tip);
        // double u_Tith_tip_ROT = u_NGVoth_tip - omega*R_Tit;
        // double beta_Tit = (180/M_PI)*std::atan((u_Tith_tip_ROT)/u_NGVoa_tip);
        // double u_Toa_tip = u_Toa;
        // double beta_Tot = (180/M_PI)*std::atan((omega*R_Tot)/u_Toa_tip);


        // double u_NGVoth_hub = (R_Tim / R_Tih) * u_Tith_STAT;// NGV outlet theta hub
        // double u_NGVoa_hub = compute_u_a(m_dot*(1+f),P_3/(R*T_4), T_4, u_NGVoth_hub, gam_h, A_Ti);
        // double beta_NGVoh = (180/M_PI)*std::atan(u_NGVoth_hub / u_NGVoa_hub);
        // double u_Tith_hub_ROT = u_NGVoth_hub - omega*R_Tih;
        // double beta_Tih = (180/M_PI)*std::atan((u_Tith_hub_ROT)/u_NGVoa_hub);
        // double u_Toa_hub = u_Toa;
        // double beta_Toh = (180/M_PI)*std::atan((omega*R_Toh)/u_Toa_hub);
        // #endif

        // Nozzle
        double Ts_6 = T_5*std::pow(Ps_7/P_5,(gam_h-1)/gam_h);
        double a_6 = std::pow(gam_h*R*Ts_6,0.5);
        double u_6 = a_6*std::pow((2/(gam_h-1))*(std::pow((P_5/Ps_7),(gam_h-1)/gam_h)-1),0.5);
        #ifdef EVAL_JET_CALC
        double R_6 = std::pow((m_dot*(1+f))/(u_6*(Ps_0/(R*Ts_6)))/M_PI,0.5);
        #endif
        double F = m_dot*((1+f)*u_6 - u_0);
        double con_F = 70-F; // Thrust at least 70N
        //Calculate objective
        double Isp = (F/(m_dot*f*9.8066));
        double eta_thermal = ((1+f)*u_6*u_6)/(2*f*h_ker);
        double inv_W_sp = std::pow(u_6*u_6*0.5,-0.5);
        double W_sp = 0.5*(u_6*u_6 - u_0*u_0);
        res = {Isp, 
            eta_thermal,
            W_sp,
            // con_beta_Ci_tip, 
            // con_beta_Co, 
            // con_beta_NGV, 
            // con_beta_Tim, 
            // con_beta_Tom, 
            // con_cordier_compressor, 
            // con_Diff_C, 
            // con_DoR, 
            con_M_Ci_tip, 
            // con_M_Co_ROT, 
            // con_M_Co_STAT, 
            // con_M_NGVo, 
            // con_M_Ti, 
            // con_M_Tom, 
            // con_phi_T, 
            // con_psi_T, 
            // con_T_width, // con_turbine_diffusion, 
            // con_turbine_outlet_width,
            // con_R_Toh,
            // con_sigma_max_Ci,
            // con_sigma_max_Co,
            // con_sigma_max_Ti,
            // con_sigma_max_To,
            F
            };
        // A physically impossible engine will usually result in a bunch of NaNs,
        // and bad inputs might give Inf due to division by zero.
        // If this happens, return a big penalty
        if (invalid_ret(res)){
            res = vector_double(6,0.0);
        }
    }catch(boost::math::evaluation_error){
        // If the compute_u_a throws an error, return high penalty immediately
        res = vector_double(6,0.0);
    }
    // Currently, this returns 6 values: imperial SFC, thermal effy, prop effy, inlet mach constraint, and a fake constraint that checks that the inlet is less than the outlet, and finally thrust
    std::vector<double> ret = {res[0], res[1], res[2], res[3], (R_Cit>0.03)?1.0:-1.0, res[4]};
    return ret;
}

void print_matrix_numpy_style(const matplot::vector_2d& matrix, int precision = 4, int width = 10) {
    std::cout << "[\n";
    for (const auto& row : matrix) {
        std::cout << "  [";
        for (size_t j = 0; j < row.size(); ++j) {
            std::cout << std::fixed << std::setw(width) << std::setprecision(precision) << row[j];
            if (j != row.size() - 1)
                std::cout << ", ";
        }
        std::cout << "],\n";
    }
    std::cout << "]\n";
}

int main(){
    using namespace matplot;

    problem_jet_calc pjc;
    // auto test = pjc.compute_contour_val(1094,1.323);
    int num_points = 100;
    int num_values = 6;

    vector_1d T_4_1d = linspace(700,1600,num_points);
    vector_1d OPR_1d = linspace(1,6,num_points);
    auto [T_4_2d, OPR_2d] = meshgrid(T_4_1d, OPR_1d);
    // (5, num_points, num_points) array of data
    std::vector<std::vector<std::vector<double>>> data(num_values, std::vector<std::vector<double>>(num_points,std::vector<double>(num_points)));
    // Loop through the grid, for every point compute the data, then store the data along the correct slice
    for(int i = 0; i<num_points; i++){
        for(int j = 0; j<num_points; j++){
            auto temp_data = pjc.compute_contour_val(T_4_1d[i],OPR_1d[j],0.4);
            for(int k=0; k<num_values; k++){
                // The indices are k,j,i here because contourf expects j,i arrays for plotting
                // Dimension j is for the y axis, i is for x axis
                data[k][j][i] = temp_data[k];
            }
        }
    }
    int max_OPR_idx = 0;
    // Start at 1 here because the first row is with an OPR of 1, so it is not physical and shouldn't be considered
    for (int i=1;i<num_points;i++){
        double temp_con_M_Cit = data[3][i][0];
        if(temp_con_M_Cit < 0){
            continue;
        }else{
            max_OPR_idx = i-1;
            break;
        }
    }
    double max_F = data[5][max_OPR_idx][num_points-1];
    double max_ISP = 0;
    double T_for_max_Isp = 0;
    for(int i=0;i<num_points;i++){
        double temp_ISP = data[0][max_OPR_idx][i];
        if (temp_ISP > max_ISP){
            max_ISP = data[0][max_OPR_idx][i];
            T_for_max_Isp = T_4_1d[i];
        }
    }

    // print_matrix_numpy_style(data[5]);
    // print_matrix_numpy_style(data[1]);
    auto f1 = figure(true);
    f1->size(1920, 1080); 
    auto ax1 = f1->add_subplot(2,2,1);
    ax1->hold(on);
    // Plot data
    auto ISP = data[0];
    auto contour_1_1 = ax1->contour(T_4_2d,OPR_2d,ISP);
    contour_1_1->line_width(2.0);
    // Plot and style constraints
    auto isocline_1_1 = ax1->contour(T_4_2d,OPR_2d,data[3], vector_1d({0.0}));
    isocline_1_1->line_width(3.0);
    isocline_1_1->color("red");
    auto isocline_1_2 = ax1->contour(T_4_2d,OPR_2d,data[4], vector_1d({0.0}));
    isocline_1_2->line_width(3.0);
    isocline_1_2->color("blue");
    // Add fripperies
    ax1->colormap();
    ax1->title("Specific Impulse, s");

    auto ax2 = f1->add_subplot(2,2,2);
    ax2->hold(on);
    // Plot data
    auto ETA_TH = data[1];
    auto contour_2_1 = ax2->contour(T_4_2d,OPR_2d,ETA_TH);
    contour_2_1 -> line_width(2.0);
    // Plot and style constraints
    auto isocline_2_1 = ax2->contour(T_4_2d,OPR_2d,data[3], vector_1d({0.0}));
    isocline_2_1->line_width(3.0);
    isocline_2_1->color("red");
    auto isocline_2_2 = ax2->contour(T_4_2d,OPR_2d,data[4], vector_1d({0.0}));
    isocline_2_2->line_width(3.0);
    isocline_2_2->color("blue");
    // Add fripperies
    ax2->colormap();
    ax2->title("Thermal Efficiency");

    auto ax3 = f1->add_subplot(2,2,3);
    ax3->hold(on);
    // Plot data
    auto INV_W_SP = data[2];
    auto contour_3_1 = ax3->contour(T_4_2d,OPR_2d,INV_W_SP);
    contour_3_1->line_width(2.0);
    // Plot and style constraints
    auto isocline_3_1 = ax3->contour(T_4_2d,OPR_2d,data[3], vector_1d({0.0}));
    isocline_3_1->line_width(3.0);
    isocline_3_1->color("red");
    auto isocline_3_2 = ax3->contour(T_4_2d,OPR_2d,data[4], vector_1d({0.0}));
    isocline_3_2->line_width(3.0);
    isocline_3_2->color("blue");
    // Add fripperies
    ax3->colormap();
    ax3->title("Specific Work");

    auto ax4 = f1->add_subplot(2,2,4);
    ax4->hold(on);
    // Plot data
    auto THRUST = data[5];
    auto contour_4_1 = ax4->contour(T_4_2d,OPR_2d,THRUST);
    contour_4_1->line_width(2.0);
    // Plot and style constraints
    auto isocline_4_1 = ax4->contour(T_4_2d,OPR_2d,data[3], vector_1d({0.0}));
    isocline_4_1->line_width(3.0);
    isocline_4_1->color("red");
    auto isocline_4_2 = ax4->contour(T_4_2d,OPR_2d,data[4], vector_1d({0.0}));
    isocline_4_2->line_width(3.0);
    isocline_4_2->color("blue");
    // Add fripperies
    ax4->colormap();
    ax4->title("Thrust");


    sleep(1);
    // Plot the entire figure
    f1->show();
    
    // vector_1d x = linspace(-10, 10, 100);
    // vector_1d y = linspace(-20, 20, 100);
    // auto [X,Y] = meshgrid(x,y);
    // vector_2d Z = transform(X,Y,[](double x, double y){return x*x + y*y;});

    // auto f1 = figure(true);
    // auto ax1 = f1->add_subplot(1,1,0);

    // ax1->contour(X,Y,Z);
    // ax1->colormap();
    // f1->show();

    // auto f2 = figure(true);
    // auto ax2 = f2->add_subplot(1,1,0);
    // ax2->hold(on);
    // ax2->contourf(X,Y,Z);
    // ax2->colormap();
    // auto isoclines = ax2->contourf(X,Y,Z);
    // isoclines->color("white");
    // isoclines->line_width(1.5);
    // f2->show();
    return 0;
}