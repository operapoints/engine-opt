#ifndef JET_CALC_H
#define JET_CALC_H

#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>


#include <boost/math/tools/roots.hpp>

using namespace pagmo;

struct problem_jet_calc{
    public:
        // Non-varying parameters
        // Define things like environment parameters and compenent efficiencies here
        static constexpr const double R = 287.0;
        static constexpr double gam_c = 1.4;
        static constexpr double gam_h = 1.36;
        static constexpr double C_pc = R*(gam_c/(gam_c-1));
        static constexpr double C_ph = R*(gam_h/(gam_h-1));

        static constexpr double Ts_0 = 298.;
        static constexpr double Ps_0 = 101325.;
        static constexpr double Ps_6 = Ps_0;
        static constexpr double u_0 = 0.;

        static constexpr double h_ker = 43e+6;

        static constexpr double eta_C = 0.75;
        static constexpr double eta_T = 0.75;

        static constexpr double sigma_C = 1.6;

        static constexpr double min_phi_T = 0.4;
        static constexpr double max_phi_T = 0.8;
        static constexpr double min_psi_T = 1.0;
        static constexpr double max_psi_T = 2.0;

        static constexpr double rho_T = 8190.;// Inconel 718
        static constexpr double rho_C = 4450.;// Ti6Al4V
        static constexpr double sigma_max_T = 600000000.;// 600MPa for Inconel 718 at 1100K, source: https://www.vdm-metals.com/fileadmin/user_upload/Downloads/Data_Sheets/Data_Sheet_VDM_Alloy_718.pdf
        static constexpr double sigma_max_C = 880000000.;// 880MPa for Ti6Al4V
        static constexpr double FOS_T = 2;
        static constexpr double FOS_C = 2;
        vector_double::size_type get_nec() const;
        vector_double::size_type get_nic() const;
        vector_double fitness(vector_double &x) const;
        std::pair<vector_double, vector_double> get_bounds() const;
    // private: 
        double compute_u_a(double m_dot,
                        double rho_t,
                        double T_t,
                        double u_th,
                        double gam,
                        double A) const;
        inline bool invalid_ret(vector_double& x) const;
        inline bool is_cordier(double sigma, double delta) const;
};

#endif