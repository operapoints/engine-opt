#include "study_physics.hpp"

#include <cmath>
#include <stdexcept>

StagePerformance study_physics(const RotorVariables& rv, double gamma,
                               Rostack_diagnostics performance_params)
{
    if (!std::isfinite(gamma) || gamma <= 1.0 ||
        !std::isfinite(rv.sigma) || rv.sigma <= 0.0 ||
        !std::isfinite(rv.nu) || rv.nu <= 0.0 || rv.nu > 1.0) {
        throw std::domain_error("Invalid gamma, solidity, or hub-to-tip ratio");
    }

    const double cth1 = rv.cth1;
    const double cth2 = rv.cth2;
    const double cm1 = rv.cm1;
    const double cm2 = rv.k * cm1;
    const double u1 = rv.u1;
    const double u2 = rv.lambda * u1;
    const double nu = rv.nu;
    const double sigma = rv.sigma;
    const double tau = 1.0 + u2 * cth2 - u1 * cth1;

    const double c1sq = cm1 * cm1 + cth1 * cth1;
    const double c2sq = cm2 * cm2 + cth2 * cth2;
    const double wth1 = cth1 - u1;
    const double wth2 = cth2 - u2;
    const double w1sq = cm1 * cm1 + wth1 * wth1;
    const double w2sq = cm2 * cm2 + wth2 * wth2;

    // Static temperatures relative to the row-inlet TOTAL temperature.
    const double theta1 = 1.0 - c1sq / 2.0;
    const double theta2 = tau - c2sq / 2.0;
    if (!std::isfinite(tau) || tau <= 0.0 ||
        !std::isfinite(theta1) || theta1 <= 0.0 ||
        !std::isfinite(theta2) || theta2 <= 0.0 ||
        !std::isfinite(w1sq) || w1sq <= 0.0 ||
        !std::isfinite(w2sq)) {
        throw std::domain_error("Invalid tip thermodynamic state");
    }

    const double Mw1 = std::sqrt(w1sq / ((gamma - 1.0) * theta1));
    const double Mw2 = std::sqrt(w2sq / ((gamma - 1.0) * theta2));
    const double M2 = std::sqrt(c2sq / ((gamma - 1.0) * theta2));
    const double w1 = std::sqrt(w1sq);
    const double w2 = std::sqrt(w2sq);
    const double DF_tip = 1.0 - w2 / w1
        + std::abs(wth1 - wth2) / (2.0 * sigma * w1);

    // Free vortex at each station, with the same nu at inlet and outlet.
    const double cth1_hub = cth1 / nu;
    const double cth2_hub = cth2 / nu;
    const double u1_hub = nu * u1;
    const double u2_hub = nu * u2;
    const double wth1_hub = cth1_hub - u1_hub;
    const double wth2_hub = cth2_hub - u2_hub;
    const double w1_hub = std::hypot(cm1, wth1_hub);
    const double w2_hub = std::hypot(cm2, wth2_hub);
    const double sigma_hub = sigma / nu; // Constant spanwise chord.

    // u_hub * ctheta_hub == u_tip * ctheta_tip: Euler work is unchanged.
    const double theta1_hub = 1.0 - (cm1 * cm1 + cth1_hub * cth1_hub) / 2.0;
    const double theta2_hub = tau - (cm2 * cm2 + cth2_hub * cth2_hub) / 2.0;
    if (!std::isfinite(theta1_hub) || theta1_hub <= 0.0 ||
        !std::isfinite(theta2_hub) || theta2_hub <= 0.0 ||
        !std::isfinite(w1_hub) || w1_hub <= 0.0 ||
        !std::isfinite(w2_hub) ||
        !std::isfinite(sigma_hub) || sigma_hub <= 0.0) {
        throw std::domain_error("Invalid hub thermodynamic state");
    }

    const double DF = 1.0 - w2_hub / w1_hub
        + std::abs(wth1_hub - wth2_hub) / (2.0 * sigma_hub * w1_hub);
    if (!std::isfinite(Mw1) || !std::isfinite(Mw2) || !std::isfinite(M2) ||
        !std::isfinite(DF) || !std::isfinite(DF_tip)) {
        throw std::domain_error("Invalid Mach number or diffusion factor");
    }

    performance_params["M1 absolute"] = std::sqrt(c1sq / ((gamma - 1.0) * theta1));
    performance_params["inlet static T / row inlet total T"] = theta1;
    performance_params["outlet static T / row inlet total T"] = theta2;
    performance_params["inlet u1 (local)"] = u1;
    performance_params["outlet u2 (local)"] = u2;
    performance_params["diffusion factor, tip"] = DF_tip;
    performance_params["diffusion factor, hub"] = DF;
    performance_params["hub solidity"] = sigma_hub;
    performance_params["hub inlet static T / row inlet total T"] = theta1_hub;
    performance_params["hub outlet static T / row inlet total T"] = theta2_hub;

    return StagePerformance{tau, Mw1, Mw2, M2, std::max(DF,DF_tip)};
}
