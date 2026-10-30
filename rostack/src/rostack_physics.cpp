#include "rostack.hpp"

#include <cmath>
#include <stdexcept>

// Physics models

StagePerformance rostack_physics_simple(const RotorVariables& rv, double gamma){
    if (!std::isfinite(gamma) || gamma <= 1.0 ||
        !std::isfinite(rv.sigma) || rv.sigma <= 0.0) {
        throw std::domain_error("Invalid gamma or solidity");
    }
    // All vars in rv are in the local normalization
    // Nondimensionalization was a bad idea but I am rolling with it
    auto cth2 = rv.cth2;
    auto u1 = rv.u1;
    auto cm1 = rv.cm1;
    auto k = rv.k;
    auto lambda = rv.lambda;
    auto sigma = rv.sigma;
    auto cth1 = rv.cth1;

    double u2 = lambda * u1;
    double cm2 = k * cm1;

    double tau = 1. + u2 * cth2 - u1 * cth1;


    double cm1sq = cm1 * cm1;
    double cth1sq = cth1 * cth1;
    double c1sq = cm1sq + cth1sq;
    double w1sq = cm1sq + (cth1 - u1)*(cth1 - u1);

    double cm2sq = cm2 * cm2;
    double cth2sq = cth2 * cth2;
    double c2sq = cm2sq + cth2sq;
    double w2sq = cm2sq + (cth2 - u2)*(cth2 - u2);

    // Both static-temperature ratios use the row inlet TOTAL temperature
    // as their reference, and subtract ABSOLUTE flow kinetic energy.
    const double theta1 = 1.0 - c1sq / 2.0;
    const double theta2 = tau - c2sq / 2.0;
    if (!std::isfinite(tau) || tau <= 0.0 ||
        !std::isfinite(theta1) || theta1 <= 0.0 ||
        !std::isfinite(theta2) || theta2 <= 0.0 ||
        !std::isfinite(w1sq) || w1sq <= 0.0 ||
        !std::isfinite(w2sq)) {
        throw std::domain_error("Invalid stage thermodynamic state");
    }

    double Mw1sq = w1sq / ((gamma - 1.) * theta1);
    double Mw2sq = w2sq / ((gamma - 1.) * theta2);

    double w1 = std::sqrt(w1sq);
    double wth1 = cth1 - u1;
    double w2 = std::sqrt(w2sq);
    double wth2 = cth2 - u2;

    double Mw1 = std::sqrt(Mw1sq);
    double Mw2 = std::sqrt(Mw2sq);

    double M2sq = c2sq / ((gamma - 1.) * theta2);
    double M2 = std::sqrt(M2sq);

    double DF = 1. - w2/w1 + std::abs(wth1 - wth2) / (2*sigma*w1);

    StagePerformance performance{
        tau,
        Mw1,
        Mw2,
        M2,
        DF
    };

    return performance;
}
