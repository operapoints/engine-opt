#include "rostack.hpp"

// Physics models

StagePerformance rostack_physics_simple(const RotorVariables& rv, double gamma){
    // All vars in rv are in the local normalization
    // Nondimensionalization was a bad idea but I am rolling with it
    auto cth2 = rv.cth2;
    auto u1 = rv.u1;
    auto cm1 = rv.cm1;
    auto k = rv.k;
    auto lambda = rv.lambda;
    auto sigma = rv.sigma;
    auto nu = rv.nu;
    auto cth1 = rv.cth1;

    double u2 = lambda * u1;
    double cm2 = k * cm1;

    double tau = 1. + u2 * cth2 - u1 - cth1;

    double cm1sq = cm1 * cm1;
    double w1sq = cm1sq + (cth1 - u1)*(cth1 - u1);

    double cm2sq = cm2 * cm2;
    double w2sq = cm2sq + (cth2 - u2)*(cth2 - u2);

    double Mw1sq = w1sq / ((gamma - 1.) * (1 - w1sq/2.));
    double Mw2sq = w2sq / ((gamma - 1.) * (tau - w2sq/2.));

    double w1 = std::sqrt(w1sq);
    double wth1 = cth1 - u1;
    double w2 = std::sqrt(w2sq);
    double wth2 = cth2 - u1;

    double DF = 1. - w1/w2 + std::abs(wth1 - wth2) / (2*sigma*w1);
}