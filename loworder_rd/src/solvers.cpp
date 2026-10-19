#include "LORd.hpp"
#include <Eigen/Dense>

Eigen::Matrix<std::complex<double>,16,1> LORd_model::mode_eigvals(LORd_model_params_t params){
    this->make_H0(params);// Construct the system Hessian and write to this->H0
    auto& H0 = this->H0;
    Eigen::Matrix<double,16,16> A = Eigen::Matrix<double,16,16>::Zero();
    const auto Am = H0.block<8,8>(0,0);
    const auto Bm = H0.block<8,8>(0,8);
    const auto Cm = H0.block<8,8>(8,0);
    const auto Dm = H0.block<8,8>(8,8);
    A.block<8,8>(0,8).setIdentity();
    Eigen::LLT<Eigen::Matrix<double,8,8>> M_factor(Dm);
    A.block<8,8>(8,0) = M_factor.solve(Am);
    A.block<8,8>(8,8) = M_factor.solve(Bm-Cm);
    Eigen::EigenSolver<Eigen::Matrix<double,16,16>> solver(A);
    Eigen::Matrix<std::complex<double>,16,1> lambda = solver.eigenvalues();
    return lambda;
}