#include "LORd.hpp"

std::vector<double> LORd_model::mode_frequencies(LORd_model_params_t params){
    this->make_H0(params);// Construct the system Hessian and write to this->H0
    auto& H0 = this->H0;
}