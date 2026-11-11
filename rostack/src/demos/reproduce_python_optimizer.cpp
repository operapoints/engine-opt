// Reproduce the attached python_optimizer.py model using pagmo.
// Standalone build (C++20, pagmo built with NLopt):
// g++ -std=c++20 -O3 -pthread reproduce_python_optimizer.cpp -lpagmo -o reproduce_python_optimizer
// ./reproduce_python_optimizer                         # M=1.5, work, independent speeds
// ./reproduce_python_optimizer 1.2                     # Reference Mach limit
// ./reproduce_python_optimizer 1.5 shock_pressure       # Reoptimize shock-proxy pressure
// ./reproduce_python_optimizer 1.5 work equal           # Equal speed magnitudes
// ./reproduce_python_optimizer --original              # Original Python limits only
// ./reproduce_python_optimizer --compare               # Run original AND added Mach limits
// Existing CMake project (after find_package(Pagmo REQUIRED)):
// add_executable(reproduce_python_optimizer src/demos/reproduce_python_optimizer.cpp)
// target_compile_features(reproduce_python_optimizer PRIVATE cxx_std_20)
// target_link_libraries(reproduce_python_optimizer PRIVATE Pagmo::pagmo)
// Do not also compile rostack_optimizer.cpp into this target: this file already
// contains the optimizer routine and does not require rostack.hpp.
//
// Original-limits Python reference (--original), M=1.5, independent speeds, work:
// co:      146.07340404 kJ/kg, tau=1.48602678, shock-proxy PR=3.14382641
// counter: 169.98225408 kJ/kg, tau=1.56557816, shock-proxy PR=3.77334900
// C++ equations match these fixed-design evaluations. The optimizer differs
// from Python's multistart SLSQP, so independently compare converged results.
// Model-only build for equation comparisons, with no pagmo dependencies:
// g++ -std=c++20 -O2 -DPYTHON_MODEL_ONLY reproduce_python_optimizer.cpp -o python_model
// ./python_model --evaluate counter 1.5 487.71549083456426 352.9847784484794 202.19126874637783
// --evaluate arguments after Mach: U1 U2 S1 [S2], all in m/s, independent speeds.
//
// This is intentionally a separate UDP, not the earlier Rostack comparison.
// Co = R-S-R-S; counter = R-R. Both finish axially. Dmax=.50; inlet Mrel=1.5.
// Default additional limits: rotor outlet Mrel<=1.2, every absolute station M<=0.8.
// Stator relative Mach equals absolute Mach. Downstream absolute inlet limits
// are enforced by the preceding outlet limit. --original disables only these
// additional limits; the Python geometry/work/temperature constraints remain.
// Speeds are independent by default. T1=288 K is STATIC and Cx=150 m/s.
// "Tip Mach" describes the representative streamline; no spanwise model exists.
// Shock multiplier 1.15 is the Python's uncalibrated assumption, not new evidence.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct PythonStudy {
    bool co = false;
    bool equal_speeds = false;
    bool shock_pressure_objective = false;
    bool additional_mach_caps = true;
    double mach_limit = 1.5;
    double outlet_relative_mach_limit = 1.2;
    double absolute_mach_limit = 0.8;
    double Cx = 150.0;
    double T1 = 288.0; // Inlet STATIC temperature, not total.
    double gamma = 1.4;
    double R = 287.0;
    double solidity = 1.5;
    double Dmax = 0.50;
    double shock_multiplier = 1.15;

    double cp() const { return gamma * R / (gamma - 1.0); }
    double Tt1() const { return T1 + Cx * Cx / (2.0 * cp()); }
    double U1_bound() const { return std::sqrt(mach_limit * mach_limit * gamma * R * T1 - Cx * Cx); }
    double U2_bound() const {
        const double hottest_T2 = T1 + U1_bound() * U1_bound() / cp();
        return std::sqrt(mach_limit * mach_limit * gamma * R * hottest_T2 - Cx * Cx);
    }
    double scale() const { return std::max(U1_bound(), U2_bound()); }
    std::size_t dimension() const { return (equal_speeds ? 2u : 3u) + (co ? 1u : 0u); }
    std::size_t base_constraint_count() const { return co ? 10u : 7u; }
    std::size_t get_nic() const {
        return base_constraint_count() + (additional_mach_caps ? (co ? 7u : 5u) : 0u);
    }
    std::string get_name() const {
        return co ? "Python reproduction: R1-S1-R2-S2" : "Python reproduction: counterrotating R1-R2";
    }
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const {
        const double a = U1_bound() / scale(), b = U2_bound() / scale();
        if (equal_speeds) {
            if (co) return {{1e-8, 0.0, 0.0}, {a, a, a}};
            return {{1e-8, 0.0}, {a, a}};
        }
        if (co) return {{1e-8, 1e-8, 0.0, 0.0}, {a, b, a, b}};
        return {{1e-8, 1e-8, 0.0}, {a, b, a}};
    }
    static double diffusion(double vin, double vout, double turning, double sigma) {
        return 1.0 - vout / vin + std::abs(turning) / (2.0 * sigma * vin);
    }
    double shock_recovery(double mach) const {
        const double m2 = std::max(mach, 1.0) * std::max(mach, 1.0);
        return std::pow((gamma + 1.0) * m2 / ((gamma - 1.0) * m2 + 2.0), gamma / (gamma - 1.0))
             * std::pow((gamma + 1.0) / (2.0 * gamma * m2 - (gamma - 1.0)), 1.0 / (gamma - 1.0));
    }
    double rotor_recovery(double mach) const {
        return std::exp(0.5 * (std::log(shock_recovery(mach))
                             + std::log(shock_recovery(shock_multiplier * mach))));
    }
    struct Evaluation {
        double U1, U2, S1, S2;
        double work1, work2, T_R1out, T_R2in, T_R2out, Ttout;
        double M1, M2, D1, D2;
        double Mrel1out, Mrel2out;
        double Mabs_inlet, Mabs1out, Mabs2in, Mabs2out;
        double Mabs_stator1out = std::numeric_limits<double>::quiet_NaN();
        double Mabs_stator2out = std::numeric_limits<double>::quiet_NaN();
        double DS1 = std::numeric_limits<double>::quiet_NaN();
        double DS2 = std::numeric_limits<double>::quiet_NaN();
        double MS1 = std::numeric_limits<double>::quiet_NaN();
        double MS2 = std::numeric_limits<double>::quiet_NaN();
        double rotor_recovery1, rotor_recovery2;
        double stator_recovery1 = 1.0, stator_recovery2 = 1.0;
        double log_recovery, log_pr_zero, log_pr_shock;
        std::array<double, 17> margins{}; // Python convention: feasible >= 0.
    };
    Evaluation evaluate(const std::vector<double>& x) const {
        if (x.size() != dimension()) throw std::invalid_argument("Wrong decision-vector dimension");
        const double tiny = 1e-8;
        Evaluation e{};
        const double velocity_scale = scale();
        e.U1 = x[0] * velocity_scale;
        if (equal_speeds) {
            e.U2 = e.U1;
            e.S1 = x[1] * velocity_scale;
            e.S2 = co ? x[2] * velocity_scale : e.S1;
        } else {
            e.U2 = x[1] * velocity_scale;
            e.S1 = x[2] * velocity_scale;
            e.S2 = co ? x[3] * velocity_scale : e.S1;
        }
        // U2 is a positive magnitude; counterrotor signed speed is -U2.
        // Co: each rotor adds S_i from zero inlet swirl; its stator removes S_i.
        // Counter: R2 receives S1 and removes it to exactly zero exit swirl.
        e.work1 = e.U1 * e.S1;
        e.work2 = e.U2 * e.S2;
        e.T_R1out = T1 + (e.work1 - e.S1 * e.S1 / 2.0) / cp();
        e.T_R2in = co ? T1 + e.work1 / cp() : e.T_R1out;
        e.Ttout = Tt1() + (e.work1 + e.work2) / cp();
        e.T_R2out = e.Ttout - (Cx * Cx + (co ? e.S2 * e.S2 : 0.0)) / (2.0 * cp());
        const double W1in = std::hypot(Cx, e.U1);
        const double W1out = std::hypot(Cx, e.S1 - e.U1);
        const double W2in = std::hypot(Cx, co ? e.U2 : e.S1 + e.U2);
        const double W2out = std::hypot(Cx, co ? e.S2 - e.U2 : e.U2);
        e.M1 = W1in / std::sqrt(gamma * R * T1);
        e.M2 = W2in / std::sqrt(gamma * R * std::max(e.T_R2in, tiny));
        e.Mrel1out = W1out / std::sqrt(gamma * R * std::max(e.T_R1out, tiny));
        e.Mrel2out = W2out / std::sqrt(gamma * R * std::max(e.T_R2out, tiny));
        e.Mabs_inlet = Cx / std::sqrt(gamma * R * T1);
        e.Mabs1out = std::hypot(Cx, e.S1) / std::sqrt(gamma * R * std::max(e.T_R1out, tiny));
        e.Mabs2in = std::hypot(Cx, co ? 0.0 : e.S1) / std::sqrt(gamma * R * std::max(e.T_R2in, tiny));
        e.Mabs2out = std::hypot(Cx, co ? e.S2 : 0.0) / std::sqrt(gamma * R * std::max(e.T_R2out, tiny));
        e.D1 = diffusion(W1in, W1out, e.S1, solidity);
        e.D2 = diffusion(W2in, W2out, e.S2, solidity);
        e.margins = {mach_limit - e.M1, mach_limit - e.M2,
                     Dmax - e.D1, Dmax - e.D2, (e.U1 - e.S1) / velocity_scale,
                     (e.T_R1out - tiny) / T1, (e.T_R2out - tiny) / T1, 0.0, 0.0, 0.0};
        if (co) {
            e.DS1 = diffusion(std::hypot(Cx, e.S1), Cx, e.S1, solidity);
            e.DS2 = diffusion(std::hypot(Cx, e.S2), Cx, e.S2, solidity);
            e.margins[7] = Dmax - e.DS1;
            e.margins[8] = Dmax - e.DS2;
            e.margins[9] = (e.U2 - e.S2) / velocity_scale;
            e.MS1 = std::hypot(Cx, e.S1) / std::sqrt(gamma * R * std::max(e.T_R1out, tiny));
            e.MS2 = std::hypot(Cx, e.S2) / std::sqrt(gamma * R * std::max(e.T_R2out, tiny));
            e.stator_recovery1 = shock_recovery(e.MS1);
            e.stator_recovery2 = shock_recovery(e.MS2);
            e.Mabs_stator1out = e.Mabs2in;
            const double final_static_temperature = e.Ttout - Cx * Cx / (2.0 * cp());
            e.Mabs_stator2out = Cx / std::sqrt(gamma * R * std::max(final_static_temperature, tiny));
        }
        if (additional_mach_caps) {
            std::size_t j = base_constraint_count();
            e.margins[j++] = outlet_relative_mach_limit - e.Mrel1out;
            e.margins[j++] = outlet_relative_mach_limit - e.Mrel2out;
            e.margins[j++] = absolute_mach_limit - e.Mabs_inlet;
            e.margins[j++] = absolute_mach_limit - e.Mabs1out;
            e.margins[j++] = absolute_mach_limit - e.Mabs2out;
            if (co) {
                e.margins[j++] = absolute_mach_limit - e.Mabs_stator1out;
                e.margins[j++] = absolute_mach_limit - e.Mabs_stator2out;
            }
        }
        e.rotor_recovery1 = rotor_recovery(e.M1);
        e.rotor_recovery2 = rotor_recovery(e.M2);
        e.log_recovery = std::log(e.rotor_recovery1) + std::log(e.rotor_recovery2)
                       + std::log(e.stator_recovery1) + std::log(e.stator_recovery2);
        e.log_pr_zero = gamma / (gamma - 1.0) * std::log(e.Ttout / Tt1());
        e.log_pr_shock = e.log_pr_zero + e.log_recovery;
        return e;
    }
    std::vector<double> fitness(const std::vector<double>& x) const {
        const auto e = evaluate(x);
        std::vector<double> f;
        f.reserve(1 + get_nic());
        f.push_back(shock_pressure_objective ? -e.log_pr_shock : -(e.work1 + e.work2) / 100000.0);
        for (std::size_t i = 0; i < get_nic(); ++i) f.push_back(-e.margins[i]);
        for (const double value : f) {
            if (!std::isfinite(value)) return std::vector<double>(1 + get_nic(), 1e6);
        }
        return f;
    }
    std::vector<double> coordinates(double U1, double U2, double S1, double S2) const {
        const double s = scale();
        if (equal_speeds) {
            if (co) return {U1 / s, S1 / s, S2 / s};
            return {U1 / s, S1 / s};
        }
        if (co) return {U1 / s, U2 / s, S1 / s, S2 / s};
        return {U1 / s, U2 / s, S1 / s};
    }
    std::vector<std::string> constraint_names() const {
        std::vector<std::string> names{"rotor1_Mrel_in", "rotor2_Mrel_in", "rotor1_diffusion", "rotor2_diffusion",
            "rotor1_swirl_le_speed", "rotor1_static_temperature", "rotor2_static_temperature"};
        if (co) names.insert(names.end(), {"stator1_diffusion", "stator2_diffusion", "rotor2_swirl_le_speed"});
        if (additional_mach_caps) {
            names.insert(names.end(), {"rotor1_Mrel_out", "rotor2_Mrel_out", "module_Mabs_in",
                                     "rotor1_Mabs_out", "rotor2_Mabs_out"});
            if (co) names.insert(names.end(), {"stator1_Mabs_out", "stator2_Mabs_out"});
        }
        return names;
    }
    std::map<std::string, double> report(const std::vector<double>& x) const {
        const auto e = evaluate(x);
        std::map<std::string, double> r{
            {"mach_limit", mach_limit}, {"rotor1_speed_m_s", e.U1}, {"rotor2_speed_m_s", e.U2},
            {"additional_mach_caps_enabled", additional_mach_caps ? 1.0 : 0.0},
            {"rotor_outlet_relative_mach_limit", outlet_relative_mach_limit},
            {"absolute_mach_limit", absolute_mach_limit},
            {"rotor2_signed_speed_m_s", co ? e.U2 : -e.U2},
            {"rotor1_exit_swirl_m_s", e.S1}, {"rotor2_exit_swirl_m_s", co ? e.S2 : 0.0},
            {"delivered_exit_swirl_m_s", 0.0},
            {"rotor1_inlet_relative_mach", e.M1}, {"rotor2_inlet_relative_mach", e.M2},
            {"rotor1_outlet_relative_mach", e.Mrel1out}, {"rotor2_outlet_relative_mach", e.Mrel2out},
            {"module_inlet_absolute_mach", e.Mabs_inlet},
            {"rotor1_outlet_absolute_mach", e.Mabs1out},
            {"rotor2_inlet_absolute_mach", e.Mabs2in},
            {"rotor2_outlet_absolute_mach", e.Mabs2out},
            {"stator1_outlet_absolute_mach", e.Mabs_stator1out},
            {"stator2_outlet_absolute_mach", e.Mabs_stator2out},
            {"shock_mach_multiplier", shock_multiplier},
            {"rotor1_impingement_mach", shock_multiplier * e.M1},
            {"rotor2_impingement_mach", shock_multiplier * e.M2},
            {"rotor1_shock_recovery", e.rotor_recovery1}, {"rotor2_shock_recovery", e.rotor_recovery2},
            {"rotor1_shock_entropy_J_kgK", -R * std::log(e.rotor_recovery1)},
            {"rotor2_shock_entropy_J_kgK", -R * std::log(e.rotor_recovery2)},
            {"rotor1_diffusion_factor", e.D1}, {"rotor2_diffusion_factor", e.D2},
            {"stator1_diffusion_factor", e.DS1}, {"stator2_diffusion_factor", e.DS2},
            {"stator1_inlet_mach", e.MS1}, {"stator2_inlet_mach", e.MS2},
            {"rotor1_work_kj_kg", e.work1 / 1000.0}, {"rotor2_work_kj_kg", e.work2 / 1000.0},
            {"total_work_kj_kg", (e.work1 + e.work2) / 1000.0},
            {"tau_overall", e.Ttout / Tt1()},
            {"pressure_ratio_zero_losses", std::exp(e.log_pr_zero)},
            {"pressure_ratio_shock_proxy", std::exp(e.log_pr_shock)},
            {"shock_recovery_product", std::exp(e.log_recovery)},
            {"outlet_total_temperature_K", e.Ttout},
            {"minimum_constraint_margin", *std::min_element(e.margins.begin(), e.margins.begin() + get_nic())}
        };
        const double signed_u2 = co ? e.U2 : -e.U2;
        const double ci2 = co ? 0.0 : e.S1;
        const double ce2 = co ? e.S2 : 0.0;
        r["rotor2_relative_inlet_angle_deg"] = std::atan2(ci2 - signed_u2, Cx) * 180.0 / 3.14159265358979323846;
        r["rotor2_relative_outlet_angle_deg"] = std::atan2(ce2 - signed_u2, Cx) * 180.0 / 3.14159265358979323846;
        r["rotor2_relative_turning_deg"] = std::abs(r["rotor2_relative_outlet_angle_deg"] - r["rotor2_relative_inlet_angle_deg"]);
        const auto names = constraint_names();
        double active_count = 0.0;
        for (std::size_t i = 0; i < get_nic(); ++i) {
            const double residual = -e.margins[i];
            const bool violated = !std::isfinite(residual) || residual > 1e-7;
            const bool active = !violated && std::abs(residual) <= 1e-6;
            r["slack: " + names[i]] = e.margins[i];
            r["residual: " + names[i]] = residual;
            r["active: " + names[i]] = active ? 1.0 : 0.0;
            r["violated: " + names[i]] = violated ? 1.0 : 0.0;
            active_count += active ? 1.0 : 0.0;
        }
        r["active_constraint_count"] = active_count;
        r["maximum_constraint_violation"] = std::max(0.0, -r.at("minimum_constraint_margin"));
        return r;
    }
};

void print_report(const std::map<std::string, double>& data)
{
    std::cout << std::setprecision(17) << "{\n";
    std::size_t n = 0;
    for (const auto& [name, value] : data) {
        std::cout << "  " << std::quoted(name) << ": ";
        if (std::isfinite(value)) std::cout << value;
        else std::cout << "null";
        std::cout << (++n == data.size() ? "\n" : ",\n");
    }
    std::cout << "}\n";
}

void print_constraint_inspection(const PythonStudy& study, const std::vector<double>& x)
{
    const auto e = study.evaluate(x);
    const auto names = study.constraint_names();
    const auto flags = std::cout.flags();
    const auto precision = std::cout.precision();
    std::cout << "Constraint inspection: residual <= 0; feasibility tolerance 1e-7; active threshold 1e-6\n";
    std::cout << std::left << std::setw(32) << "Constraint" << std::right << std::setw(18) << "Residual" << "  Status\n";
    for (std::size_t i = 0; i < study.get_nic(); ++i) {
        const double residual = -e.margins[i];
        const char* status = !std::isfinite(residual) || residual > 1e-7 ? "VIOLATED"
            : std::abs(residual) <= 1e-6 ? "ACTIVE" : "slack";
        std::cout << std::left << std::setw(32) << names[i] << std::right
                  << std::scientific << std::setprecision(8) << std::setw(18) << residual << "  " << status << '\n';
    }
    std::cout.flags(flags);
    std::cout.precision(precision);
}

#ifndef PYTHON_MODEL_ONLY

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/algorithms/nlopt.hpp>
#include <pagmo/archipelago.hpp>
#include <pagmo/population.hpp>
#include <pagmo/topologies/fully_connected.hpp>
#include <pagmo/topology.hpp>
#include <pagmo/utils/constrained.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


struct OptimizationResult {
    pagmo::vector_double x;
    pagmo::vector_double fitness;
    bool feasible = false;
    double elapsed_seconds = 0.0;
};
namespace {
struct ConstraintDiagnostics {
    std::size_t violated = 0;
    double max_violation = 0.0; // Excess beyond the configured tolerance.
    double max_inequality = -std::numeric_limits<double>::infinity();
};

ConstraintDiagnostics constraint_diagnostics(const pagmo::vector_double& f,
                                             const pagmo::problem& p) {
    ConstraintDiagnostics d;
    const auto& tol = p.get_c_tol();
    for (std::size_t j = 0; j < tol.size(); ++j) {
        const double c = f[p.get_nobj() + j];
        if (j >= p.get_nec()) d.max_inequality = std::max(d.max_inequality, c);
        const double residual = j < p.get_nec() ? std::abs(c) : c;
        const double violation = std::isfinite(c)
            ? std::max(0.0, residual - tol[j])
            : std::numeric_limits<double>::infinity();
        if (violation > 0.0) ++d.violated;
        d.max_violation = std::max(d.max_violation, violation);
    }
    return d;
}

void print_fitness_diagnostics(const std::string& label,
                               const pagmo::vector_double& f,
                               const pagmo::problem& p) {
    const auto d = constraint_diagnostics(f, p);
    std::cout << label << ": objective=" << f[0]
              << ", feasible=" << (p.feasibility_f(f) ? "yes" : "no")
              << ", violated=" << d.violated
              << ", max violation beyond tolerance=" << d.max_violation
              << ", largest inequality residual=" << d.max_inequality << '\n';
}

// Agreement is measured in objective space among feasible island champions only.
// It is a convergence proxy, not a certificate of global optimality; migration
// means these are not independent statistical samples.
struct ObjectiveAgreement {
    std::size_t count = 0;
    std::size_t near_best = 0;
    std::size_t top_count = 0;
    double median = 0.0;
    double iqr = 0.0;
    double top_spread_percent = 0.0;
};

ObjectiveAgreement objective_agreement(std::vector<double> values) {
    ObjectiveAgreement result;
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double v) { return !std::isfinite(v); }), values.end());
    if (values.empty()) return result;
    std::sort(values.begin(), values.end());
    result.count = values.size();
    const auto quantile = [&](double fraction) {
        const double position = fraction * static_cast<double>(values.size() - 1);
        const auto lo = static_cast<std::size_t>(position);
        const auto hi = std::min(lo + 1, values.size() - 1);
        return values[lo] + (position - lo) * (values[hi] - values[lo]);
    };
    result.median = quantile(0.5);
    result.iqr = quantile(0.75) - quantile(0.25);
    // Relative gaps also work for negative objectives. Near zero, use a small
    // absolute scale floor to avoid division by zero.
    const double scale = std::max(std::abs(values.front()), 1e-12);
    const double near_best_gap = 0.001 * scale; // Within 0.1% of the best.
    for (double value : values) {
        if (value - values.front() <= near_best_gap) ++result.near_best;
    }
    if (values.size() >= 4) {
        result.top_count = (values.size() + 3) / 4; // Best quarter, rounded up.
        result.top_spread_percent = 100.0
            * (values[result.top_count - 1] - values.front()) / scale;
    }
    return result;
}

void print_objective_agreement(const std::vector<pagmo::vector_double>& fitness,
                               const pagmo::problem& p, std::size_t island_count) {
    std::vector<double> feasible_objectives;
    for (std::size_t i = 0; i < std::min(island_count, fitness.size()); ++i) {
        if (p.feasibility_f(fitness[i])) feasible_objectives.push_back(fitness[i][0]);
    }
    const auto a = objective_agreement(std::move(feasible_objectives));
    std::cout << "  Objective agreement (feasible islands): ";
    if (a.count == 0) {
        std::cout << "n/a (none feasible)\n";
        return;
    }
    std::cout << "median=" << a.median << ", IQR=" << a.iqr;
    if (a.top_count > 0) {
        std::cout << ", best-quarter spread=" << a.top_spread_percent
                  << "% (n=" << a.top_count << ')';
    } else {
        std::cout << ", best-quarter spread=n/a (fewer than 4 feasible islands)";
    }
    std::cout << ", within 0.1% of best=" << a.near_best << '/' << a.count << '\n';
}

void print_progress_summary(const std::string& label,
                            const std::vector<pagmo::vector_double>& fitness,
                            const pagmo::problem& p, std::size_t island_count,
                            const std::string& extra = "") {
    if (fitness.empty()) return;
    const auto order = pagmo::sort_population_con(fitness, p.get_nec(), p.get_c_tol());
    const auto best = order[0];
    std::size_t feasible_count = 0;
    double best_feasible = std::numeric_limits<double>::infinity();
    double largest_objective = -std::numeric_limits<double>::infinity();
    for (const auto& f : fitness) {
        largest_objective = std::max(largest_objective, f[0]);
        if (p.feasibility_f(f)) {
            ++feasible_count;
            best_feasible = std::min(best_feasible, f[0]);
        }
    }
    std::cout << label << " — feasible champions: " << feasible_count << '/' << fitness.size()
              << ", best feasible objective: ";
    if (feasible_count > 0) std::cout << best_feasible;
    else std::cout << "n/a";
    std::cout << ", largest champion objective: " << largest_objective << extra << '\n';
    const std::string best_label = "  Best champion by constrained ranking (island "
        + std::to_string(best + 1) + ")";
    print_fitness_diagnostics(best_label, fitness[best], p);
    print_objective_agreement(fitness, p, island_count);
}
} // namespace

OptimizationResult optimize_udp(pagmo::problem problem)
{
    const auto start = std::chrono::steady_clock::now();
    problem.set_c_tol(1e-7);
    std::cout << std::setprecision(17) << problem;

    pagmo::algorithm global_algorithm{pagmo::gaco(
        100,     // generations per evolution
        31,      // archive size
        1.0,     // q
        100.0,   // oracle
        0.01,    // accuracy
        100,     // convergence-speed threshold
        7,       // generation mark
        100,     // improvement stop
        100000,  // evaluation stop
        0.0,     // focus
        true     // memory
    )};
    pagmo::topology topology{pagmo::fully_connected{}};
    pagmo::archipelago archipelago{topology, 32u, global_algorithm, problem, 6000u};

    // Populations are initialized by pagmo; no supplied incumbent or seed.
    for (int round = 0; round < 8; ++round) {
        archipelago.evolve(1);
        archipelago.wait_check();
        const auto fitness = archipelago.get_champions_f();
        print_progress_summary("Round " + std::to_string(round + 1) + " / 8",
                               fitness, problem, fitness.size());
    }

    pagmo::algorithm local_algorithm{pagmo::nlopt{"cobyla"}};
    auto* nlopt = local_algorithm.extract<pagmo::nlopt>();
    nlopt->set_maxeval(2000);
    nlopt->set_xtol_rel(1e-10);
    nlopt->set_ftol_rel(1e-10);
    nlopt->set_xtol_abs(1e-12);
    local_algorithm.set_verbosity(0);

    auto refined_x = archipelago.get_champions_x();
    std::vector<pagmo::vector_double> refined_f;
    refined_f.reserve(refined_x.size());
    std::vector<pagmo::population> local_populations;
    local_populations.reserve(refined_x.size());
    for (const auto& x : refined_x) {
        local_populations.emplace_back(problem, 0u);
        local_populations.back().push_back(x);
        refined_f.push_back(local_populations.back().champion_f());
    }

    std::vector<bool> active(refined_x.size(), true);
    for (int pass = 0; pass < 20; ++pass) {
        std::size_t attempted = 0;
        std::size_t errors = 0;
        std::string first_error;
        for (std::size_t i = 0; i < local_populations.size(); ++i) {
            if (!active[i]) continue;
            ++attempted;
            auto& population = local_populations[i];
            const double previous_objective = population.champion_f()[0];
            try {
                population = local_algorithm.evolve(population);
            } catch (const std::exception& e) {
                ++errors;
                active[i] = false;
                if (first_error.empty()) {
                    first_error = "island " + std::to_string(i + 1) + ": " + e.what();
                    std::replace(first_error.begin(), first_error.end(), '\n', ' ');
                    std::replace(first_error.begin(), first_error.end(), '\r', ' ');
                }
                continue; // Preserve this island's previous result.
            }
            refined_x[i] = population.champion_x();
            refined_f[i] = population.champion_f();
            if (std::abs(previous_objective - refined_f[i][0])
                < 1e-8 * std::abs(previous_objective)) {
                active[i] = false;
            }
        }
        if (attempted == 0) break;
        const auto remaining = std::count(active.begin(), active.end(), true);
        std::string extra = ", polished=" + std::to_string(attempted)
            + ", remaining active=" + std::to_string(remaining)
            + ", errors=" + std::to_string(errors);
        if (errors > 0) extra += ", first error=" + first_error;
        print_progress_summary("Local polish pass " + std::to_string(pass + 1) + " / 20",
                               refined_f, problem, refined_x.size(), extra);
        if (remaining == 0) break;
    }

    const auto ranking = pagmo::sort_population_con(
        refined_f, problem.get_nec(), problem.get_c_tol());
    const auto best = ranking.front();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    OptimizationResult result{refined_x[best], refined_f[best],
                              problem.feasibility_f(refined_f[best]), elapsed};

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
              << "No. refined champions: " << refined_f.size() << '\n'
              << "Rank 1: island " << best + 1 << '\n';
    print_fitness_diagnostics("Final candidate", result.fitness, problem);
    if (!result.feasible) {
        std::cout << "No feasible optimum found; returning the best constrained-ranked candidate.\n";
    }
    std::ostringstream duration;
    duration << std::fixed << std::setprecision(3) << elapsed;
    std::cout << "Optimization time: " << duration.str() << " s\n";
    return result;
}

#endif

int main(int argc, char** argv) try
{
    bool original = false, compare = false;
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--original") original = true;
        else if (arg == "--compare") compare = true;
        else args.push_back(arg);
    }
    if (original && compare) throw std::invalid_argument("Choose --original OR --compare");
    if (!args.empty() && args[0] == "--evaluate") {
        if (args.size() < 6 || args.size() > 7 || compare) {
            throw std::invalid_argument("Usage: --evaluate co|counter Mach U1 U2 S1 [S2] [--original]");
        }
        PythonStudy study;
        const std::string architecture = args[1];
        if (architecture != "co" && architecture != "counter") throw std::invalid_argument("Unknown architecture");
        study.co = architecture == "co";
        study.additional_mach_caps = !original;
        if (study.co && args.size() != 7) throw std::invalid_argument("Co evaluation requires S2");
        study.mach_limit = std::stod(args[2]);
        const double U1 = std::stod(args[3]), U2 = std::stod(args[4]), S1 = std::stod(args[5]);
        const double S2 = args.size() == 7 ? std::stod(args[6]) : 0.0;
        if (!std::isfinite(study.mach_limit) || study.mach_limit <= study.Cx/std::sqrt(study.gamma*study.R*study.T1))
            throw std::invalid_argument("Mach limit must exceed inlet axial Mach");
        print_report(study.report(study.coordinates(U1,U2,S1,S2)));
        return 0;
    }
#ifndef PYTHON_MODEL_ONLY
    const double limit = !args.empty() ? std::stod(args[0]) : 1.5;
    const std::string objective = args.size() > 1 ? args[1] : "work";
    const std::string speeds = args.size() > 2 ? args[2] : "independent";
    if (args.size() > 3 || (objective != "work" && objective != "shock_pressure") ||
        (speeds != "independent" && speeds != "equal")) {
        throw std::invalid_argument("Usage: program [Mach [work|shock_pressure [independent|equal]]] [--original|--compare]");
    }
    PythonStudy settings;
    if (!std::isfinite(limit) || limit <= settings.Cx/std::sqrt(settings.gamma*settings.R*settings.T1))
        throw std::invalid_argument("Mach limit must exceed inlet axial Mach");
    const std::vector<bool> modes = compare ? std::vector<bool>{false, true} : std::vector<bool>{!original};
    std::array<double,2> advantage{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    for (bool added_caps : modes) {
        std::cout << "\nLIMIT SET: " << (added_caps ? "Python plus prior-study Mach caps" : "Original Python") << '\n';
        std::array<std::map<std::string,double>,2> reports;
        std::array<bool,2> feasible{};
        for (std::size_t i = 0; i < 2; ++i) {
            PythonStudy udp;
            udp.co = i == 0;
            udp.additional_mach_caps = added_caps;
            udp.mach_limit = limit;
            udp.equal_speeds = speeds == "equal";
            udp.shock_pressure_objective = objective == "shock_pressure";
            std::cout << "\n============================================================\n"
                      << udp.get_name() << "\nAllowed inlet Mrel = " << limit
                      << ", diffusion <= " << udp.Dmax << ", objective = " << objective << ", speeds = " << speeds << '\n';
            if (added_caps) {
                std::cout << "Rotor outlet Mrel <= " << udp.outlet_relative_mach_limit
                          << "; absolute Mach at every station <= " << udp.absolute_mach_limit << '\n';
            }
            const auto best = optimize_udp(pagmo::problem{udp});
            feasible[i] = best.feasible;
            reports[i] = udp.report(best.x);
            print_report(reports[i]);
            print_constraint_inspection(udp, best.x);
        }
        std::cout << "\nComparison (co / counter):\n";
        for (const auto* name : {"total_work_kj_kg", "tau_overall", "pressure_ratio_zero_losses", "pressure_ratio_shock_proxy"}) {
            std::cout << name << ": " << reports[0].at(name) << " / " << reports[1].at(name) << '\n';
        }
        if (feasible[0] && feasible[1]) {
            const double percent = 100.0 * (reports[1].at("total_work_kj_kg") / reports[0].at("total_work_kj_kg") - 1.0);
            advantage[added_caps ? 1 : 0] = percent;
            std::cout << "Counterrotating work advantage: " << percent << "%\n";
        } else {
            std::cout << "At least one result is infeasible; comparison is not a feasible-design ranking.\n";
        }
    }
    if (compare && std::isfinite(advantage[0]) && std::isfinite(advantage[1])) {
        std::cout << "\nWork advantage, original -> added caps: " << advantage[0] << "% -> " << advantage[1]
                  << "% (change " << advantage[1] - advantage[0] << " percentage points)\n";
    }
    return 0;
#else
    std::cerr << "Model-only build: use --evaluate co|counter Mach U1 U2 S1 [S2] [--original].\n";
    return 1;
#endif
}
catch (const std::exception& e)
{
    std::cerr << "reproduce_python_optimizer: " << e.what() << '\n';
    return 1;
}
