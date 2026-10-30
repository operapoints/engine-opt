#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

// C++20. No pagmo dependency is needed: std::vector<double> is compatible
// with pagmo::vector_double. All configuration is copied into the UDP.
inline constexpr std::size_t NUM_ROSTACK_VARS = 8;
inline constexpr std::size_t NUM_PERFORMANCE_PARAMS = 5;
using Rostack_bounds =
    std::array<std::array<double, 2>, NUM_ROSTACK_VARS>;
using Rostack_live_vars = std::array<bool, NUM_ROSTACK_VARS>;

struct RotorVariables {
    double cth2, u1, cm1, k, lambda, sigma, nu, cth1;
};

struct StagePerformance {
    double tau; // Absolute outlet / inlet total temperature, for THIS row.
    double mach_rel_in;
    double mach_rel_out;
    double mach_abs_out;
    double diffusion;
};

struct FlowState {
    double tau_cumulative; // Current absolute total T / stack inlet total T.
    double cm;            // Velocity / sqrt(Cp * stack inlet total T).
    double ctheta;        // Velocity / sqrt(Cp * stack inlet total T).
};

struct StageEvaluation {
    StagePerformance performance;
    FlowState outlet;
    double radius_ratio; // R2/R1, including for stators.
    double u1_global;    // Blade speeds / sqrt(Cp * stack inlet total T).
    double u2_global;
};

struct StackPerformance {
    std::vector<StagePerformance> stages;
    double tau_total = 1.0;
};

// Supply your existing row equations through this hook. ALL four velocity
// inputs are normalized by sqrt(Cp * THIS ROW'S inlet absolute total T).
// k, lambda, sigma and nu are unchanged dimensionless parameters.
// Output Mach numbers use local static sound speed; tau is a LOCAL ratio.
// In the physics equations, the row's outlet blade speed is u2 = lambda * u1.
// The function must not mutate shared state if evaluated concurrently.
using Rostack_physics = StagePerformance (*)(const RotorVariables&, double gamma);

class Rostack_rotor {
public:
    Rostack_rotor(Rostack_bounds bounds, Rostack_live_vars live_vars,
                 bool inherits_inlet = false, bool inherits_speed = false)
        : bounds_(std::move(bounds)), live_vars_(live_vars),
          inherits_inlet_(inherits_inlet),
          is_stator_(!live_vars[1] && bounds_[1][0] == 0.0),
          inherits_speed_(inherits_speed && !is_stator_)
    {
        if (inherits_speed_) live_vars_[1] = false;
        if (inherits_inlet_) {
            live_vars_[2] = false; // cm1 comes from preceding row.
            live_vars_[7] = false; // cth1 comes from preceding row.
        }
        for (std::size_t i = 0; i < NUM_ROSTACK_VARS; ++i) {
            // Bounds for inherited inputs are ignored, including their minima.
            if (inherits_inlet_ && (i == 2 || i == 7)) continue;
            if (inherits_speed_ && i == 1) continue;
            if (!std::isfinite(bounds_[i][0]) ||
                (live_vars_[i] && (!std::isfinite(bounds_[i][1]) ||
                                  bounds_[i][0] > bounds_[i][1]))) {
                throw std::invalid_argument("Invalid rotor variable bounds");
            }
            n_opt_vars_ += live_vars_[i];
        }
    }

    std::size_t get_n_opt_vars() const noexcept { return n_opt_vars_; }
    bool is_stator() const noexcept { return is_stator_; }

    void append_bounds(std::vector<double>& lower,
                       std::vector<double>& upper) const
    {
        for (std::size_t i = 0; i < NUM_ROSTACK_VARS; ++i) {
            if (live_vars_[i]) {
                lower.push_back(bounds_[i][0]);
                upper.push_back(bounds_[i][1]);
            }
        }
    }

    StageEvaluation evaluate(std::span<const double> x,
                             const std::optional<FlowState>& incoming,
                             Rostack_physics physics, double gamma,
                             std::optional<double> shaft_u1_global = std::nullopt) const
    {
        if (!physics) throw std::logic_error("No rotor physics function supplied");
        if (x.size() != n_opt_vars_ || inherits_inlet_ != incoming.has_value()) {
            throw std::invalid_argument("Rotor input layout mismatch");
        }

        std::array<double, NUM_ROSTACK_VARS> global{};
        std::size_t j = 0;
        for (std::size_t i = 0; i < NUM_ROSTACK_VARS; ++i) {
            global[i] = live_vars_[i] ? x[j++] : bounds_[i][0];
        }
        const double tau_in = incoming ? incoming->tau_cumulative : 1.0;
        if (incoming) {
            global[2] = incoming->cm;
            global[7] = incoming->ctheta;
        }
        if (inherits_speed_) {
            if (!shaft_u1_global) {
                throw std::logic_error("Missing inherited spool speed");
            }
            global[1] = *shaft_u1_global;
        }
        if (!std::isfinite(tau_in) || tau_in <= 0.0) {
            throw std::domain_error("Non-positive or non-finite inlet total T");
        }
        for (double value : global) {
            if (!std::isfinite(value)) {
                throw std::domain_error("Non-finite rotor input");
            }
        }
        if (global[4] <= 0.0) {
            throw std::domain_error("Radius ratio must be positive");
        }
        const double u2_global = global[1] * global[4];
        if (!std::isfinite(u2_global)) {
            throw std::domain_error("Non-finite outlet blade speed");
        }

        // Global-normalized V -> row-inlet-normalized V. Cp is constant.
        const double factor = 1.0 / std::sqrt(tau_in);
        const RotorVariables local{
            global[0] * factor, global[1] * factor, global[2] * factor,
            global[3], global[4], global[5], global[6], global[7] * factor
        };
        const auto p = physics(local, gamma);
        if (!std::isfinite(p.tau) || p.tau <= 0.0 ||
            !std::isfinite(p.mach_rel_in) || p.mach_rel_in < 0.0 ||
            !std::isfinite(p.mach_rel_out) || p.mach_rel_out < 0.0 ||
            !std::isfinite(p.mach_abs_out) || p.mach_abs_out < 0.0 ||
            !std::isfinite(p.diffusion)) {
            throw std::domain_error("Invalid rotor performance");
        }
        // Keep the handoff on the GLOBAL normalization basis.
        const FlowState outlet{tau_in * p.tau, global[3] * global[2], global[0]};
        if (!std::isfinite(outlet.tau_cumulative) ||
            outlet.tau_cumulative <= 0.0 || !std::isfinite(outlet.cm)) {
            throw std::domain_error("Invalid rotor outlet state");
        }
        return {p, outlet, global[4], global[1], u2_global};
    }

private:
    Rostack_bounds bounds_;
    Rostack_live_vars live_vars_;
    bool inherits_inlet_;
    bool is_stator_;
    bool inherits_speed_;
    std::size_t n_opt_vars_ = 0;
};

// A configuration object containing consecutive components on one shaft.
// The first NON-STATOR supplies the independent u1 (fixed or optimized).
// All later non-stators have u1 overridden; stators always retain u1 = 0.
// Their lambda still changes the radius at the next component's inlet.
// Validation and the effective variable layout are finalized by add_spool().
class Rostack_spool {
public:
    void add_rotor(Rostack_bounds bounds, Rostack_live_vars live_vars)
    {
        components_.push_back({std::move(bounds), live_vars});
    }

    std::size_t stage_count() const noexcept { return components_.size(); }

private:
    struct Component {
        Rostack_bounds bounds;
        Rostack_live_vars live_vars;
    };
    std::vector<Component> components_;
    friend class Rostack_model;
};

class Rostack_model {
public:
    Rostack_model() = default;
    explicit Rostack_model(Rostack_physics physics, double gamma = 1.4)
        : physics_(physics), gamma_(gamma)
    {
        if (!physics_ || !std::isfinite(gamma_) || gamma_ <= 1.0) {
            throw std::invalid_argument("Invalid physics function or gamma");
        }
    }

    // Copies the spool configuration and prepares a flat evaluation layout.
    // Downstream u1 bounds/flags are ignored within the same spool, except
    // that fixed zero identifies stators BEFORE any overrides are applied.
    // cm1/cth1 propagate across spool boundaries as well as within spools.
    void add_spool(const Rostack_spool& spool)
    {
        if (spool.components_.empty()) {
            throw std::invalid_argument("Cannot add an empty spool");
        }
        std::vector<Rostack_rotor> prepared;
        prepared.reserve(spool.stage_count());
        bool has_rotor = false;
        for (const auto& c : spool.components_) {
            const bool is_stator = !c.live_vars[1] && c.bounds[1][0] == 0.0;
            prepared.emplace_back(c.bounds, c.live_vars,
                                  !rotors_.empty() || !prepared.empty(),
                                  has_rotor && !is_stator);
            if (!is_stator) has_rotor = true;
        }

        // Allocate before changing the logical model. The prepared rotors
        // contain only value members, so the following moves cannot throw.
        rotors_.reserve(rotors_.size() + prepared.size());
        offsets_.reserve(offsets_.size() + prepared.size());
        spool_rows_.reserve(spool_rows_.size() + 1);
        for (auto& rotor : prepared) {
            offsets_.push_back(offsets_.back() + rotor.get_n_opt_vars());
            rotors_.push_back(std::move(rotor));
        }
        spool_rows_.push_back(rotors_.size());
    }

    // Compatibility convenience: each call creates a separate one-row spool,
    // preserving independently specified blade speeds from the old interface.
    // Use add_spool() to put several components on a shared shaft.
    void add_rotor(Rostack_bounds bounds, Rostack_live_vars live_vars)
    {
        Rostack_spool spool;
        spool.add_rotor(std::move(bounds), live_vars);
        add_spool(spool);
    }

    std::size_t stage_count() const noexcept { return rotors_.size(); }
    std::size_t spool_count() const noexcept { return spool_rows_.size() - 1; }
    std::size_t get_n_opt_vars() const noexcept { return offsets_.back(); }

    std::pair<std::vector<double>, std::vector<double>> get_bounds() const
    {
        std::vector<double> lower, upper;
        lower.reserve(get_n_opt_vars());
        upper.reserve(get_n_opt_vars());
        for (const auto& rotor : rotors_) rotor.append_bounds(lower, upper);
        return {std::move(lower), std::move(upper)};
    }

    StackPerformance evaluate(std::span<const double> x) const
    {
        if (!physics_ || rotors_.empty()) {
            throw std::logic_error("Configure physics and at least one row");
        }
        if (x.size() != get_n_opt_vars()) {
            throw std::invalid_argument("Wrong stack decision-vector size");
        }
        StackPerformance result;
        result.stages.reserve(stage_count());
        std::optional<FlowState> flow;
        for (std::size_t spool = 0; spool < spool_count(); ++spool) {
            // Shaft speeds use GLOBAL normalization, never the row-local basis.
            std::optional<double> shaft_u1_global;
            for (std::size_t i = spool_rows_[spool];
                 i < spool_rows_[spool + 1]; ++i) {
                const auto row = rotors_[i].evaluate(
                    x.subspan(offsets_[i], offsets_[i + 1] - offsets_[i]),
                    flow, physics_, gamma_, shaft_u1_global);
                result.stages.push_back(row.performance);
                flow = row.outlet;

                if (shaft_u1_global) {
                    // Includes radius changes THROUGH intervening stators.
                    *shaft_u1_global *= row.radius_ratio;
                } else if (!rotors_[i].is_stator()) {
                    shaft_u1_global = row.u2_global;
                }
                if (shaft_u1_global && !std::isfinite(*shaft_u1_global)) {
                    throw std::domain_error("Non-finite spool blade speed");
                }
            }
        }
        result.tau_total = flow->tau_cumulative;
        return result;
    }

private:
    std::vector<Rostack_rotor> rotors_;
    std::vector<std::size_t> offsets_{0};
    std::vector<std::size_t> spool_rows_{0}; // Flat row ranges per spool.
    Rostack_physics physics_ = nullptr;
    double gamma_ = 1.4;
};

using StageMetric = double StagePerformance::*;
enum class Relation { upper, lower, equal };

struct StageConstraint {
    std::size_t stage;
    StageMetric metric;
    Relation relation;
    double limit;
    double scale = 1.0;

    double residual(const StackPerformance& result) const
    {
        const double value = result.stages[stage].*metric;
        return (relation == Relation::lower ? limit - value : value - limit)
               / scale;
    }
};

inline StageConstraint upper(std::size_t row, StageMetric metric,
                             double limit, double scale = 1.0)
{ return {row, metric, Relation::upper, limit, scale}; }

inline StageConstraint lower(std::size_t row, StageMetric metric,
                             double limit, double scale = 1.0)
{ return {row, metric, Relation::lower, limit, scale}; }

inline StageConstraint equal(std::size_t row, StageMetric metric,
                             double target, double scale = 1.0)
{ return {row, metric, Relation::equal, target, scale}; }

// Illustrative objective ONLY. Replace with your chosen physical design goal.
struct MaximizeTemperatureRatio {
    double operator()(const StackPerformance& result) const
    { return -result.tau_total; }
};

// Objective must be default/copy/move constructible and callable as a const
// object with StackPerformance, returning one scalar to minimize.
// To change configuration, construct a new problem rather than mutating a UDP.
template<class Objective = MaximizeTemperatureRatio>
class Rostack_problem {
public:
    // Structural default constructor required by pagmo; configure before use.
    Rostack_problem() = default;

    Rostack_problem(Rostack_model model,
                    std::vector<StageConstraint> constraints,
                    Objective objective = {})
        : model_(std::move(model)), objective_(std::move(objective)),
          bounds_(model_.get_bounds())
    {
        if (model_.stage_count() == 0 || bounds_.first.empty()) {
            throw std::invalid_argument("Optimization needs rows and live variables");
        }
        for (const auto& c : constraints) {
            if (c.stage >= model_.stage_count() || !valid_metric(c.metric) ||
                !std::isfinite(c.limit) || !std::isfinite(c.scale) ||
                c.scale <= 0.0 ||
                (c.relation != Relation::equal && c.relation != Relation::upper &&
                 c.relation != Relation::lower)) {
                throw std::invalid_argument("Invalid stage constraint");
            }
            (c.relation == Relation::equal ? equalities_ : inequalities_)
                .push_back(c);
        }
    }

    std::size_t get_nobj() const noexcept { return 1; }
    std::size_t get_nec() const noexcept { return equalities_.size(); }
    std::size_t get_nic() const noexcept { return inequalities_.size(); }
    std::pair<std::vector<double>, std::vector<double>> get_bounds() const
    { return bounds_; }

    std::vector<double> fitness(const std::vector<double>& x) const
    {
        const auto result = model_.evaluate(std::span<const double>{x});
        std::vector<double> f;
        f.reserve(1 + get_nec() + get_nic());
        f.push_back(objective_(result));
        for (const auto& c : equalities_) f.push_back(c.residual(result));
        for (const auto& c : inequalities_) f.push_back(c.residual(result));
        for (double value : f) {
            if (!std::isfinite(value)) {
                throw std::domain_error("Non-finite objective or constraint");
            }
        }
        return f;
    }

private:
    static bool valid_metric(StageMetric m) noexcept
    {
        return m == &StagePerformance::tau ||
               m == &StagePerformance::mach_rel_in ||
               m == &StagePerformance::mach_rel_out ||
               m == &StagePerformance::mach_abs_out ||
               m == &StagePerformance::diffusion;
    }

    Rostack_model model_;
    Objective objective_;
    std::pair<std::vector<double>, std::vector<double>> bounds_;
    std::vector<StageConstraint> equalities_, inequalities_;
};


// Failure policy: invalid physical states throw std::domain_error. Optimizers
// will generally stop on an exception. Restrict the search to valid states,
// or add an explicit, fixed-length finite penalty policy for your equations.
// Do not blanket-catch logic/configuration errors as infeasible designs.
// This header supplies neither analytic derivatives nor Boost serialization.
