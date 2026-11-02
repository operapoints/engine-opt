#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// C++20. No pagmo dependency is needed: std::vector<double> is compatible
// with pagmo::vector_double. All configuration is copied into the UDP.
inline constexpr std::size_t NUM_ROSTACK_VARS = 8;
inline constexpr std::size_t NUM_PERFORMANCE_PARAMS = 5;
using Rostack_bounds =
    std::array<std::array<double, 2>, NUM_ROSTACK_VARS>;
using Rostack_live_vars = std::array<bool, NUM_ROSTACK_VARS>;
using Rostack_values = std::map<std::string, double>; // Lexicographic ordering.

// An optional, per-call diagnostics sink. operator[] writes into the inspection
// map when enabled, and into a disposable scalar otherwise (no map allocation).
// Use it for diagnostic assignments only, never as working storage for physics.
// Do not retain this object or references obtained from it after the callback.
class Rostack_diagnostics {
public:
    Rostack_diagnostics() = default;
    explicit Rostack_diagnostics(Rostack_values& values) : values_(&values) {}
    explicit operator bool() const noexcept { return values_ != nullptr; }
    double& operator[](std::string_view key)
    {
        return values_ ? (*values_)[std::string(key)] : discarded_;
    }
private:
    Rostack_values* values_ = nullptr;
    double discarded_ = 0.0;
};

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
    bool valid = true;
    std::string error;
    std::optional<std::size_t> failed_stage; // Zero-based global row index.
};

struct Rostack_inspection {
    Rostack_values design;           // Resolved inputs, including fixed/inherited.
    Rostack_values constraints;      // Scaled residuals, exactly as in fitness().
    Rostack_values performance_data; // Standard outputs plus physics diagnostics.
    std::vector<double> fitness;     // Populated by the UDP, not the bare model.
    StackPerformance result;         // May contain only completed rows on failure.
    bool valid = true;               // Physical/objective validity, not feasibility.
    std::string error;
    std::vector<std::string> invalid_values; // Nonfinite entries in the maps.
};

namespace rostack_detail {
struct StageTrace {
    Rostack_values design;
    Rostack_values performance;
};

inline void insert_unique(Rostack_values& values, const std::string& key, double value)
{
    if (values.emplace(key, value).second) return;
    for (std::size_t n = 2; ; ++n) {
        if (values.emplace(key + " #" + std::to_string(n), value).second) return;
    }
}

inline void find_invalid_values(Rostack_inspection& report)
{
    report.invalid_values.clear();
    const auto scan = [&](std::string_view category, const Rostack_values& values) {
        for (const auto& [key, value] : values) {
            if (!std::isfinite(value)) {
                report.invalid_values.push_back(std::string(category) + ": " + key);
            }
        }
    };
    scan("Design", report.design);
    scan("Constraints", report.constraints);
    scan("Performance", report.performance_data);
}
} // namespace rostack_detail

// Supply your existing row equations through this hook. ALL four velocity
// inputs are normalized by sqrt(Cp * THIS ROW'S inlet absolute total T).
// k, lambda, sigma and nu are unchanged dimensionless parameters.
// Output Mach numbers use local static sound speed; tau is a LOCAL ratio.
// In the physics equations, the row's outlet blade speed is u2 = lambda * u1.
// The function must not mutate shared state if evaluated concurrently.
using Rostack_physics = StagePerformance (*)(const RotorVariables&, double gamma);
// Optional new callback form. Existing two-argument physics functions still work.
using Rostack_inspectable_physics = StagePerformance (*)(
    const RotorVariables&, double gamma, Rostack_diagnostics performance_params);

class Rostack_rotor {
public:
    Rostack_rotor(Rostack_bounds bounds, Rostack_live_vars live_vars,
                 std::string stage_name = {},
                 bool inherits_inlet = false, bool inherits_speed = false)
        : bounds_(std::move(bounds)), live_vars_(live_vars),
          name_(std::move(stage_name)),
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
    const std::string& name() const noexcept { return name_; }

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
                             std::optional<double> shaft_u1_global = std::nullopt,
                             Rostack_inspectable_physics inspectable_physics = nullptr,
                             rostack_detail::StageTrace* trace = nullptr) const
    {
        if (!physics && !inspectable_physics) {
            throw std::logic_error("No rotor physics function supplied");
        }
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
        if (trace) {
            constexpr std::array<const char*, NUM_ROSTACK_VARS> labels{
                "outlet cth2", "inlet u1", "inlet cm1", "k", "lambda",
                "sigma", "nu", "inlet cth1"
            };
            for (std::size_t i = 0; i < global.size(); ++i) {
                trace->design[labels[i]] = global[i];
            }
            trace->design["outlet u2"] = global[1] * global[4];
            trace->design["outlet cm2"] = global[2] * global[3];
            trace->performance["inlet cumulative tau"] = tau_in;
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
        if (!std::isfinite(local.cth2) || !std::isfinite(local.u1) ||
            !std::isfinite(local.cm1) || !std::isfinite(local.cth1)) {
            throw std::domain_error("Non-finite locally normalized velocity");
        }
        const auto p = inspectable_physics
            ? inspectable_physics(local, gamma,
                  trace ? Rostack_diagnostics{trace->performance} : Rostack_diagnostics{})
            : physics(local, gamma);
        if (trace) {
            // Canonical outputs take precedence over custom diagnostics using
            // the same names, so displayed constrained metrics stay consistent.
            trace->performance["tau"] = p.tau;
            trace->performance["M1 relative"] = p.mach_rel_in;
            trace->performance["M2 relative"] = p.mach_rel_out;
            trace->performance["M2 absolute"] = p.mach_abs_out;
            trace->performance["diffusion factor"] = p.diffusion;
            trace->performance["outlet cumulative tau"] = tau_in * p.tau;
        }
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
    std::string name_;
    bool inherits_inlet_;
    bool is_stator_;
    bool inherits_speed_;
    std::size_t n_opt_vars_ = 0;
    friend class Rostack_spool;
};

// A configuration object containing consecutive components on one shaft.
// The first NON-STATOR supplies the independent u1 (fixed or optimized).
// All later non-stators have u1 overridden; stators always retain u1 = 0.
// Their lambda still changes the radius at the next component's inlet.
// Validation and the effective variable layout are finalized by add_spool().
class Rostack_spool {
public:
    void add_rotor(Rostack_bounds bounds, Rostack_live_vars live_vars,
                   std::string stage_name = {})
    {
        components_.push_back({std::move(bounds), live_vars, std::move(stage_name)});
    }

    void add_rotor(const Rostack_rotor& rotor)
    {
        components_.push_back({rotor.bounds_, rotor.live_vars_, rotor.name_});
    }

    std::size_t stage_count() const noexcept { return components_.size(); }

private:
    struct Component {
        Rostack_bounds bounds;
        Rostack_live_vars live_vars;
        std::string name;
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

    explicit Rostack_model(Rostack_inspectable_physics physics, double gamma = 1.4)
        : inspectable_physics_(physics), gamma_(gamma)
    {
        if (!inspectable_physics_ || !std::isfinite(gamma_) || gamma_ <= 1.0) {
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
        std::set<std::string> names;
        for (const auto& rotor : rotors_) names.insert(rotor.name());
        bool has_rotor = false;
        for (const auto& c : spool.components_) {
            const bool is_stator = !c.live_vars[1] && c.bounds[1][0] == 0.0;
            const std::string name = c.name.empty()
                ? "Spool " + std::to_string(spool_count() + 1) +
                  " Stage " + std::to_string(prepared.size() + 1)
                : c.name;
            if (!names.insert(name).second) {
                throw std::invalid_argument("Duplicate stage name: " + name);
            }
            prepared.emplace_back(c.bounds, c.live_vars, name,
                                  !rotors_.empty() || !prepared.empty(),
                                  has_rotor && !is_stator);
            if (!is_stator) has_rotor = true;
        }

        // Allocate before changing the logical model. The prepared rotors
        // contain noexcept-movable value members, so these moves cannot throw.
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
    void add_rotor(Rostack_bounds bounds, Rostack_live_vars live_vars,
                   std::string stage_name = {})
    {
        Rostack_spool spool;
        spool.add_rotor(std::move(bounds), live_vars, std::move(stage_name));
        add_spool(spool);
    }

    std::size_t stage_count() const noexcept { return rotors_.size(); }
    std::size_t spool_count() const noexcept { return spool_rows_.size() - 1; }
    std::size_t get_n_opt_vars() const noexcept { return offsets_.back(); }
    const std::string& stage_name(std::size_t stage) const
    { return rotors_.at(stage).name(); }

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
        return evaluate_impl(x, nullptr);
    }

    // A runtime bool cannot change a C++ return type. This overload always
    // returns a report; false skips all map construction. The one-argument
    // overload above retains the original StackPerformance return type.
    // A bare model cannot report optimization constraints: use the UDP's
    // evaluate(x, true) overload when you want all three maps and fitness.
    Rostack_inspection evaluate(std::span<const double> x, bool inspect_output) const
    {
        Rostack_inspection report;
        report.result = evaluate_impl(x, inspect_output ? &report : nullptr);
        report.valid = report.result.valid;
        report.error = report.result.error;
        if (inspect_output) rostack_detail::find_invalid_values(report);
        return report;
    }

private:
    StackPerformance evaluate_impl(std::span<const double> x,
                                   Rostack_inspection* report) const
    {
        if ((!physics_ && !inspectable_physics_) || rotors_.empty()) {
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
                std::optional<rostack_detail::StageTrace> trace;
                if (report) trace.emplace();
                const auto fail = [&](const std::exception& e) {
                    result.valid = false;
                    result.failed_stage = i;
                    result.error = rotors_[i].name() + ": " + e.what();
                    result.tau_total = std::numeric_limits<double>::quiet_NaN();
                };
                try {
                    const auto row = rotors_[i].evaluate(
                        x.subspan(offsets_[i], offsets_[i + 1] - offsets_[i]),
                        flow, physics_, gamma_, shaft_u1_global,
                        inspectable_physics_, trace ? &*trace : nullptr);
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
                } catch (const std::domain_error& e) {
                    fail(e);
                } catch (const std::overflow_error& e) {
                    fail(e);
                }
                if (report) {
                    const std::string prefix = rotors_[i].name() + ", ";
                    for (const auto& [key, value] : trace->design) {
                        rostack_detail::insert_unique(report->design, prefix + key, value);
                    }
                    for (const auto& [key, value] : trace->performance) {
                        rostack_detail::insert_unique(report->performance_data, prefix + key, value);
                    }
                }
                if (!result.valid) return result;
            }
        }
        result.tau_total = flow->tau_cumulative;
        if (report) report->performance_data["Stack, total tau"] = result.tau_total;
        return result;
    }

private:
    std::vector<Rostack_rotor> rotors_;
    std::vector<std::size_t> offsets_{0};
    std::vector<std::size_t> spool_rows_{0}; // Flat row ranges per spool.
    Rostack_physics physics_ = nullptr;
    Rostack_inspectable_physics inspectable_physics_ = nullptr;
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
                    Objective objective = {}, double invalid_penalty = 1e6)
        : model_(std::move(model)), objective_(std::move(objective)),
          bounds_(model_.get_bounds()), invalid_penalty_(invalid_penalty)
    {
        if (!std::isfinite(invalid_penalty_) || invalid_penalty_ <= 0.0) {
            throw std::invalid_argument("Invalid-candidate penalty must be finite and positive");
        }
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
        return fitness_from_result(result, nullptr);
    }

    // Use this on the UDP for design + constraints + performance. It evaluates
    // physics once. Constraint entries are signed SCALED RESIDUALS, not raw
    // measured values; each key indicates the constraint relation.
    Rostack_inspection evaluate(std::span<const double> x,
                                bool inspect_output = true) const
    {
        auto report = model_.evaluate(x, inspect_output);
        report.fitness = fitness_from_result(report.result, &report.error);
        report.valid = report.error.empty();
        if (inspect_output) {
            std::size_t index = 1;
            const auto append = [&](const StageConstraint& c) {
                const std::string key = model_.stage_name(c.stage) + ", " +
                    std::string(metric_name(c.metric)) + " [" +
                    (c.relation == Relation::equal ? "equal" :
                     c.relation == Relation::upper ? "upper" : "lower") + "]";
                rostack_detail::insert_unique(report.constraints, key,
                                               report.fitness[index++]);
            };
            for (const auto& c : equalities_) append(c);
            for (const auto& c : inequalities_) append(c);
            rostack_detail::find_invalid_values(report);
        }
        return report;
    }

private:
    std::vector<double> penalty() const
    {
        // Fitness dimension must remain fixed even when a candidate is invalid.
        return std::vector<double>(1 + get_nec() + get_nic(), invalid_penalty_);
    }

    std::vector<double> fitness_from_result(const StackPerformance& result,
                                            std::string* error) const
    {
        if (!result.valid) {
            if (error) *error = result.error;
            return penalty();
        }
        try {
            std::vector<double> f;
            f.reserve(1 + get_nec() + get_nic());
            f.push_back(objective_(result));
            for (const auto& c : equalities_) f.push_back(c.residual(result));
            for (const auto& c : inequalities_) f.push_back(c.residual(result));
            for (double value : f) {
                if (!std::isfinite(value)) {
                    throw std::domain_error("Non-finite objective or constraint residual");
                }
            }
            if (error) error->clear();
            return f;
        } catch (const std::domain_error& e) {
            if (error) *error = std::string("Objective/constraints: ") + e.what();
        } catch (const std::overflow_error& e) {
            if (error) *error = std::string("Objective/constraints: ") + e.what();
        }
        return penalty();
    }

    static std::string_view metric_name(StageMetric m) noexcept
    {
        if (m == &StagePerformance::tau) return "tau";
        if (m == &StagePerformance::mach_rel_in) return "M1 relative";
        if (m == &StagePerformance::mach_rel_out) return "M2 relative";
        if (m == &StagePerformance::mach_abs_out) return "M2 absolute";
        if (m == &StagePerformance::diffusion) return "diffusion factor";
        return "unknown metric"; // Rejected during construction.
    }

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
    double invalid_penalty_ = 1e6;
};

// Pretty-print numeric maps. Nonfinite diagnostics are displayed as null and
// listed afterwards; they remain available in the report for debugging. An
// ancillary diagnostic NaN alone does not change physical fitness validity.
inline void print_inspection(std::ostream& out, const Rostack_inspection& report)
{
    const auto previous_precision = out.precision();
    out << std::setprecision(17);
    if (!report.valid) out << "Invalid evaluation: " << report.error << '\n';
    const auto print_map = [&](const char* title, const Rostack_values& values) {
        out << title << ":\n{\n";
        std::size_t i = 0;
        for (const auto& [key, value] : values) {
            out << "    " << std::quoted(key) << ": ";
            if (std::isfinite(value)) out << value;
            else out << "null";
            if (++i != values.size()) out << ',';
            out << '\n';
        }
        out << "}\n";
    };
    print_map("Best design", report.design);
    print_map("Constraints", report.constraints);
    print_map("Performance Data", report.performance_data);
    if (!report.invalid_values.empty()) {
        out << "Nonfinite reported values:\n";
        for (const auto& key : report.invalid_values) out << "    " << key << '\n';
    }
    out.precision(previous_precision);
}

// Failure policy: the rotor/physics may throw domain_error or overflow_error;
// model evaluation converts these to valid=false, and the UDP returns a finite
// penalty of the correct fitness dimension. Wrong vector sizes, configuration
// errors and other programming exceptions still throw. Choose the penalty to
// be unfavorable relative to the scale of your objective and constraints.
// This header supplies neither analytic derivatives nor Boost serialization.
