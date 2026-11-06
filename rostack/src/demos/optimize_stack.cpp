#include "rostack.hpp"
#include "rostack_demo_physics.hpp"
#include "rostack_optimizer.hpp"

#include <iostream>
#include <utility>
#include <vector>

int main()
{
    // Bounds/live order: cth2, u1, cm1, k, lambda, sigma, nu, cth1.
    // First-spool stator stays stationary; second-spool rotor has independent speed.
    Rostack_model stack{rostack_demo_physics, 1.4};
    Rostack_spool first_spool;
    first_spool.add_rotor({{{0.0, 1.5}, {1e-6, 1.5}, {0.4, 0.4}, {1.0, 1.0}, {1.0, 1.0}, {1.5, 1.5}, {0.5, 0.5}, {0.0, 0.0}}}, {true, true, false, false, false, false, false, false}, "Stage 1 compressor");
    first_spool.add_rotor({{{0.0, 1.5}, {0.0, 0.0}, {0.4, 0.4}, {1.0, 1.0}, {1.0, 1.0}, {1.5, 1.5}, {0.5, 0.5}, {0.0, 0.0}}}, {true, false, false, false, false, false, false, false}, "Interstage diffuser");
    stack.add_spool(first_spool);
    Rostack_spool second_spool;
    second_spool.add_rotor({{{-1.5, 1.5}, {-1.5, 1.5}, {0.4, 0.4}, {1.0, 1.0}, {1.0, 1.0}, {1.5, 1.5}, {0.5, 0.5}, {0.0, 0.0}}}, {true, true, false, false, false, false, false, false}, "Stage 2 compressor");
    stack.add_spool(second_spool);

    // Global component indices: 0 = rotor 1, 1 = diffuser, 2 = rotor 2.
    // Rotor relative Mach <= 1.2; all absolute Mach <= 0.8; diffusion <= 0.55.
    // Diffuser relative Mach equals absolute Mach because its blade speed is zero.
    std::vector<StageConstraint> constraints;
    constraints.push_back(upper(0, &StagePerformance::diffusion, 0.55));
    constraints.push_back(upper(0, &StagePerformance::mach_rel_in, 1.2));
    constraints.push_back(upper(0, &StagePerformance::mach_rel_out, 1.2));
    constraints.push_back(upper(0, &StagePerformance::mach_abs_out, 0.8));
    constraints.push_back(upper(1, &StagePerformance::diffusion, 0.55));
    constraints.push_back(upper(1, &StagePerformance::mach_rel_in, 0.8));
    constraints.push_back(upper(1, &StagePerformance::mach_rel_out, 0.8));
    constraints.push_back(upper(1, &StagePerformance::mach_abs_out, 0.8));
    constraints.push_back(upper(2, &StagePerformance::diffusion, 0.55));
    constraints.push_back(upper(2, &StagePerformance::mach_rel_in, 1.2));
    constraints.push_back(upper(2, &StagePerformance::mach_rel_out, 1.2));
    constraints.push_back(upper(2, &StagePerformance::mach_abs_out, 0.8));

    Rostack_problem<> udp{std::move(stack), std::move(constraints)};
    const auto best = optimize_udp(udp);
    print_inspection(std::cout, udp.evaluate(best.x, true));
}