#pragma once

#include "rostack.hpp"

// Tip velocities use the row-inlet total-temperature normalization.
// rv.nu is the hub/tip radius ratio at BOTH stations (0 < nu <= 1).
// Meridional velocity and inlet total temperature are uniform across the span.
// Free-vortex swirl: ctheta_hub = ctheta_tip / nu.
// Constant spanwise chord: rv.sigma is TIP solidity; sigma_hub = sigma_tip / nu.
// Returns tip Mach numbers and HUB diffusion factor. Applies to stators too.
// Throws std::domain_error for invalid states; the Rostack UDP handles penalties.
// Compile study_physics.cpp once into the library or executable using this header.
StagePerformance study_physics(const RotorVariables& rv, double gamma,
                               Rostack_diagnostics performance_params);
