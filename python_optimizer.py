"""Two-rotor co/counter-rotating compressor study with common loading limits.

Install: python -m pip install numpy scipy matplotlib
Run:     python compressor_mach_sweep.py

The input is an ALLOWABLE rotor-relative inlet tip Mach, not a prescribed
pre-shock Mach. Speed and swirl are optimized for each input. All calculations
refer to one representative constant-radius streamline with constant Cx.
Rotor shock loss uses the mean normal-shock entropy rise at inlet Mach M and
an estimated suction-surface impingement Mach k*M. Default k=1.15 is a rough
design assumption, not a calibrated correlation. This represents suction-
surface acceleration, not boundary-layer suction/aspiration. Separation,
profile, clearance, endwall and mixing losses are omitted. Attached flow is
an assumption. No independent structural speed cap is imposed.

Sources: Euler work/ideal-gas energy and entropy relations; Lieblein diffusion
factor; https://www.grc.nasa.gov/www/k-12/airplane/normal.html ; Kirchner,
Aerodynamic Design of an Aspirated Counter-Rotating Compressor (2002),
section 2.3.2 (mean entropy method; the default multiplier here is user chosen).
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


def shock_pt_recovery(mach, gamma=1.4):
    """Normal-shock pt_out/pt_in; proxy recovery = 1 for Mach <= 1.

    The subsonic extension means 'no inlet-normal-shock loss in this proxy',
    not that subsonic-inlet blades cannot contain local supersonic pockets.
    """
    m = np.asarray(mach, dtype=float)
    if not np.isfinite(gamma) or gamma <= 1:
        raise ValueError("gamma must be finite and greater than 1.")
    if np.any(~np.isfinite(m)) or np.any(m < 0):
        raise ValueError("Mach numbers must be finite and nonnegative.")
    m2 = np.maximum(m, 1.0)**2
    recovery = (
        ((gamma + 1)*m2 / ((gamma - 1)*m2 + 2))**(gamma/(gamma - 1))
        * ((gamma + 1)/(2*gamma*m2 - (gamma - 1)))**(1/(gamma - 1))
    )
    return float(recovery) if m.ndim == 0 else recovery


def _diffusion(v_in, v_out, turning, solidity):
    return 1 - v_out/v_in + abs(turning)/(2*solidity*v_in)


def rotor_shock_pt_recovery(mach, gamma=1.4, shock_mach_multiplier=1.15):
    """Kirchner-style mean entropy at M and M_imp = k*M.

    Since ds = -R*ln(pt_recovery), averaging entropy gives the geometric
    mean of the two recoveries, NOT a mean Mach or two successive shocks:
        recovery_row = exp(0.5*(ln(recovery(M)) + ln(recovery(k*M)))).
    k=1 exactly recovers the previous inlet-Mach-only rotor proxy.
    """
    if not np.isfinite(shock_mach_multiplier) or shock_mach_multiplier < 1:
        raise ValueError("shock_mach_multiplier must be finite and >= 1.")
    m = np.asarray(mach, dtype=float)
    inlet_recovery = shock_pt_recovery(m, gamma)
    impingement_recovery = shock_pt_recovery(shock_mach_multiplier*m, gamma)
    recovery = np.exp(0.5*(np.log(inlet_recovery) + np.log(impingement_recovery)))
    return float(recovery) if m.ndim == 0 else recovery


@lru_cache(maxsize=2048)
def _solve(mach_limit, architecture, Cx, T1, gamma, R, solidity, Dmax,
           equal_speeds, objective, shock_mach_multiplier):
    """Optimize one operating limit; cached to reuse reference solutions."""
    cp = gamma*R/(gamma - 1)
    Tt1 = T1 + Cx**2/(2*cp)
    co = architecture == "co"
    U1_bound = np.sqrt(mach_limit**2*gamma*R*T1 - Cx**2)
    # Since S1 <= U1, this also bounds all feasible second-rotor speeds.
    hottest_T2 = T1 + U1_bound**2/cp
    U2_bound = np.sqrt(mach_limit**2*gamma*R*hottest_T2 - Cx**2)
    scale = max(U1_bound, U2_bound)
    tiny = 1e-8
    if equal_speeds:
        bounds = [(tiny, U1_bound/scale), (0, U1_bound/scale)]
        if co:
            bounds.append((0, U1_bound/scale))
    else:
        bounds = [(tiny, U1_bound/scale), (tiny, U2_bound/scale),
                  (0, U1_bound/scale)]
        if co:
            bounds.append((0, U2_bound/scale))

    def evaluate(x):
        values = np.asarray(x)*scale
        if equal_speeds:
            U1 = U2 = values[0]
            S1 = values[1]
            S2 = values[2] if co else S1
        else:
            U1, U2, S1 = values[:3]
            S2 = values[3] if co else S1

        # U2 is a magnitude; the counter-rotating rotor has signed speed -U2.
        w1, w2 = U1*S1, U2*S2  # J/kg
        T_R1out = T1 + (w1 - S1**2/2)/cp
        T_R2in = T1 + w1/cp if co else T_R1out
        Ttout = Tt1 + (w1 + w2)/cp
        T_R2out = Ttout - (Cx**2 + (S2**2 if co else 0))/(2*cp)

        W1in = np.hypot(Cx, U1)
        W1out = np.hypot(Cx, S1 - U1)
        W2in = np.hypot(Cx, U2 if co else S1 + U2)
        W2out = np.hypot(Cx, S2 - U2 if co else U2)
        M1 = W1in/np.sqrt(gamma*R*T1)
        # Clipping only keeps infeasible trial points numerically evaluable.
        M2 = W2in/np.sqrt(gamma*R*max(T_R2in, tiny))
        D1 = _diffusion(W1in, W1out, S1, solidity)
        D2 = _diffusion(W2in, W2out, S2, solidity)
        DS1 = _diffusion(np.hypot(Cx, S1), Cx, S1, solidity) if co else np.nan
        DS2 = _diffusion(np.hypot(Cx, S2), Cx, S2, solidity) if co else np.nan

        margins = [mach_limit - M1, mach_limit - M2,
                   Dmax - D1, Dmax - D2, (U1 - S1)/scale,
                   (T_R1out - tiny)/T1, (T_R2out - tiny)/T1]
        if co:
            margins.extend([Dmax - DS1, Dmax - DS2, (U2 - S2)/scale])

        rotor_recovery = rotor_shock_pt_recovery(
            [M1, M2], gamma, shock_mach_multiplier)
        if co:
            MS1 = np.hypot(Cx, S1)/np.sqrt(gamma*R*max(T_R1out, tiny))
            MS2 = np.hypot(Cx, S2)/np.sqrt(gamma*R*max(T_R2out, tiny))
            stator_recovery = shock_pt_recovery([MS1, MS2], gamma)
        else:
            MS1 = MS2 = np.nan
            stator_recovery = np.ones(2)

        # Add shock entropy rises, then apply them in the absolute frame.
        log_recovery = np.log(rotor_recovery).sum() + np.log(stator_recovery).sum()
        log_pr_zero = gamma/(gamma - 1)*np.log(Ttout/Tt1)
        log_pr_shock = log_pr_zero + log_recovery
        result = {
            "mach_limit": mach_limit,
            "rotor1_speed_m_s": U1,
            "rotor2_speed_m_s": U2,
            "rotor1_exit_swirl_m_s": S1,
            "rotor2_exit_swirl_m_s": S2 if co else 0.0,
            "rotor1_inlet_relative_mach": M1,
            "rotor2_inlet_relative_mach": M2,
            "shock_mach_multiplier": shock_mach_multiplier,
            "rotor1_impingement_mach": shock_mach_multiplier*M1,
            "rotor2_impingement_mach": shock_mach_multiplier*M2,
            "rotor1_shock_recovery": rotor_recovery[0],
            "rotor2_shock_recovery": rotor_recovery[1],
            "rotor1_shock_entropy_J_kgK": -R*np.log(rotor_recovery[0]),
            "rotor2_shock_entropy_J_kgK": -R*np.log(rotor_recovery[1]),
            "rotor1_diffusion_factor": D1,
            "rotor2_diffusion_factor": D2,
            "stator1_diffusion_factor": DS1,
            "stator2_diffusion_factor": DS2,
            "stator1_inlet_mach": MS1,
            "stator2_inlet_mach": MS2,
            "rotor1_work_kj_kg": w1/1000,
            "rotor2_work_kj_kg": w2/1000,
            "total_work_kj_kg": (w1 + w2)/1000,
            "pressure_ratio_zero_losses": np.exp(log_pr_zero),
            "pressure_ratio_shock_proxy": np.exp(log_pr_shock),
            "shock_recovery_product": np.exp(log_recovery),
            "outlet_total_temperature_K": Ttout,
            "minimum_constraint_margin": min(margins),
        }
        return result, np.asarray(margins), log_pr_shock

    def loss(x):
        result, _, log_pr = evaluate(x)
        return -result["total_work_kj_kg"]/100 if objective == "work" else -log_pr

    constraints = {"type": "ineq", "fun": lambda x: evaluate(x)[1]}
    rng = np.random.default_rng(145)
    starts = []
    for _ in range(12):
        U1 = rng.uniform(.35, .95)*U1_bound
        U2 = U1 if equal_speeds else rng.uniform(.25, .95)*U2_bound
        S1 = rng.uniform(.10, .55)*U1
        S2 = rng.uniform(.10, .55)*U2
        values = [U1, S1] if equal_speeds else [U1, U2, S1]
        if co:
            values.append(S2)
        starts.append(np.asarray(values)/scale)

    best = None
    for start in starts:
        candidate = minimize(loss, start, method="SLSQP", bounds=bounds,
                             constraints=constraints,
                             options={"ftol": 1e-10, "maxiter": 500})
        if candidate.success and min(evaluate(candidate.x)[1]) >= -1e-7:
            if best is None or candidate.fun < best.fun:
                best = candidate
    if best is None:
        raise RuntimeError(f"No converged feasible solution: {architecture}, M={mach_limit}.")
    return {key: float(value) for key, value in evaluate(best.x)[0].items()}


def compressor_performance(mach_limit, architecture="counter", *,
                           reference_mach=1.2, Cx=150.0, T1=288.0,
                           gamma=1.4, R=287.0, solidity=1.5, Dmax=0.50,
                           equal_speeds=False, objective="work",
                           shock_mach_multiplier=1.15):
    """Optimize the velocity triangles and return performance as a dictionary.

    mach_limit: scalar or 1-D array of allowed rotor-relative inlet tip Mach.
    architecture: 'co' or 'counter'. Both have axial inlet/final outlet.
    reference_mach: reference design for all 'additional_*' work outputs.
    Cx [m/s], T1 [K]: axial speed and inlet STATIC temperature.
    solidity, Dmax: common values for every existing rotor/stator row.
    equal_speeds: if True, impose U1 = |U2|; otherwise vary independently.
    shock_mach_multiplier: k >= 1 in estimated rotor M_imp = k*M_rel,in.
        Rotor entropy is the mean of normal-shock entropy rises at inlet
        and impingement Mach. Default 1.15 is a user-chosen first estimate;
        set 1.0 for the original model. It changes shock losses, not the
        inlet Mach constraint or the velocity-triangle equations. With
        objective='shock_pressure', reoptimization may select new triangles.
        Stators retain the original absolute-inlet-Mach shock proxy.
    objective: 'work' maximizes total Euler work; 'shock_pressure' maximizes
        pressure ratio under the selected rotor shock proxy. Both pressure
        ratios are evaluated on the SAME selected velocity triangles.

    Speeds are positive magnitudes; work values are kJ/kg. Additional work
    is each selected rotor's work minus its work at reference_mach, NOT a
    separate maximization of each rotor. Missing stator quantities are NaN.
    Array input gives arrays of the same shape for every returned value.
    Profile, separation, clearance, mixing, and endwall losses are omitted.
    Numerical multistart optimization does not prove global optimality.
    """
    if architecture not in ("co", "counter"):
        raise ValueError("architecture must be 'co' or 'counter'.")
    if objective not in ("work", "shock_pressure"):
        raise ValueError("objective must be 'work' or 'shock_pressure'.")
    if not isinstance(equal_speeds, (bool, np.bool_)):
        raise ValueError("equal_speeds must be True or False.")
    settings = np.asarray([Cx, T1, gamma, R, solidity, Dmax, reference_mach,
                           shock_mach_multiplier], dtype=float)
    if np.any(~np.isfinite(settings)):
        raise ValueError("All settings must be finite.")
    if Cx <= 0 or T1 <= 0 or R <= 0 or gamma <= 1 or solidity <= 0 or not 0 < Dmax < 1:
        raise ValueError("Require Cx,T1,R,solidity > 0, gamma > 1, and 0 < Dmax < 1.")
    if shock_mach_multiplier < 1:
        raise ValueError("shock_mach_multiplier must be >= 1.")
    mach = np.asarray(mach_limit, dtype=float)
    axial_mach = Cx/np.sqrt(gamma*R*T1)
    if mach.ndim > 1 or mach.size == 0 or np.any(~np.isfinite(mach)):
        raise ValueError("mach_limit must be a finite scalar or nonempty 1-D array.")
    if np.any(mach <= axial_mach) or reference_mach <= axial_mach:
        raise ValueError("Mach limits must exceed the inlet axial Mach number.")
    args = (architecture, float(Cx), float(T1), float(gamma), float(R),
            float(solidity), float(Dmax), bool(equal_speeds), objective,
            float(shock_mach_multiplier))
    reference = _solve(float(reference_mach), *args)
    rows = []
    for value in mach.reshape(-1):
        result = dict(_solve(float(value), *args))
        for name in ("rotor1", "rotor2", "total"):
            key = f"{name}_work_kj_kg"
            result[f"reference_{key}"] = reference[key]
            result[f"additional_{key}"] = result[key] - reference[key]
        result["reference_mach"] = float(reference_mach)
        rows.append(result)
    if mach.ndim == 0:
        return rows[0]
    return {key: np.array([row[key] for row in rows]) for key in rows[0]}


def plot_mach_sweep(mach_limits=None, *, show=True, save_path=None, **settings):
    """Plot both architectures; return (figure, {'co': results, 'counter': results})."""
    import matplotlib.pyplot as plt
    if mach_limits is None:
        mach_limits = np.linspace(1.2, 1.5, 61)
    mach_limits = np.asarray(mach_limits, dtype=float)
    if mach_limits.ndim != 1 or len(mach_limits) < 2:
        raise ValueError("Plotting requires at least two Mach limits in a 1-D array.")
    if np.any(np.diff(mach_limits) <= 0):
        raise ValueError("Plotting Mach limits must be strictly increasing.")
    data = {arch: compressor_performance(mach_limits, arch, **settings)
            for arch in ("co", "counter")}
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = {"co": "#2455a4", "counter": "#ba5716"}
    labels = {"co": "Co-rotating", "counter": "Counter-rotating"}
    for arch, result in data.items():
        color, label = colors[arch], labels[arch]
        axes[0, 0].plot(mach_limits, result["pressure_ratio_zero_losses"],
                        "--", color=color, lw=2, label=f"{label}: zero losses")
        axes[0, 0].plot(mach_limits, result["pressure_ratio_shock_proxy"],
                        color=color, lw=2.2, label=f"{label}: shock proxy")
        axes[0, 1].plot(mach_limits, result["total_work_kj_kg"],
                        color=color, lw=2.2, label=label)
        axes[1, 0].plot(mach_limits, result["additional_rotor1_work_kj_kg"],
                        color=color, lw=2.2, label=label)
        axes[1, 1].plot(mach_limits, result["additional_rotor2_work_kj_kg"],
                        color=color, lw=2.2, label=label)
    titles = ["Overall total-pressure ratio", "Total Euler work capacity",
              "Additional rotor 1 work", "Additional rotor 2 work"]
    reference = settings.get("reference_mach", 1.2)
    ylabels = ["Pressure ratio", "Work [kJ/kg]",
               f"Work increase from M = {reference:g} [kJ/kg]",
               f"Work increase from M = {reference:g} [kJ/kg]"]
    for ax, title, ylabel in zip(axes.flat, titles, ylabels):
        ax.set(title=title, xlabel="Allowed rotor-relative inlet tip Mach", ylabel=ylabel,
               xlim=(mach_limits[0], mach_limits[-1]))
        ax.grid(True, alpha=.23)
        ax.legend(fontsize=8, frameon=False)
    speed_label = "equal speed magnitudes" if settings.get("equal_speeds", False) else "independent rotor speeds"
    fig.suptitle(
        f"Common loading: solidity = {settings.get('solidity', 1.5):g}, "
        f"D limit = {settings.get('Dmax', .5):g}; {speed_label}\n"
        f"Rotor shock model: mean entropy at M and "
        f"{settings.get('shock_mach_multiplier', 1.15):g} M", fontsize=12)
    if save_path is not None:
        fig.savefig(save_path, dpi=200)
    if show:
        plt.show()
    return fig, data


if __name__ == "__main__":
    figure, results = plot_mach_sweep(
        save_path=Path(__file__).with_name("compressor_mach_sweep.png"))
    for arch, result in results.items():
        print(f"\n{arch}: M={result['mach_limit'][0]:.2f} to {result['mach_limit'][-1]:.2f}")
        for key in ("pressure_ratio_zero_losses", "pressure_ratio_shock_proxy",
                    "rotor1_work_kj_kg", "rotor2_work_kj_kg",
                    "additional_rotor1_work_kj_kg", "additional_rotor2_work_kj_kg"):
            print(f"  {key}: {result[key][0]:.4f} -> {result[key][-1]:.4f}")
