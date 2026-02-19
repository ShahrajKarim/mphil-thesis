# theoretical_model/initial_model.py
#
# Comparative statics: symmetric Nash effort e*(M)
# Two firms choose effort -> clinical gate -> if both pass, price stage (procurement auction)
#
# Units:
# - M is in millions of prescriptions/patients (volume)
# - prices/costs/tau(M) are $ per prescription
# - profits are in million-$
# - k is in million-$ units

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.stats import norm


_EPS_EFFORT = 1e-10


# ============================================================
# 1) Parameters + wedge shapes
# ============================================================

def tau_log1p(M: float) -> float:
    return float(np.log1p(M))


def tau_sqrt(M: float) -> float:
    return float(np.sqrt(M))


def tau_linear(M: float) -> float:
    return float(M)


def tau_saturating(M: float) -> float:
    return float(M / (1.0 + M))


TAU_SHAPES: Dict[str, Callable[[float], float]] = {
    "log1p(M)": tau_log1p,
    "sqrt(M)": tau_sqrt,
    "M": tau_linear,
    "M/(1+M)": tau_saturating,
}


@dataclass(frozen=True)
class Params:
    # effort cost: 0.5 * k * e^2 (million-$)
    k: float = 0.5

    # clinical tech: a = log(e) + eps, eps ~ N(0, sigma^2)
    sigma: float = 1.0
    a_bar: float = 0.0

    # monopoly price cap ($/rx)
    p_mono: float = 3.0

    # cost distribution in price stage ($/rx)
    c_low: float = 0.0
    c_high: float = 1.0

    # wedge intensity (per-unit wedge in $/rx)
    alpha: float = 0.25

    # wedge shape: tau(M) = alpha * tau_shape(M)
    tau_shape: Callable[[float], float] = tau_log1p


def tau(M_millions: float, p: Params) -> float:
    return float(p.alpha * p.tau_shape(float(M_millions)))


# ============================================================
# 2) Clinical pass probability
# ============================================================

def Q(e: float, p: Params) -> float:
    e_safe = max(float(e), _EPS_EFFORT)
    z = (np.log(e_safe) - p.a_bar) / p.sigma
    return float(norm.cdf(z))


# ============================================================
# 3) Payoffs
# ============================================================

def expected_cost_uniform(p: Params) -> float:
    return 0.5 * (p.c_low + p.c_high)


def duopoly_markup_two_bidders_uniform(p: Params) -> float:
    # Two-bidder procurement auction: expected winner markup (in $/rx)
    return (p.c_high - p.c_low) / 6.0


def V_monopoly(M_millions: float, p: Params) -> float:
    # profit (million-$) = M * (p_mono - E[c] - tau(M))
    margin = p.p_mono - expected_cost_uniform(p) - tau(M_millions, p)
    return float(M_millions * margin)


def V_duopoly(M_millions: float, p: Params) -> float:
    """
    Key change (for a visible break):
    duopoly expected per-unit profit is the auction markup minus the wedge.
    Once tau(M) is big enough, the duopoly state becomes near-zero.
    """
    per_unit = duopoly_markup_two_bidders_uniform(p) - tau(M_millions, p)
    per_unit = max(0.0, per_unit)
    return float(M_millions * per_unit)


# ============================================================
# 4) Objective, best response, symmetric fixed point
# ============================================================

def psi(e_self: float, M_millions: float, q_rival: float, p: Params) -> float:
    e_safe = max(float(e_self), _EPS_EFFORT)
    q_self = Q(e_safe, p)

    VM = V_monopoly(M_millions, p)
    VD = V_duopoly(M_millions, p)

    expected_payoff = q_self * (1.0 - q_rival) * VM + q_self * q_rival * VD
    effort_cost = 0.5 * p.k * (e_safe ** 2)

    return float(expected_payoff - effort_cost)


def argmax_psi_given_qrival(
    M_millions: float,
    q_rival: float,
    p: Params,
    e_bounds: Tuple[float, float] = (1e-8, 30.0),
) -> float:
    lo, hi = e_bounds

    res = minimize_scalar(
        lambda e: -psi(e, M_millions, q_rival, p),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": 1e-6},
    )

    e_hat = float(res.x)

    if psi(e_hat, M_millions, q_rival, p) < 0.0:
        return 0.0

    return e_hat


def best_response(
    M_millions: float,
    e_rival: float,
    p: Params,
    e_bounds: Tuple[float, float] = (1e-8, 30.0),
) -> float:
    q_rival = Q(e_rival, p)
    return argmax_psi_given_qrival(M_millions, q_rival, p, e_bounds=e_bounds)


def solve_symmetric_nash(
    M_millions: float,
    p: Params,
    e_init: float = 1.0,
    max_iter: int = 200,
    tol: float = 1e-5,
) -> float:
    e_old = max(float(e_init), 0.0)

    for _ in range(max_iter):
        e_new = best_response(M_millions, e_old, p)

        if abs(e_new - e_old) <= tol:
            return float(e_new)

        e_old = 0.8 * e_old + 0.2 * e_new

    return float(e_old)


# ============================================================
# 5) Comparative statics helper
# ============================================================

def solve_curve(M_grid: np.ndarray, p: Params, e_init: float = 1.0) -> np.ndarray:
    e_grid = np.zeros_like(M_grid, dtype=float)
    e_guess = float(e_init)

    for i, M in enumerate(M_grid):
        e_star = solve_symmetric_nash(float(M), p, e_init=e_guess)
        e_grid[i] = e_star
        e_guess = max(e_star, 1e-6)

    return e_grid


# ============================================================
# 6) Main: 2x2 robustness plots
# ============================================================

def main() -> None:
    # baseline chosen so a break is plausible within M in [0, 25]
    base = Params(
        k=0.6,
        sigma=1.0,
        a_bar=0.0,
        p_mono=3.0,      # smaller monopoly margin => wedge can matter
        c_low=0.0,
        c_high=1.0,
        alpha=0.25,      # stronger wedge
        tau_shape=tau_log1p,
    )

    M_grid = np.linspace(0.1, 25.0, 180)

    # A) alpha sweep (include values that cross the duopoly “kill point”)
    alpha_grid = [0.05, 0.15, 0.25, 0.40, 0.60]

    # B) sigma sweep
    sigma_grid = [0.6, 1.0, 1.6]

    # C) cost supports (changes the auction markup (c_high-c_low)/6)
    cost_supports = [(0.0, 0.6), (0.0, 1.0), (0.0, 1.6)]

    # D) wedge shapes
    tau_shape_grid = ["log1p(M)", "sqrt(M)", "M/(1+M)", "M"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    axA, axB, axC, axD = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # ---- Panel A: vary alpha ----
    for a in alpha_grid:
        p = Params(**{**base.__dict__, "alpha": a, "tau_shape": TAU_SHAPES["M"]})
        e_grid = solve_curve(M_grid, p)
        axA.plot(M_grid, e_grid, linewidth=2.2, label=f"alpha = {a:g} (tau=M)")

    axA.set_title("Vary rebate intensity (alpha)")
    axA.set_xlabel("Market size M (millions)")
    axA.set_ylabel("Optimal effort $e^*(M)$")
    axA.grid(alpha=0.3)
    axA.legend()

    # ---- Panel B: vary sigma ----
    for s in sigma_grid:
        p = Params(**{**base.__dict__, "sigma": s})
        e_grid = solve_curve(M_grid, p)
        axB.plot(M_grid, e_grid, linewidth=2.2, label=f"sigma = {s:g}")

    axB.set_title("Vary clinical noise (sigma)")
    axB.set_xlabel("Market size M (millions)")
    axB.set_ylabel("Optimal effort $e^*(M)$")
    axB.grid(alpha=0.3)
    axB.legend()

    # ---- Panel C: vary cost dispersion ----
    for (cl, ch) in cost_supports:
        p = Params(**{**base.__dict__, "c_low": cl, "c_high": ch})
        e_grid = solve_curve(M_grid, p)
        axC.plot(M_grid, e_grid, linewidth=2.2, label=f"cost ~ U[{cl:g}, {ch:g}]")

    axC.set_title("Vary cost dispersion (auction markup)")
    axC.set_xlabel("Market size M (millions)")
    axC.set_ylabel("Optimal effort $e^*(M)$")
    axC.grid(alpha=0.3)
    axC.legend()

    # ---- Panel D: vary rebate functional form ----
    for name in tau_shape_grid:
        p = Params(**{**base.__dict__, "tau_shape": TAU_SHAPES[name]})
        e_grid = solve_curve(M_grid, p)
        axD.plot(M_grid, e_grid, linewidth=2.2, label=name)

    axD.set_title("Vary rebate functional form")
    axD.set_xlabel("Market size M (millions)")
    axD.set_ylabel("Optimal effort $e^*(M)$")
    axD.grid(alpha=0.3)
    axD.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()