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

from typing import Dict, Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.stats import norm


# ============================================================
# 1) Parameters
# ============================================================

# effort cost: 0.5 * k * e^2 (million-$)
k = 0.5

# clinical tech: a = log(e) + eps, eps ~ N(0, sigma^2)
sigma = 1.0
a_bar = 0.0

# monopoly price cap ($ per prescription)
p_mono = 10.0

# cost distribution in price stage ($ per prescription)
c_low, c_high = 0.0, 1.0

# wedge intensity (per-unit wedge in $ per prescription)
alpha = 0.1

_EPS_EFFORT = 1e-10


def tau(M_millions: float) -> float:
    """
    Per-unit wedge tau(M) in $/rx.
    M is in millions, so log1p(M) avoids mechanical blow-ups from units.
    """
    return alpha * np.log1p(M_millions)


# ============================================================
# 2) Clinical pass probability
# ============================================================

def Q(e: float) -> float:
    e_safe = max(e, _EPS_EFFORT)
    z = (np.log(e_safe) - a_bar) / sigma
    return float(norm.cdf(z))


# ============================================================
# 3) Payoffs
# ============================================================

def expected_cost_uniform() -> float:
    return 0.5 * (c_low + c_high)


def V_monopoly(M_millions: float) -> float:
    # profit (million-$) = M * (p_mono - E[c] - tau(M))
    margin = p_mono - expected_cost_uniform() - tau(M_millions)
    return M_millions * margin


def profit_per_unit_duopoly_uniform_two_bidders() -> float:
    # Two-bidder procurement auction: expected winner markup = (c_high - c_low)/6
    return (c_high - c_low) / 6.0


def V_duopoly(M_millions: float) -> float:
    return M_millions * profit_per_unit_duopoly_uniform_two_bidders()


# ============================================================
# 4) Objective, best response, symmetric fixed point
# ============================================================

def psi(e_self: float, M_millions: float, q_rival: float) -> float:
    """
    Expected profit (million-$), given rival pass probability q_rival:

      Psi(e) = q*(1-q_rival)*V_M + q*q_rival*V_D - 0.5*k*e^2
    """
    e_safe = max(e_self, _EPS_EFFORT)
    q_self = Q(e_safe)

    VM = V_monopoly(M_millions)
    VD = V_duopoly(M_millions)

    expected_payoff = q_self * (1.0 - q_rival) * VM + q_self * q_rival * VD
    effort_cost = 0.5 * k * (e_safe ** 2)

    return float(expected_payoff - effort_cost)


def argmax_psi_given_qrival(
    M_millions: float,
    q_rival: float,
    e_bounds: Tuple[float, float] = (1e-8, 30.0)
) -> float:
    """
    Direct maximisation of Psi(e | q_rival fixed).
    Returns 0 if the maximised value is negative (exit).
    """
    lo, hi = e_bounds

    res = minimize_scalar(
        lambda e: -psi(e, M_millions, q_rival),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": 1e-6}
    )

    e_hat = float(res.x)

    # Exit option
    if psi(e_hat, M_millions, q_rival) < 0:
        return 0.0

    return e_hat


def best_response(
    M_millions: float,
    e_rival: float,
    e_bounds: Tuple[float, float] = (1e-8, 30.0)
) -> float:
    q_rival = Q(e_rival)
    return argmax_psi_given_qrival(M_millions, q_rival, e_bounds=e_bounds)


def solve_symmetric_nash(
    M_millions: float,
    e_init: float = 1.0,
    max_iter: int = 200,
    tol: float = 1e-5
) -> float:
    """
    Symmetric fixed point: e = BR(M, e).
    """
    e_old = max(e_init, 0.0)

    for _ in range(max_iter):
        e_new = best_response(M_millions, e_old)

        if abs(e_new - e_old) <= tol:
            return float(e_new)

        # mild damping for stability
        e_old = 0.8 * e_old + 0.2 * e_new

    return float(e_old)


# ============================================================
# 5) Comparative statics + plots
# ============================================================

def main() -> None:

    # ---- A) Solve e*(M) on a grid ----
    M_grid = np.linspace(0.1, 25.0, 160)

    e_grid = np.zeros_like(M_grid)
    e_guess = 1.0

    for i, M in enumerate(M_grid):
        e_star = solve_symmetric_nash(M, e_init=e_guess)
        e_grid[i] = e_star
        e_guess = max(e_star, 1e-6)

    # ---- B) Pick cases: small / peak / large (automatically) ----
    interior_idx = np.where(e_grid > 1e-6)[0]
    if len(interior_idx) == 0:
        raise RuntimeError("No interior solutions found. Reduce alpha or increase p_mono.")

    peak_i = interior_idx[np.argmax(e_grid[interior_idx])]
    M_peak = float(M_grid[peak_i])

    M_small = float(np.quantile(M_grid[interior_idx], 0.10))
    M_large = float(np.quantile(M_grid[interior_idx], 0.90))

    M_cases = [M_small, M_peak, M_large]

    # Solve equilibrium at the chosen cases
    case_e = {M: solve_symmetric_nash(M, e_init=1.0) for M in M_cases}

    # ---- C) Colors consistent across panels ----
    colors: Dict[float, str] = {
        M_small: "C0",
        M_peak:  "C2",
        M_large: "C3",
    }

    # ---- D) Build 2-panel figure ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: effort vs market size
    ax1 = axes[0]
    ax1.plot(M_grid, e_grid, linewidth=2.5, label="Optimal effort")

    for M in M_cases:
        ax1.scatter(M, case_e[M], s=90, color=colors[M], label=f"M = {M:.2f}")

    ax1.set_title("Effort vs Market Size")
    ax1.set_xlabel("Market size M (millions)")
    ax1.set_ylabel("Optimal effort $e^*(M)$")
    ax1.grid(alpha=0.3)
    ax1.legend()

    # Panel 2: profit hills Psi(e; M | q_rival fixed at equilibrium)
    ax2 = axes[1]

    e_max_plot = max(8.0, 1.25 * float(np.max(e_grid)))
    e_plot = np.linspace(0.0, e_max_plot, 300)

    for M in M_cases:
        e_star = case_e[M]
        q_rival = Q(e_star)

        psi_vals = [psi(max(e, _EPS_EFFORT), M, q_rival) for e in e_plot]
        ax2.plot(e_plot, psi_vals, linewidth=2.5, color=colors[M], label=f"M = {M:.2f}")

        # Mark the argmax of the plotted hill (should match equilibrium in symmetry)
        e_hat = argmax_psi_given_qrival(M, q_rival, e_bounds=(1e-8, 30.0))
        psi_hat = psi(max(e_hat, _EPS_EFFORT), M, q_rival)
        ax2.scatter(e_hat, psi_hat, s=90, color=colors[M])

        if abs(e_hat - e_star) > 1e-3:
            print(f"[warning] M={M:.2f}: equilibrium e={e_star:.4f} but hill argmax e={e_hat:.4f}")

    ax2.set_title(r"Profit function $\Psi(e; M)$")
    ax2.set_xlabel("Effort level (e)")
    ax2.set_ylabel("Expected profit (million-$)")
    ax2.grid(alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()