# theoretical_model/initial_model.py
#
# Pipeline overview:
# 1. Set primitives + parameter values.
# 2. Define clinical gate: Q(e) and Q'(e) under a_i = log(e_i) + eps_i.
# 3. Define payoffs by state:
#       - Monopoly: i passes, j fails  -> earns M * (p_mono - E[c] - tau(M))
#       - Duopoly: both pass           -> first-price procurement auction (lowest bid wins)
#       - Failure: i fails             -> 0
# 4. Solve symmetric interior effort e*(M) using a bracketing search + Brent.
# 5. Plot e*(M) over a grid of M and (optionally) save outputs for Quarto.

# --- libraries --- #
import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import brentq
from scipy.stats import norm

# --- parameters --- #

# effort cost: C(e) = 0.5 * k * e^2
k = 2.0

# clinical noise + threshold
sigma = 1.0
a_bar = 0.5

# monopoly price (or cap) used in monopoly state
p_mono = 2.5

# statutory wedge tau(M): increasing in market size
alpha = 0.008


def tau(M: float) -> float:
    # increasing wedge with diminishing effect
    return alpha * np.log(1.0 + M)


# marginal cost support in the auction stage
c_low = 0.0
c_high = 1.0

# numerical guard for log(e)
_EPS = 1e-12


# --- clinical gate: Q(e) and Q'(e) --- #

def Q(e: float) -> float:
    # Pr(log(e) + eps >= a_bar), eps ~ N(0, sigma^2)
    e = max(e, _EPS)
    z = (np.log(e) - a_bar) / sigma
    return float(norm.cdf(z))


def Q_prime(e: float) -> float:
    # derivative dQ/de = phi(z) * (1/(e*sigma))
    e = max(e, _EPS)
    z = (np.log(e) - a_bar) / sigma
    return float(norm.pdf(z) * (1.0 / (e * sigma)))


# --- payoffs by state --- #

def V_monopoly(M: float) -> float:
    # Ex-ante monopoly value:
    #   V_M(M) = M * (p_mono - E[c] - tau(M))
    avg_cost = 0.5 * (c_low + c_high)
    margin = p_mono - avg_cost - tau(M)
    return M * margin


def duopoly_profit_per_unit_two_bidders_uniform() -> float:
    # With two bidders, costs iid Uniform[c_low, c_high], first-price procurement auction:
    # expected profit per unit market size is (c_high - c_low) / 6.
    #
    # IMPORTANT NOTE ABOUT tau(M):
    # If the wedge enters payoff as: pi = M(b - c - tau(M)),
    # then in a procurement auction bidders optimally shift bids upward by tau(M),
    # so expected *markup* is unchanged under no binding price cap.
    #
    # If you want margin compression in the duopoly state, impose a cap/constraint
    # (e.g. bids cannot exceed some reference price), or let tau affect feasibility.
    return (c_high - c_low) / 6.0


def V_duopoly(M: float) -> float:
    # Expected value in the "both pass" state (auction).
    return M * duopoly_profit_per_unit_two_bidders_uniform()


# --- symmetric effort condition --- #

def foc_symmetric(e: float, M: float) -> float:
    # Symmetric FOC used in your slides:
    #
    # q  = Q(e)
    # qp = Q'(e)
    #
    # Expected payoff conditional on i passing the clinical gate:
    #   (1 - q) * V_M(M)   [rival fails -> monopoly]
    # + q       * V_D(M)   [rival passes -> auction]
    #
    # FOC:
    #   qp * [ (1-q) V_M + q V_D ] - k e = 0
    q = Q(e)
    qp = Q_prime(e)

    VM = V_monopoly(M)
    VD = V_duopoly(M)

    expected_value_given_pass = (1.0 - q) * VM + q * VD
    return qp * expected_value_given_pass - k * e


def solve_e_star(
    M: float,
    e_min: float = 1e-8,
    e_max: float = 100.0,
    n_grid: int = 250
) -> float:
    # Solve foc_symmetric(e, M) = 0 by:
    # 1) scanning a grid for a sign change
    # 2) bracketing root with brentq
    #
    # Conservative corner handling:
    # - if no sign change is found, return 0.0

    # Optional: if monopoly value <= 0, set effort to 0
    # (you can remove this if you want effort driven purely by the duopoly state)
    if V_monopoly(M) <= 0:
        return 0.0

    grid = np.geomspace(e_min, e_max, n_grid)
    vals = np.array([foc_symmetric(float(e), M) for e in grid])

    sgn = np.sign(vals)
    idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]

    if idx.size == 0:
        return 0.0

    a = float(grid[idx[0]])
    b = float(grid[idx[0] + 1])

    try:
        return float(brentq(lambda x: foc_symmetric(float(x), M), a, b))
    except ValueError:
        return 0.0


# --- comparative statics --- #

def run_comparative_statics(
    M_min: float = 0.1,
    M_max: float = 15.0,
    n_M: int = 120
) -> tuple[np.ndarray, np.ndarray]:
    M_values = np.linspace(M_min, M_max, n_M)
    e_stars = np.array([solve_e_star(float(M)) for M in M_values])
    return M_values, e_stars


# --- outputs --- #

def ensure_dir(path: str) -> None:
    if path == "":
        return
    os.makedirs(path, exist_ok=True)


def main() -> None:
    # run
    M_values, e_stars = run_comparative_statics()

    # diagnostics
    print("Comparative statics complete.")
    print("M range:", float(M_values.min()), "to", float(M_values.max()))
    print("e*(M) min/max:", float(e_stars.min()), float(e_stars.max()))
    print("Share of zeros:", float(np.mean(e_stars == 0.0)))

    # plot
    plt.figure(figsize=(9, 5))
    plt.plot(M_values, e_stars, linewidth=2)
    plt.xlabel("Market size (M)")
    plt.ylabel(r"Optimal effort $e^*(M)$")
    plt.title("Comparative statics: effort vs market size")
    plt.grid(alpha=0.3)

    # save (optional but useful for Quarto)
    out_dir = "theoretical_model/outputs"
    ensure_dir(out_dir)

    fig_path = os.path.join(out_dir, "effort_vs_market_size.png")
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    print("Saved figure:", fig_path)

    # also save the data so you can plot in Quarto/R if you want
    csv_path = os.path.join(out_dir, "effort_vs_market_size.csv")
    np.savetxt(
        csv_path,
        np.column_stack([M_values, e_stars]),
        delimiter=",",
        header="M,e_star",
        comments=""
    )
    print("Saved data:", csv_path)

    plt.show()


if __name__ == "__main__":
    main()