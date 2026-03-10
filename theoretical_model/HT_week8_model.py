"""
Comparative Statics: ACA Expansion & Pharmaceutical Innovation
================================================================
Model: HT Week 8 — Mohammad-Shah Karim

Four panels, each showing Power / Log / Linear tau forms overlaid.
Run with: python comparative_statics.py
Saves: comparative_statics.png in the same directory.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import norm
from scipy.optimize import brentq

# ─────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────
P = dict(
    sigma   = 1.0,   # Std dev of R&D noise ε_i. Higher σ → more uncertain innovation process → Q(e) less responsive to effort
    a_bar   = 0.5,   # Regulatory approval threshold ā. Drug passes P&T iff a_i ≥ ā. Higher ā → harder approval → lower Q(e) → lower e*
    p_bar   = 10.0,  # Price cap p̄ imposed by state on monopolist. Monopolist always sets p* = p̄ (binding constraint). Sets ceiling on V^M
    c_bar   = 2.0,   # Expected marginal production cost c̄. Enters V^M as M(p̄ - c̄ - τ) and sets lower bound of duopoly auction integral
    k_base  = 1.0,   # Effort cost curvature in C(e) = ½ke². Higher k → FOC satisfied at lower e*. Varies across dynamic states: k(V^M) < k(V^D) < k(0)
    M_base  = 100.0, # Baseline Medicaid market size (enrollee population). Central variable shifted by ACA expansion. Enters V^M and V^D as revenue multiplier
    alpha   = 0.02,  # Scale of rebate function τ(M). Sets baseline level of state monopsony power independently of market size
    beta    = 0.75,  # Elasticity of rebate w.r.t. M: d(ln τ)/d(ln M) = β. β < 1 → sublinear (volume effect dominates). β > 1 → incentive trap operates
    c_lo    = 0.5,   # Lower bound of uniform cost distribution F(c) in duopoly auction
    c_hi    = 3.5,   # Upper bound of uniform cost distribution F(c). Width [c_lo, c_hi] governs auction competitiveness and size of V^D
)

# ─────────────────────────────────────────────
# TAU FORMS
# ─────────────────────────────────────────────
TAU_FORMS = {
    r"Power ($\tau = \alpha M^{\beta}$)":  lambda M: P["alpha"] * M ** P["beta"],
    r"Log ($\tau = \alpha \ln M$)":        lambda M: P["alpha"] * np.log(np.maximum(M, 1e-6)),
    r"Linear ($\tau = \alpha M$)":         lambda M: P["alpha"] * M,
}

COLORS = {
    r"Power ($\tau = \alpha M^{\beta}$)":  "#2563EB",
    r"Log ($\tau = \alpha \ln M$)":        "#DC2626",
    r"Linear ($\tau = \alpha M$)":         "#16A34A",
}

DASH = {
    r"Power ($\tau = \alpha M^{\beta}$)":  "-",
    r"Log ($\tau = \alpha \ln M$)":        "--",
    r"Linear ($\tau = \alpha M$)":         "-.",
}

# ─────────────────────────────────────────────
# MARKET VALUES
# ─────────────────────────────────────────────
def V_mono(M, tau_fn):
    return M * (P["p_bar"] - P["c_bar"] - tau_fn(M))

def V_duo(M):
    integral = (P["c_hi"] - P["c_bar"])**2 / (2 * (P["c_hi"] - P["c_lo"]))
    return M * integral

# ─────────────────────────────────────────────
# Q(e) AND Q'(e)
# ─────────────────────────────────────────────
def Q(e):
    return 1.0 - norm.cdf((P["a_bar"] - np.log(np.maximum(e, 1e-9))) / P["sigma"])

def Q_prime(e):
    return norm.pdf((P["a_bar"] - np.log(np.maximum(e, 1e-9))) / P["sigma"]) / (P["sigma"] * np.maximum(e, 1e-9))

# ─────────────────────────────────────────────
# EQUILIBRIUM SOLVER (FIXED)
# ─────────────────────────────────────────────
def foc(e, VM, VD, k):
    return Q_prime(e) * ((1 - Q(e)) * VM + Q(e) * VD) - k * e

def solve_effort(VM, VD, k=None):
    if k is None:
        k = P["k_base"]
    if VM <= 0:
        return np.nan
        
    # Calculate the analytical peak of the Marginal Benefit curve.
    # We only look for stable equilibria to the right of this peak.
    e_peak = np.exp(P["a_bar"] - P["sigma"]**2)
    
    # If the FOC is already negative at the peak, the prize is too small to justify effort
    if foc(e_peak, VM, VD, k) < 0:
        return np.nan
        
    try:
        # Bracket from the peak to a sufficiently large number
        return brentq(foc, e_peak, 1000.0, args=(VM, VD, k), xtol=1e-8)
    except Exception:
        return np.nan

# ─────────────────────────────────────────────
# STYLE
# ─────────────────────────────────────────────
STYLE = dict(fig_bg="#F8FAFC", ax_bg="#FFFFFF", grid="#E2E8F0", spine="#CBD5E1")

def style_ax(ax, title, xlabel, ylabel):
    ax.set_facecolor(STYLE["ax_bg"])
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(True, color=STYLE["grid"], linewidth=0.8, linestyle="--", alpha=0.7)
    for sp in ax.spines.values():
        sp.set_color(STYLE["spine"])
        sp.set_linewidth(0.8)

# ─────────────────────────────────────────────
# DATA GENERATION
# ─────────────────────────────────────────────
M_grid     = np.linspace(10, 400, 300)
pbar_grid  = np.linspace(3.5, 20, 200)
delta      = 0.7

# Panel 1: e* vs M
data_p1 = {}
for label, tau_fn in TAU_FORMS.items():
    data_p1[label] = np.array([solve_effort(V_mono(M, tau_fn), V_duo(M)) for M in M_grid])

# Panel 2: e* vs tau'(M) (FIXED AXIS LOGIC)
# Sweep the target derivative directly so all curves align on the x-axis
tp_grid = np.linspace(0.001, 0.08, 200)
data_p2 = {}

for label, tau_fn in TAU_FORMS.items():
    e_arr = []
    for tp in tp_grid:
        # Back out the required alpha to hit the target tau'(M) at M_base
        if "Linear" in label:
            P["alpha"] = tp
        elif "Log" in label:
            P["alpha"] = tp * P["M_base"]
        elif "Power" in label:
            P["alpha"] = tp / (P["beta"] * P["M_base"]**(P["beta"] - 1))
            
        VM = V_mono(P["M_base"], tau_fn)
        VD = V_duo(P["M_base"])
        e = solve_effort(VM, VD)
        e_arr.append(e)
        
    data_p2[label] = (tp_grid, np.array(e_arr))

P["alpha"] = 0.02 # Reset alpha to baseline

# Panel 3: e* vs rent gap — sweep p_bar at M_base
data_p3 = {}
for label, tau_fn in TAU_FORMS.items():
    gap_arr, e_arr = [], []
    for pb in pbar_grid:
        P["p_bar"] = pb
        VM  = V_mono(P["M_base"], tau_fn)
        VD  = V_duo(P["M_base"])
        gap = VM - VD
        e   = solve_effort(VM, VD)
        gap_arr.append(gap); e_arr.append(e)
    data_p3[label] = (np.array(gap_arr), np.array(e_arr))
P["p_bar"] = 10.0 # Reset p_bar to baseline

# Panel 4: Period 2 effort by prior state (Power form)
tau_power = TAU_FORMS[r"Power ($\tau = \alpha M^{\beta}$)"]
e2_mono = np.array([solve_effort(V_mono(M, tau_power), V_duo(M), k=P["k_base"] * delta**2) for M in M_grid])
e2_duo  = np.array([solve_effort(V_mono(M, tau_power), V_duo(M), k=P["k_base"] * delta)    for M in M_grid])
e2_fail = np.array([solve_effort(V_mono(M, tau_power), V_duo(M), k=P["k_base"])             for M in M_grid])

# ─────────────────────────────────────────────
# FIGURE PLOTTING
# ─────────────────────────────────────────────
fig = plt.figure(figsize=(14, 10), facecolor=STYLE["fig_bg"])
fig.suptitle(
    "Comparative Statics: ACA Expansion and Pharmaceutical Innovation",
    fontsize=13, fontweight="bold", y=0.98
)

gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32,
                       left=0.07, right=0.97, top=0.93, bottom=0.07)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

lw = 2.2

# Panel 1
style_ax(ax1, "Panel 1 — $e^*$ vs Market Size $M$",
         "Market Size $M$", "Equilibrium Effort $e^*$")
for label in TAU_FORMS:
    ax1.plot(M_grid, data_p1[label], lw=lw, color=COLORS[label], ls=DASH[label], label=label)
ax1.axvline(P["M_base"], color="#94A3B8", lw=1, ls=":", label="$M_0 = 100$")
ax1.legend(fontsize=7.5, framealpha=0.9)

# Panel 2
style_ax(ax2, r"Panel 2 — Rebate Sensitivity: $e^*$ vs $\tau'(M)$",
         r"Rebate Sensitivity $\tau'(M)$", "Equilibrium Effort $e^*$")
for label in TAU_FORMS:
    tp_grid_plot, e_arr = data_p2[label]
    ax2.plot(tp_grid_plot, e_arr, lw=lw, color=COLORS[label], ls=DASH[label], label=label)
ax2.legend(fontsize=7.5, framealpha=0.9)

# Panel 3
style_ax(ax3, r"Panel 3 — Rent Gap: $e^*$ vs $(V^M - V^D)$",
         r"Rent Gap $V^M - V^D$", "Equilibrium Effort $e^*$")
for label in TAU_FORMS:
    gap_arr, e_arr = data_p3[label]
    ax3.plot(gap_arr, e_arr, lw=lw, color=COLORS[label], ls=DASH[label], label=label)
ax3.legend(fontsize=7.5, framealpha=0.9)

# Panel 4
style_ax(ax4, r"Panel 4 — Dynamic Ranking: $e^*_2$ by Prior State" + "\n(Power $\\tau$ form)",
         "Market Size $M$", r"Period 2 Effort $e^*_2$")
c_mono, c_duo, c_fail = "#7C3AED", "#0891B2", "#DC2626"
ax4.plot(M_grid, e2_mono, lw=lw, color=c_mono, ls="-",  label="Prior: Monopoly win")
ax4.plot(M_grid, e2_duo,  lw=lw, color=c_duo,  ls="--", label="Prior: Duopoly win")
ax4.plot(M_grid, e2_fail, lw=lw, color=c_fail, ls="-.", label="Prior: Failure")
ax4.fill_between(M_grid, e2_mono, e2_duo,  alpha=0.07, color=c_mono)
ax4.fill_between(M_grid, e2_duo,  e2_fail, alpha=0.07, color=c_duo)
ax4.legend(fontsize=7.5, framealpha=0.9)

# Footer
fig.text(
    0.5, 0.01,
    r"Baseline: $\sigma=1.0$,  $\bar{a}=0.5$,  $\bar{p}=10$,  $\bar{c}=2$,  "
    r"$k=1.5$,  $\alpha=0.02$,  $\beta=0.75$,  $c \sim U[0.5,\,3.5]$",
    ha="center", fontsize=8, color="#64748B"
)

# ─────────────────────────────────────────────
# SAVE & SHOW
# ─────────────────────────────────────────────
# Define the output directory and create it if it doesn't exist
output_dir = "theoretical_model"
os.makedirs(output_dir, exist_ok=True)

# Define the full file path
output_path = os.path.join(output_dir, "HT_week8_model.png")

plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor=STYLE["fig_bg"])
print(f"Saved: {output_path}")
plt.show()