"""
Comparative Statics: ACA Expansion & Pharmaceutical Innovation
================================================================
Model: TT Week 1 — Mohammad-Shah Karim

Class-indexed extension of the HT8 model. Rebate elasticity beta_c
varies across therapeutic classes through beta_c = beta(n_c), where
n_c is the number of clinically interchangeable molecules in class c.

Six panels:
  1. e* vs M (three tau forms) — same as HT8
  2. e* vs tau'(M)              — same as HT8
  3. e* vs rent gap V^M - V^D   — same as HT8
  4. Period-2 effort by prior   — same as HT8
  5. NEW: e* vs M_c for several beta_c values (low/mid/high n_c)
  6. NEW: Class-specific threshold M_c* as a function of beta_c

Run with: python theoretical_model/TT_week1_model.py
Saves: theoretical_model/TT_week1_model.png
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import norm
from scipy.optimize import brentq

# ─────────────────────────────────────────────
# PARAMETERS (baseline class)
# ─────────────────────────────────────────────
P = dict(
    sigma   = 1.0,
    a_bar   = 0.5,
    p_bar   = 10.0,
    c_bar   = 2.0,
    k_base  = 1.0,
    M_base  = 100.0,
    alpha   = 0.02,   # alpha_c at baseline
    beta    = 0.75,   # beta_c at baseline
    c_lo    = 0.5,
    c_hi    = 3.5,
)

# ─────────────────────────────────────────────
# TAU FORMS
# ─────────────────────────────────────────────
TAU_FORMS = {
    r"Power ($\tau_c = \alpha_c M_c^{\beta_c}$)":  lambda M: P["alpha"] * M ** P["beta"],
    r"Log ($\tau_c = \alpha_c \ln M_c$)":          lambda M: P["alpha"] * np.log(np.maximum(M, 1e-6)),
    r"Linear ($\tau_c = \alpha_c M_c$)":           lambda M: P["alpha"] * M,
}

COLORS = {
    r"Power ($\tau_c = \alpha_c M_c^{\beta_c}$)":  "#2563EB",
    r"Log ($\tau_c = \alpha_c \ln M_c$)":          "#DC2626",
    r"Linear ($\tau_c = \alpha_c M_c$)":           "#16A34A",
}

DASH = {
    r"Power ($\tau_c = \alpha_c M_c^{\beta_c}$)":  "-",
    r"Log ($\tau_c = \alpha_c \ln M_c$)":          "--",
    r"Linear ($\tau_c = \alpha_c M_c$)":           "-.",
}

# ─────────────────────────────────────────────
# MARKET VALUES (now class-indexed in interpretation)
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
# EQUILIBRIUM SOLVER
# ─────────────────────────────────────────────
def foc(e, VM, VD, k):
    return Q_prime(e) * ((1 - Q(e)) * VM + Q(e) * VD) - k * e

def solve_effort(VM, VD, k=None):
    if k is None:
        k = P["k_base"]
    if VM <= 0:
        return np.nan
    e_peak = np.exp(P["a_bar"] - P["sigma"]**2)
    if foc(e_peak, VM, VD, k) < 0:
        return np.nan
    try:
        return brentq(foc, e_peak, 1000.0, args=(VM, VD, k), xtol=1e-8)
    except Exception:
        return np.nan

# ─────────────────────────────────────────────
# CLASS-SPECIFIC THRESHOLD (closed form, power tau)
# ─────────────────────────────────────────────
def M_star(alpha_c, beta_c, p_bar=None, c_bar=None):
    """Class-specific threshold M_c* under power tau: dV^M/dM = 0."""
    if p_bar is None:
        p_bar = P["p_bar"]
    if c_bar is None:
        c_bar = P["c_bar"]
    rent = p_bar - c_bar
    if rent <= 0:
        return np.nan
    return (rent / (alpha_c * (1.0 + beta_c))) ** (1.0 / beta_c)

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

# --- Panel 1: e* vs M (three tau forms) ---
data_p1 = {}
for label, tau_fn in TAU_FORMS.items():
    data_p1[label] = np.array([solve_effort(V_mono(M, tau_fn), V_duo(M)) for M in M_grid])

# --- Panel 2: e* vs tau'(M) ---
tp_grid = np.linspace(0.001, 0.08, 200)
data_p2 = {}
for label, tau_fn in TAU_FORMS.items():
    e_arr = []
    for tp in tp_grid:
        if "Linear" in label:
            P["alpha"] = tp
        elif "Log" in label:
            P["alpha"] = tp * P["M_base"]
        elif "Power" in label:
            P["alpha"] = tp / (P["beta"] * P["M_base"]**(P["beta"] - 1))
        VM = V_mono(P["M_base"], tau_fn)
        VD = V_duo(P["M_base"])
        e_arr.append(solve_effort(VM, VD))
    data_p2[label] = (tp_grid, np.array(e_arr))
P["alpha"] = 0.02  # reset

# --- Panel 3: e* vs rent gap ---
data_p3 = {}
for label, tau_fn in TAU_FORMS.items():
    gap_arr, e_arr = [], []
    for pb in pbar_grid:
        P["p_bar"] = pb
        VM  = V_mono(P["M_base"], tau_fn)
        VD  = V_duo(P["M_base"])
        gap_arr.append(VM - VD)
        e_arr.append(solve_effort(VM, VD))
    data_p3[label] = (np.array(gap_arr), np.array(e_arr))
P["p_bar"] = 10.0  # reset

# --- Panel 4: Period 2 effort by prior state (power tau) ---
tau_power = TAU_FORMS[r"Power ($\tau_c = \alpha_c M_c^{\beta_c}$)"]
e2_mono = np.array([solve_effort(V_mono(M, tau_power), V_duo(M), k=P["k_base"] * delta**2) for M in M_grid])
e2_duo  = np.array([solve_effort(V_mono(M, tau_power), V_duo(M), k=P["k_base"] * delta)    for M in M_grid])
e2_fail = np.array([solve_effort(V_mono(M, tau_power), V_duo(M), k=P["k_base"])             for M in M_grid])

# --- Panel 5: e* vs M for several class types (varying beta_c via n_c) ---
# Three stylised classes: low / mid / high substitutability n_c
classes_p5 = [
    ("Low $n_c$ (e.g. H, S): $\\beta_c = 0.55$", 0.55, "#16A34A", "-"),
    ("Mid $n_c$ (typical):  $\\beta_c = 0.75$",  0.75, "#2563EB", "--"),
    ("High $n_c$ (e.g. N): $\\beta_c = 1.05$",   1.05, "#DC2626", "-."),
]
data_p5 = {}
for label, beta_c, col, ls in classes_p5:
    P["beta"] = beta_c
    e_arr = np.array([solve_effort(V_mono(M, tau_power), V_duo(M)) for M in M_grid])
    data_p5[label] = (e_arr, col, ls, beta_c)
P["beta"] = 0.75  # reset

# --- Panel 6: M_c* as a function of beta_c ---
beta_grid = np.linspace(0.30, 1.40, 200)
M_star_grid = np.array([M_star(P["alpha"], b) for b in beta_grid])

# ─────────────────────────────────────────────
# FIGURE PLOTTING
# ─────────────────────────────────────────────
fig = plt.figure(figsize=(14, 14), facecolor=STYLE["fig_bg"])
fig.suptitle(
    "Comparative Statics: Class-Indexed Pharmaceutical Innovation",
    fontsize=13, fontweight="bold", y=0.99
)

gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.32,
                       left=0.07, right=0.97, top=0.95, bottom=0.06)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])
ax5 = fig.add_subplot(gs[2, 0])
ax6 = fig.add_subplot(gs[2, 1])

lw = 2.2

# Panel 1
style_ax(ax1, "Panel 1 — $e_c^*$ vs Class Market Size $M_c$",
         "Class Market Size $M_c$", "Equilibrium Effort $e_c^*$")
for label in TAU_FORMS:
    ax1.plot(M_grid, data_p1[label], lw=lw, color=COLORS[label], ls=DASH[label], label=label)
ax1.axvline(P["M_base"], color="#94A3B8", lw=1, ls=":", label="$M_0 = 100$")
ax1.legend(fontsize=7.5, framealpha=0.9)

# Panel 2
style_ax(ax2, r"Panel 2 — Rebate Sensitivity: $e_c^*$ vs $\tau_c'(M_c)$",
         r"Rebate Sensitivity $\tau_c'(M_c)$", "Equilibrium Effort $e_c^*$")
for label in TAU_FORMS:
    tp_grid_plot, e_arr = data_p2[label]
    ax2.plot(tp_grid_plot, e_arr, lw=lw, color=COLORS[label], ls=DASH[label], label=label)
ax2.legend(fontsize=7.5, framealpha=0.9)

# Panel 3
style_ax(ax3, r"Panel 3 — Rent Gap: $e_c^*$ vs $(V_c^M - V_c^D)$",
         r"Class Rent Gap $V_c^M - V_c^D$", "Equilibrium Effort $e_c^*$")
for label in TAU_FORMS:
    gap_arr, e_arr = data_p3[label]
    ax3.plot(gap_arr, e_arr, lw=lw, color=COLORS[label], ls=DASH[label], label=label)
ax3.legend(fontsize=7.5, framealpha=0.9)

# Panel 4
style_ax(ax4, r"Panel 4 — Dynamic Ranking: $e_{c,2}^*$ by Prior State" + "\n(Power $\\tau_c$ form)",
         "Class Market Size $M_c$", r"Period 2 Effort $e_{c,2}^*$")
c_mono, c_duo, c_fail = "#7C3AED", "#0891B2", "#DC2626"
ax4.plot(M_grid, e2_mono, lw=lw, color=c_mono, ls="-",  label="Prior: Monopoly win")
ax4.plot(M_grid, e2_duo,  lw=lw, color=c_duo,  ls="--", label="Prior: Duopoly win")
ax4.plot(M_grid, e2_fail, lw=lw, color=c_fail, ls="-.", label="Prior: Failure")
ax4.fill_between(M_grid, e2_mono, e2_duo,  alpha=0.07, color=c_mono)
ax4.fill_between(M_grid, e2_duo,  e2_fail, alpha=0.07, color=c_duo)
ax4.legend(fontsize=7.5, framealpha=0.9)

# Panel 5 — NEW: e_c* vs M_c across stylised classes (varying beta_c)
style_ax(ax5,
         "Panel 5 — Class Heterogeneity: $e_c^*$ vs $M_c$ for varying $\\beta_c = \\beta(n_c)$",
         "Class Market Size $M_c$", "Equilibrium Effort $e_c^*$")
for label, (e_arr, col, ls, beta_c) in data_p5.items():
    ax5.plot(M_grid, e_arr, lw=lw, color=col, ls=ls, label=label)
    # Mark class-specific threshold M_c*
    Mc_star = M_star(P["alpha"], beta_c)
    if not np.isnan(Mc_star) and 10 <= Mc_star <= 400:
        ax5.axvline(Mc_star, color=col, lw=1.0, ls=":", alpha=0.55)
ax5.legend(fontsize=7.5, framealpha=0.9, loc="best")

# Panel 6 — NEW: M_c* as a function of beta_c
style_ax(ax6,
         "Panel 6 — Class-Specific Threshold $M_c^*$ as a Function of $\\beta_c$",
         r"Rebate Elasticity $\beta_c = \beta(n_c)$",
         r"Threshold $M_c^*$ (log scale)")
ax6.semilogy(beta_grid, M_star_grid, lw=lw, color="#7C3AED")
for label, _, col, _ in classes_p5:
    beta_c = float(label.split("=")[-1].strip().replace("$", ""))
    Mc_star = M_star(P["alpha"], beta_c)
    ax6.scatter([beta_c], [Mc_star], s=45, color=col, zorder=5,
                label=f"$\\beta_c = {beta_c:.2f}$  $\\Rightarrow$  $M_c^* \\approx {Mc_star:,.0f}$")
ax6.axhline(P["M_base"], color="#94A3B8", lw=1, ls=":", label="$M_0 = 100$")
ax6.legend(fontsize=7.5, framealpha=0.9, loc="best")

# Footer
fig.text(
    0.5, 0.005,
    r"Baseline: $\sigma=1.0$,  $\bar{a}=0.5$,  $\bar{p}=10$,  $\bar{c}=2$,  "
    r"$k=1.0$,  $\alpha_c=0.02$,  $\beta_c=0.75$ (baseline),  $c \sim U[0.5,\,3.5]$. "
    r"Panels 5–6 use the power $\tau_c$ form.",
    ha="center", fontsize=8, color="#64748B"
)

# ─────────────────────────────────────────────
# SAVE
# ─────────────────────────────────────────────
output_dir = "theoretical_model"
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "TT_week1_model.png")

plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor=STYLE["fig_bg"])
print(f"Saved: {output_path}")
plt.show()
