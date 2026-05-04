#!/usr/bin/env python3
"""
Plot results from multi-component joint calibration test.
Compares LQ target, analytical calibrated model, and MC simulation SF curves.
"""

import os
import sys
import csv
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "./multicomponent_calibration_output"
OUT_DIR = os.path.join(DATA_DIR, "plots")
os.makedirs(OUT_DIR, exist_ok=True)

summary_file = os.path.join(DATA_DIR, "mc_calibration_summary.csv")
rows = []
with open(summary_file) as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append(r)

cell_lines = [r["CellLine"] for r in rows]
anal_errs = [float(r["AvgLog10Error_Analytical"]) for r in rows]
sim_errs = [float(r["AvgLog10Error_Sim"]) for r in rows]

def sanitize(name):
    out = ""
    for c in name:
        if c in r'/\:*?"<>|':
            out += '_'
        else:
            out += c
    return out

def load_cell_csv(cell_line):
    fname = os.path.join(DATA_DIR, f"mc_calib_{sanitize(cell_line)}.csv")
    doses, sf_lq, sf_cal, sf_sim, sf_sim_sd = [], [], [], [], []
    with open(fname) as f:
        reader = csv.DictReader(f)
        for r in reader:
            doses.append(float(r["Dose"]))
            sf_lq.append(float(r["SF_LQ"]))
            sf_cal.append(float(r["SF_Calibrated"]))
            sf_sim.append(float(r["SF_Simulated"]))
            sf_sim_sd.append(float(r["SF_Simulated_SD"]))
    return np.array(doses), np.array(sf_lq), np.array(sf_cal), np.array(sf_sim), np.array(sf_sim_sd)

# ---- Figure 1: SF curves for all cell lines ----
n_cells = len(rows)
n_cols = min(4, n_cells)
n_rows_grid = math.ceil(n_cells / n_cols)

fig, axes = plt.subplots(n_rows_grid, n_cols, figsize=(4.5*n_cols, 4*n_rows_grid))
if n_rows_grid == 1:
    axes = np.array([axes])
axes = axes.flatten()

for idx, r in enumerate(rows):
    ax = axes[idx]
    cl = r["CellLine"]
    alpha = float(r["Alpha_LQ"])
    beta = float(r["Beta_LQ"])
    Nc = float(r["Nc"])
    wp = float(r["omega_p"])
    lp = float(r["lambda_p"])
    sim_err = float(r["AvgLog10Error_Sim"])
    anal_err = float(r["AvgLog10Error_Analytical"])

    doses, sf_lq, sf_cal, sf_sim, sf_sim_sd = load_cell_csv(cl)

    ax.semilogy(doses, sf_lq, 'k-', linewidth=2, label='LQ target')
    ax.semilogy(doses, sf_cal, 'b--', linewidth=1.5, label='Analytical fit')

    valid = sf_sim > 0
    if np.any(valid):
        ax.errorbar(doses[valid], sf_sim[valid], yerr=sf_sim_sd[valid],
                     fmt='ro', markersize=5, capsize=3, label='MC simulation')

    ax.set_xlabel("Dose (Gy)")
    ax.set_ylabel("Surviving Fraction")
    ax.set_ylim(bottom=1e-4, top=1.5)
    ax.set_title(f"{cl}\n(a={alpha:.3f}, b={beta:.3f})\n"
                 f"Nc={Nc:.0f}, wp={wp:.1f}, lp={lp:.3f}\n"
                 f"Anal={anal_err:.3f}, Sim={sim_err:.3f}",
                 fontsize=9)
    ax.legend(fontsize=7, loc='lower left')
    ax.grid(True, alpha=0.3)

for idx in range(n_cells, len(axes)):
    axes[idx].set_visible(False)

fig.suptitle("Multi-Component Calibration: SF Curves (1000 cells, 5 reps)",
             fontsize=14, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "sf_curves_all.png"), dpi=150, bbox_inches='tight')
print(f"Saved: {os.path.join(OUT_DIR, 'sf_curves_all.png')}")

# ---- Figure 2: Error comparison bar chart ----
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

x = np.arange(len(cell_lines))
width = 0.35

ax1.bar(x - width/2, anal_errs, width, label='Analytical', color='steelblue')
ax1.bar(x + width/2, sim_errs, width, label='Simulation', color='salmon')
ax1.set_xlabel("Cell Line")
ax1.set_ylabel("Avg |dlog10(SF)|")
ax1.set_title("Error: Analytical vs Simulation")
ax1.set_xticks(x)
ax1.set_xticklabels(cell_lines, rotation=45, ha='right', fontsize=9)
ax1.legend()
ax1.grid(True, axis='y', alpha=0.3)

ratios = [s/a if a > 0 else 0 for s, a in zip(sim_errs, anal_errs)]
ax2.bar(x, ratios, color='mediumpurple')
ax2.set_xlabel("Cell Line")
ax2.set_ylabel("Sim Error / Analytical Error")
ax2.set_title("Simulation-to-Analytical Error Ratio")
ax2.set_xticks(x)
ax2.set_xticklabels(cell_lines, rotation=45, ha='right', fontsize=9)
ax2.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
ax2.grid(True, axis='y', alpha=0.3)

overall_anal = np.mean(anal_errs)
overall_sim = np.mean(sim_errs)
fig2.suptitle(f"Multi-Component Calibration Errors\n"
              f"Overall: Analytical={overall_anal:.3f}, Simulation={overall_sim:.3f}",
              fontsize=13, fontweight='bold')
fig2.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "error_comparison.png"), dpi=150, bbox_inches='tight')
print(f"Saved: {os.path.join(OUT_DIR, 'error_comparison.png')}")

# ---- Figure 3: Calibrated parameter distribution ----
fig3, axes3 = plt.subplots(2, 3, figsize=(15, 8))

params_to_plot = [
    ("Nc", [float(r["Nc"]) for r in rows]),
    ("omega_p", [float(r["omega_p"]) for r in rows]),
    ("lambda_p", [float(r["lambda_p"]) for r in rows]),
    ("e2", [float(r["e2"]) for r in rows]),
    ("e3", [float(r["e3"]) for r in rows]),
    ("a", [float(r["a"]) for r in rows]),
]

for i, (pname, pvals) in enumerate(params_to_plot):
    ax = axes3.flatten()[i]
    ax.bar(x, pvals, color='teal', alpha=0.7)
    ax.set_xlabel("Cell Line")
    ax.set_ylabel(pname)
    ax.set_title(f"Calibrated {pname}")
    ax.set_xticks(x)
    ax.set_xticklabels(cell_lines, rotation=45, ha='right', fontsize=8)
    ax.grid(True, axis='y', alpha=0.3)

fig3.suptitle("Multi-Component Calibration: Parameter Distribution",
              fontsize=14, fontweight='bold')
fig3.tight_layout()
fig3.savefig(os.path.join(OUT_DIR, "parameter_distribution.png"), dpi=150, bbox_inches='tight')
print(f"Saved: {os.path.join(OUT_DIR, 'parameter_distribution.png')}")

# ---- Figure 4: Residuals (log10 space) for selected cell lines ----
best_3 = sorted(range(len(sim_errs)), key=lambda i: sim_errs[i])[:3]
worst_3 = sorted(range(len(sim_errs)), key=lambda i: sim_errs[i])[-3:]
selection = best_3 + worst_3

fig4, axes4 = plt.subplots(2, 3, figsize=(15, 8))
axes4 = axes4.flatten()

for plot_idx, cell_idx in enumerate(selection):
    ax = axes4[plot_idx]
    cl = rows[cell_idx]["CellLine"]
    doses, sf_lq, sf_cal, sf_sim, sf_sim_sd = load_cell_csv(cl)

    eps = 1e-12
    resid_anal = np.log10(sf_cal + eps) - np.log10(sf_lq + eps)
    resid_sim = np.where(sf_sim > 0,
                         np.log10(sf_sim + eps) - np.log10(sf_lq + eps),
                         np.nan)

    ax.plot(doses, resid_anal, 'b-o', markersize=5, label='Analytical')
    valid = ~np.isnan(resid_sim)
    if np.any(valid):
        ax.plot(doses[valid], resid_sim[valid], 'r-s', markersize=5, label='Simulation')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel("Dose (Gy)")
    ax.set_ylabel("log10(SF_pred) - log10(SF_LQ)")
    label = "BEST" if plot_idx < 3 else "WORST"
    ax.set_title(f"{cl} ({label}, err={sim_errs[cell_idx]:.3f})", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

fig4.suptitle("Residuals: Best 3 vs Worst 3 Cell Lines",
              fontsize=14, fontweight='bold')
fig4.tight_layout()
fig4.savefig(os.path.join(OUT_DIR, "residuals_best_worst.png"), dpi=150, bbox_inches='tight')
print(f"Saved: {os.path.join(OUT_DIR, 'residuals_best_worst.png')}")

print("\nAll plots generated successfully.")
print(f"\nSummary:")
print(f"  Analytical avg error: {overall_anal:.3f}")
print(f"  Simulation avg error: {overall_sim:.3f}")
print(f"  Ratio (sim/anal):     {overall_sim/overall_anal:.1f}x")
print(f"\n  Best simulation cell line:  {cell_lines[best_3[0]]} (err={sim_errs[best_3[0]]:.3f})")
print(f"  Worst simulation cell line: {cell_lines[worst_3[-1]]} (err={sim_errs[worst_3[-1]]:.3f})")
