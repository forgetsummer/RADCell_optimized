#!/usr/bin/env python3
"""
Analyze multi-component energy model grid scan results.

Reads the combined grid_scan_results.csv and per-run summary CSVs,
generates comparison plots and identifies the best (Nc, omega_p, lambda_p).

Usage:
    python3 analyze_multicomponent_results.py <grid_scan_output_dir>
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import glob

def load_grid_results(results_dir):
    """Load the combined grid scan CSV."""
    csv_path = os.path.join(results_dir, 'grid_scan_results.csv')
    if not os.path.exists(csv_path):
        print(f"ERROR: {csv_path} not found")
        sys.exit(1)

    data = []
    with open(csv_path) as f:
        header = f.readline().strip().split(',')
        for line in f:
            vals = line.strip().split(',')
            if len(vals) < 8:
                continue
            row = {
                'Nc': float(vals[0]),
                'omega_p': float(vals[1]),
                'lambda_p': float(vals[2]),
                'CellLine': vals[3],
                'Alpha_LQ': float(vals[4]),
                'Beta_LQ': float(vals[5]),
                'AvgLog10Error': float(vals[6]),
                'MaxLog10Error': float(vals[7]),
            }
            data.append(row)
    return data

def load_per_run_csvs(results_dir):
    """Load per-cell-line CSVs from each subdirectory for SF curves."""
    runs = {}
    for subdir in sorted(glob.glob(os.path.join(results_dir, '*/'))):
        dirname = os.path.basename(os.path.normpath(subdir))
        csvs = glob.glob(os.path.join(subdir, 'mc_validation_*.csv'))
        if not csvs:
            continue
        run_data = {}
        for csv_path in csvs:
            fname = os.path.basename(csv_path)
            if fname == 'mc_validation_summary.csv':
                continue
            cell_name = fname.replace('mc_validation_', '').replace('.csv', '')
            try:
                arr = np.genfromtxt(csv_path, delimiter=',', skip_header=1)
                if arr.ndim == 1:
                    arr = arr.reshape(1, -1)
                run_data[cell_name] = arr
            except Exception:
                continue
        if run_data:
            runs[dirname] = run_data
    return runs

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 analyze_multicomponent_results.py <grid_scan_output_dir>")
        sys.exit(1)

    results_dir = sys.argv[1]
    output_dir = os.path.join(results_dir, 'analysis')
    os.makedirs(output_dir, exist_ok=True)

    print("Loading grid scan results...")
    data = load_grid_results(results_dir)
    if not data:
        print("No data found.")
        sys.exit(1)

    # Aggregate: for each (Nc, omega_p, lambda_p), compute mean error across cell lines
    combos = {}
    for row in data:
        key = (row['Nc'], row['omega_p'], row['lambda_p'])
        if key not in combos:
            combos[key] = []
        combos[key].append(row['AvgLog10Error'])

    combo_means = {}
    for key, errors in combos.items():
        combo_means[key] = np.mean(errors)

    # Sort and find best
    sorted_combos = sorted(combo_means.items(), key=lambda x: x[1])
    best_key, best_error = sorted_combos[0]

    print(f"\n{'='*60}")
    print(f"  GRID SCAN RESULTS")
    print(f"{'='*60}")
    print(f"  Total combinations: {len(combo_means)}")
    print(f"  Best: Nc={best_key[0]}, omega_p={best_key[1]}, lambda_p={best_key[2]}")
    print(f"  Best mean |Δlog10(SF)|: {best_error:.4f}")

    # Check if baseline exists
    baseline_key = None
    for key in combo_means:
        if key[0] >= 9999:
            baseline_key = key
            break
    if baseline_key:
        baseline_error = combo_means[baseline_key]
        improvement = baseline_error - best_error
        print(f"\n  Baseline mean |Δlog10(SF)|: {baseline_error:.4f}")
        print(f"  Improvement: {improvement:.4f} ({improvement/baseline_error*100:.1f}%)")

    # Print top 10
    print(f"\n  Top 10 parameter combinations:")
    print(f"  {'Nc':>6} {'omega_p':>8} {'lambda_p':>9} {'Mean Error':>11}")
    print(f"  {'-'*38}")
    for (nc, op, lp), err in sorted_combos[:10]:
        marker = " <-- BEST" if (nc, op, lp) == best_key else ""
        print(f"  {nc:6.0f} {op:8.2f} {lp:9.4f} {err:11.4f}{marker}")

    # Save results table
    table_path = os.path.join(output_dir, 'grid_scan_ranking.csv')
    with open(table_path, 'w') as f:
        f.write("Rank,Nc,omega_p,lambda_p,MeanAvgLog10Error\n")
        for rank, ((nc, op, lp), err) in enumerate(sorted_combos, 1):
            f.write(f"{rank},{nc},{op},{lp},{err:.6f}\n")
    print(f"\n  Ranking saved to: {table_path}")

    # --- Plot 1: Heatmap of error vs omega_p and Nc (averaged over lambda_p) ---
    nc_vals = sorted(set(k[0] for k in combo_means if k[0] < 9999))
    op_vals = sorted(set(k[1] for k in combo_means))
    lp_vals = sorted(set(k[2] for k in combo_means if k[0] < 9999))

    if len(nc_vals) >= 2 and len(op_vals) >= 2:
        fig, axes = plt.subplots(1, len(lp_vals), figsize=(5*len(lp_vals), 4), squeeze=False)
        for idx, lp in enumerate(lp_vals):
            ax = axes[0][idx]
            grid = np.full((len(nc_vals), len(op_vals)), np.nan)
            for i, nc in enumerate(nc_vals):
                for j, op in enumerate(op_vals):
                    key = (nc, op, lp)
                    if key in combo_means:
                        grid[i, j] = combo_means[key]

            im = ax.imshow(grid, aspect='auto', origin='lower',
                          extent=[op_vals[0]-0.25, op_vals[-1]+0.25,
                                  nc_vals[0]-7.5, nc_vals[-1]+7.5])
            ax.set_xlabel('omega_p')
            ax.set_ylabel('Nc')
            ax.set_title(f'lambda_p = {lp}')
            ax.set_xticks(op_vals)
            ax.set_yticks(nc_vals)
            plt.colorbar(im, ax=ax, label='Mean |Δlog10(SF)|')

            for i, nc in enumerate(nc_vals):
                for j, op in enumerate(op_vals):
                    if not np.isnan(grid[i, j]):
                        ax.text(op, nc, f'{grid[i,j]:.3f}', ha='center', va='center',
                               fontsize=8, color='white' if grid[i,j] > np.nanmean(grid) else 'black')

        plt.suptitle('Multi-Component Energy: Error Heatmap', fontsize=14)
        plt.tight_layout()
        heatmap_path = os.path.join(output_dir, 'error_heatmap.png')
        plt.savefig(heatmap_path, dpi=150, bbox_inches='tight')
        plt.savefig(heatmap_path.replace('.png', '.pdf'), bbox_inches='tight')
        plt.close()
        print(f"  Heatmap saved to: {heatmap_path}")

    # --- Plot 2: Per-cell-line error comparison (baseline vs best) ---
    cell_lines = sorted(set(row['CellLine'] for row in data))

    baseline_errors = {}
    best_errors = {}
    for row in data:
        key = (row['Nc'], row['omega_p'], row['lambda_p'])
        if baseline_key and key == baseline_key:
            baseline_errors[row['CellLine']] = row['AvgLog10Error']
        if key == best_key:
            best_errors[row['CellLine']] = row['AvgLog10Error']

    if baseline_errors and best_errors:
        common_cells = sorted(set(baseline_errors.keys()) & set(best_errors.keys()))
        if common_cells:
            fig, ax = plt.subplots(figsize=(10, 5))
            x = np.arange(len(common_cells))
            w = 0.35
            bl_vals = [baseline_errors[c] for c in common_cells]
            be_vals = [best_errors[c] for c in common_cells]
            ax.bar(x - w/2, bl_vals, w, label='Baseline', color='#4C72B0', alpha=0.8)
            ax.bar(x + w/2, be_vals, w,
                   label=f'Best (Nc={best_key[0]}, ωp={best_key[1]}, λp={best_key[2]})',
                   color='#DD8452', alpha=0.8)
            ax.set_ylabel('Avg |Δlog10(SF)|')
            ax.set_title('Per-Cell-Line Error: Baseline vs Best Multi-Component')
            ax.set_xticks(x)
            ax.set_xticklabels(common_cells, rotation=45, ha='right')
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            bar_path = os.path.join(output_dir, 'per_cellline_comparison.png')
            plt.savefig(bar_path, dpi=150, bbox_inches='tight')
            plt.savefig(bar_path.replace('.png', '.pdf'), bbox_inches='tight')
            plt.close()
            print(f"  Bar chart saved to: {bar_path}")

    # --- Plot 3: SF curves for baseline vs best for each cell line ---
    runs = load_per_run_csvs(results_dir)

    baseline_dirname = 'baseline'
    best_nc, best_op, best_lp = best_key
    best_dirname = f'nc{int(best_nc)}_op{best_op}_lp{best_lp}'

    if baseline_dirname in runs and best_dirname in runs:
        bl_run = runs[baseline_dirname]
        be_run = runs[best_dirname]
        common_cells_sf = sorted(set(bl_run.keys()) & set(be_run.keys()))

        n_cells = len(common_cells_sf)
        if n_cells > 0:
            ncols = min(3, n_cells)
            nrows = (n_cells + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows), squeeze=False)

            for idx, cell_name in enumerate(common_cells_sf):
                r, c = divmod(idx, ncols)
                ax = axes[r][c]

                bl_data = bl_run[cell_name]
                be_data = be_run[cell_name]

                doses = bl_data[:, 0]
                sf_lq = bl_data[:, 1]
                sf_bl_sim = bl_data[:, 3]
                sf_be_sim = be_data[:, 3]

                mask = sf_lq > 0
                ax.semilogy(doses[mask], sf_lq[mask], 'ro-', linewidth=2, markersize=6,
                           label='LQ (exp)', zorder=3)
                mask_bl = sf_bl_sim > 0
                ax.semilogy(doses[mask_bl], sf_bl_sim[mask_bl], 'b^--', linewidth=1.5,
                           markersize=5, label='Baseline sim', alpha=0.8)
                mask_be = sf_be_sim > 0
                ax.semilogy(doses[mask_be], sf_be_sim[mask_be], 'gs-', linewidth=1.5,
                           markersize=5, label='Er+Ep sim', alpha=0.8)

                ax.set_xlabel('Dose (Gy)')
                ax.set_ylabel('Survival Fraction')
                ax.set_title(cell_name.replace('_', '/'))
                ax.legend(fontsize=7)
                ax.grid(True, alpha=0.3)
                ax.set_xlim(-0.2, 6.5)

            for idx in range(n_cells, nrows * ncols):
                r, c = divmod(idx, ncols)
                axes[r][c].set_visible(False)

            plt.suptitle(f'SF Curves: Baseline vs Best (Nc={best_nc}, ωp={best_op}, λp={best_lp})',
                        fontsize=14)
            plt.tight_layout()
            sf_path = os.path.join(output_dir, 'sf_curves_comparison.png')
            plt.savefig(sf_path, dpi=150, bbox_inches='tight')
            plt.savefig(sf_path.replace('.png', '.pdf'), bbox_inches='tight')
            plt.close()
            print(f"  SF curves saved to: {sf_path}")

    print(f"\nAll analysis output in: {output_dir}")

if __name__ == '__main__':
    main()
