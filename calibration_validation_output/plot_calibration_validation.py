#!/usr/bin/env python3
"""
Plot calibration validation results:
1. Survival fraction curves (LQ, Calibrated, Simulated)
2. Cell state distribution over time for different doses
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import glob
import re

# Set style
plt.rcParams['figure.figsize'] = (10, 7)
plt.rcParams['font.size'] = 12

def plot_survival_curves(comparison_file='comparison.csv', output_dir='.', 
                        num_replicates=5, alpha_lq=0.75, beta_lq=0.11):
    """
    Plot survival fraction curves comparing LQ, Calibrated, and Simulated results.
    Format matches the reference figure style.
    
    Args:
        comparison_file: Path to comparison.csv file
        output_dir: Directory to save plots
        num_replicates: Number of Monte Carlo replicates
        alpha_lq: LQ model alpha parameter
        beta_lq: LQ model beta parameter
    """
    print(f"Reading survival data from: {comparison_file}")
    
    # Load data
    data = np.genfromtxt(comparison_file, delimiter=',', skip_header=1)
    
    dose = data[:, 0]
    sf_lq = data[:, 1]
    sf_calib = data[:, 2]
    sf_sim = data[:, 3]
    sf_sim_sd = data[:, 4]
    abs_log10_error = data[:, 5]
    
    # Average log-space error across all doses (includes high-dose tail)
    avg_log10_error = float(np.mean(abs_log10_error)) if len(abs_log10_error) > 0 else 0.0
    avg_factor_error = 10.0 ** avg_log10_error
    
    # Create figure with larger size for better readability
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Replace zero values with small value for log scale (but keep original for markers)
    sf_sim_plot = np.where(sf_sim > 0, sf_sim, 1e-5)  # Use small value for zero
    
    # Plot curves with exact format from reference figure:
    # 1. LQ Model: Red, solid line, circles
    ax.semilogy(dose, sf_lq, 'r-', linewidth=2.5, marker='o', 
                markersize=8, label=f'LQ Model (α={alpha_lq}, β={beta_lq}) - Target', 
                alpha=0.9, zorder=3, markeredgecolor='darkred', markeredgewidth=1)
    
    # 2. Calibrated Model: Green, dashed line, squares
    ax.semilogy(dose, sf_calib, 'g--', linewidth=2.0, marker='s', 
                markersize=7, label='Calibrated Model (analytical)', 
                alpha=0.9, zorder=2, markeredgecolor='darkgreen', markeredgewidth=1)
    
    # 3. Simulated: Blue, dotted line, triangles (plot only non-zero for line)
    mask_sim_nonzero = sf_sim > 0
    if np.any(mask_sim_nonzero):
        # Plot line only for non-zero points
        ax.semilogy(dose[mask_sim_nonzero], sf_sim[mask_sim_nonzero], 'b:', 
                    linewidth=2.0, label=f'Simulated SF ({num_replicates} replicates)', 
                    alpha=0.9, zorder=1)
        # Plot markers for non-zero points
        ax.semilogy(dose[mask_sim_nonzero], sf_sim[mask_sim_nonzero], 'b^', 
                    markersize=7, alpha=0.9, zorder=1, 
                    markeredgecolor='darkblue', markeredgewidth=1, linestyle='None')
    
    # Add error bars for simulated data (if SD > 0)
    mask_sd = (sf_sim_sd > 0) & (sf_sim > 0)
    if np.any(mask_sd):
        ax.errorbar(dose[mask_sd], sf_sim[mask_sd], yerr=sf_sim_sd[mask_sd], 
                   fmt='none', color='b', alpha=0.4, capsize=3, zorder=0)
    
    # Formatting
    ax.set_xlabel('Dose (Gy)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Survival Fraction', fontsize=14, fontweight='bold')
    
    # Title and subtitle
    ax.set_title('Survival Fraction Comparison\n(LQ Model vs Calibrated Model vs Simulation)', 
                fontsize=16, fontweight='bold', pad=10)
    
    # Add subtitle for model type
    ax.text(0.5, 0.98, 'Repair Mediated Model With Misrepair', 
           transform=ax.transAxes, fontsize=12, style='italic',
           ha='center', va='top')
    
    # Legend in upper right
    legend = ax.legend(loc='upper right', fontsize=11, framealpha=0.95, 
                      edgecolor='black', fancybox=False)
    legend.get_frame().set_linewidth(1.2)
    
    # Grid
    ax.grid(True, alpha=0.3, which='both', linestyle='-', linewidth=0.5)
    ax.set_xlim(-0.2, max(dose) * 1.05)
    ax.set_ylim(1e-4, 1.2)
    
    # Add info box (similar to reference figure)
    info_text = (f'Avg |Δlog10(SF)|: {avg_log10_error:.3f} (×{avg_factor_error:.2f})\n'
                 f'Monte Carlo: {num_replicates} replicates/dose')
    ax.text(0.02, 0.02, info_text, transform=ax.transAxes, 
           fontsize=10, verticalalignment='bottom',
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7, edgecolor='black'),
           zorder=10)
    
    plt.tight_layout()
    
    # Save plots
    png_file = os.path.join(output_dir, 'survival_fraction_comparison.png')
    pdf_file = os.path.join(output_dir, 'survival_fraction_comparison.pdf')
    svg_file = os.path.join(output_dir, 'survival_fraction_comparison.svg')
    
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.savefig(svg_file, bbox_inches='tight')
    
    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")
    print(f"Saved: {svg_file}")
    
    return fig


def plot_cell_state_evolution(state_files_pattern='state_evolution_dose_*.csv', 
                               output_dir='.', doses_to_plot=None):
    """
    Plot cell state distribution over time for different doses.
    
    Args:
        state_files_pattern: Glob pattern for state evolution CSV files
        output_dir: Directory to save plots
        doses_to_plot: List of doses to plot (None = plot all)
    """
    # Find all state evolution files
    state_files = sorted(glob.glob(state_files_pattern))
    
    if not state_files:
        print(f"Warning: No files found matching pattern: {state_files_pattern}")
        return
    
    print(f"Found {len(state_files)} state evolution files")
    
    # Extract doses from filenames
    doses = []
    for f in state_files:
        # Extract dose from filename like "state_evolution_dose_2.0Gy.csv"
        try:
            dose_str = f.split('dose_')[1].split('Gy')[0]
            doses.append(float(dose_str))
        except:
            print(f"Warning: Could not extract dose from {f}")
            continue
    
    # Filter doses if specified
    if doses_to_plot is not None:
        filtered_files = []
        filtered_doses = []
        for f, d in zip(state_files, doses):
            if d in doses_to_plot:
                filtered_files.append(f)
                filtered_doses.append(d)
        state_files = filtered_files
        doses = filtered_doses
    
    # Create figure with subplots
    n_plots = len(state_files)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_plots == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Color scheme
    colors = {'S1': 'blue', 'S2': 'green', 'S3': 'red'}
    
    for idx, (state_file, dose) in enumerate(zip(state_files, doses)):
        ax = axes[idx]
        
        # Load data
        data = np.genfromtxt(state_file, delimiter=',', skip_header=1)
        
        time_hours = data[:, 0]
        s1_fraction = data[:, 4]
        s2_fraction = data[:, 5]
        s3_fraction = data[:, 6]
        
        # Plot state fractions
        ax.plot(time_hours, s1_fraction, color=colors['S1'], linewidth=2, 
                label='S1 (Healthy)', alpha=0.8)
        ax.plot(time_hours, s2_fraction, color=colors['S2'], linewidth=2, 
                label='S2 (Repair)', alpha=0.8)
        ax.plot(time_hours, s3_fraction, color=colors['S3'], linewidth=2, 
                label='S3 (Death)', alpha=0.8)
        
        ax.set_xlabel('Time (hours)', fontsize=11)
        ax.set_ylabel('Cell State Fraction', fontsize=11)
        ax.set_title(f'Dose = {dose:.1f} Gy', fontsize=12)
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.05, 1.05)
    
    # Hide unused subplots
    for idx in range(n_plots, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Cell State Evolution Over Time', fontsize=16, y=1.02)
    plt.tight_layout()
    
    # Save plots
    png_file = os.path.join(output_dir, 'cell_state_evolution.png')
    pdf_file = os.path.join(output_dir, 'cell_state_evolution.pdf')
    svg_file = os.path.join(output_dir, 'cell_state_evolution.svg')
    
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.savefig(svg_file, bbox_inches='tight')
    
    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")
    print(f"Saved: {svg_file}")
    
    return fig


def plot_cell_state_comparison(state_files_pattern='state_evolution_dose_*.csv',
                               output_dir='.', selected_doses=[0, 1, 2, 4]):
    """
    Plot cell state fractions at final time point vs dose.
    
    Args:
        state_files_pattern: Glob pattern for state evolution CSV files
        output_dir: Directory to save plots
        selected_doses: List of doses to include
    """
    # Find all state evolution files
    state_files = sorted(glob.glob(state_files_pattern))
    
    if not state_files:
        print(f"Warning: No files found matching pattern: {state_files_pattern}")
        return
    
    # Extract doses and final state fractions
    doses = []
    final_s1 = []
    final_s2 = []
    final_s3 = []
    
    for f in state_files:
        try:
            dose_str = f.split('dose_')[1].split('Gy')[0]
            dose = float(dose_str)
            
            if selected_doses is None or dose in selected_doses:
                data = np.genfromtxt(f, delimiter=',', skip_header=1)
                doses.append(dose)
                final_s1.append(data[-1, 4])  # S1 fraction
                final_s2.append(data[-1, 5])  # S2 fraction
                final_s3.append(data[-1, 6])  # S3 fraction
        except Exception as e:
            print(f"Warning: Error processing {f}: {e}")
            continue
    
    if not doses:
        print("No valid data found")
        return
    
    # Sort by dose
    sorted_indices = np.argsort(doses)
    doses = np.array(doses)[sorted_indices]
    final_s1 = np.array(final_s1)[sorted_indices]
    final_s2 = np.array(final_s2)[sorted_indices]
    final_s3 = np.array(final_s3)[sorted_indices]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot stacked area or bars
    width = 0.6
    x_pos = np.arange(len(doses))
    
    ax.bar(x_pos, final_s1, width, label='S1 (Healthy)', color='blue', alpha=0.7)
    ax.bar(x_pos, final_s2, width, bottom=final_s1, label='S2 (Repair)', 
           color='green', alpha=0.7)
    ax.bar(x_pos, final_s3, width, bottom=final_s1+final_s2, label='S3 (Death)', 
           color='red', alpha=0.7)
    
    ax.set_xlabel('Dose (Gy)', fontsize=14)
    ax.set_ylabel('Final Cell State Fraction', fontsize=14)
    ax.set_title('Final Cell State Distribution vs Dose', fontsize=16)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f'{d:.1f}' for d in doses])
    ax.legend(loc='best', fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1)
    
    plt.tight_layout()
    
    # Save plots
    png_file = os.path.join(output_dir, 'cell_state_vs_dose.png')
    pdf_file = os.path.join(output_dir, 'cell_state_vs_dose.pdf')
    svg_file = os.path.join(output_dir, 'cell_state_vs_dose.svg')
    
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.savefig(svg_file, bbox_inches='tight')
    
    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")
    print(f"Saved: {svg_file}")
    
    return fig


if __name__ == "__main__":
    # Determine output directory (same as script location)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if not script_dir:
        script_dir = '.'
    
    print("=" * 60)
    print("Calibration Validation Plotting")
    print("=" * 60)
    
    # Read parameters from calibrated_params.txt if available
    alpha_lq = 0.75  # Default values
    beta_lq = 0.11
    num_replicates = 5  # Default (from our test)
    
    params_file = os.path.join(script_dir, 'calibrated_params.txt')
    if os.path.exists(params_file):
        with open(params_file, 'r') as f:
            content = f.read()
            # Extract alpha_LQ
            alpha_match = re.search(r'alpha_LQ\s*=\s*([\d.]+)', content)
            if alpha_match:
                alpha_lq = float(alpha_match.group(1))
            # Extract beta_LQ
            beta_match = re.search(r'beta_LQ\s*=\s*([\d.]+)', content)
            if beta_match:
                beta_lq = float(beta_match.group(1))
    
    # Plot 1: Survival fraction curves
    comparison_file = os.path.join(script_dir, 'comparison.csv')
    if os.path.exists(comparison_file):
        print("\n1. Plotting survival fraction curves...")
        plot_survival_curves(comparison_file, script_dir, 
                           num_replicates=num_replicates,
                           alpha_lq=alpha_lq, beta_lq=beta_lq)
        plt.close('all')
    else:
        print(f"\nWarning: {comparison_file} not found")
    
    # Plot 2: Cell state evolution over time
    state_pattern = os.path.join(script_dir, 'state_evolution_dose_*.csv')
    if glob.glob(state_pattern):
        print("\n2. Plotting cell state evolution over time...")
        plot_cell_state_evolution(state_pattern, script_dir, 
                                 doses_to_plot=[0, 0.5, 1, 1.5, 2, 2.5, 3, 4])
        plt.close('all')
    else:
        print(f"\nWarning: No state evolution files found matching: {state_pattern}")
    
    # Plot 3: Final cell state vs dose
    if glob.glob(state_pattern):
        print("\n3. Plotting final cell state distribution vs dose...")
        plot_cell_state_comparison(state_pattern, script_dir, 
                                  selected_doses=None)  # Plot all doses
        plt.close('all')
    
    print("\n" + "=" * 60)
    print("Plotting complete!")
    print("=" * 60)

