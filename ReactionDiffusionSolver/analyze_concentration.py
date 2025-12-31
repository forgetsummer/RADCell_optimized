#!/usr/bin/env python3
"""
============================================================================
Analysis Script: Concentration and Integral Concentration Comparison
============================================================================

Purpose:
    This script analyzes and compares the concentration and integral 
    concentration results from the DiffusionReactionSolver simulation.

Author: Generated for RADCell project analysis
        Based on work by Ruirui Liu, Oregon State University

Usage:
    python3 analyze_concentration.py [concentration_folder] [time_step]
    
    Examples:
        python3 analyze_concentration.py                    # Uses defaults
        python3 analyze_concentration.py ./concentration/   # Specify folder
        python3 analyze_concentration.py ./concentration/ 500  # Specify time step

Output:
    - Console statistics summary
    - Visualization plots saved as PNG files
============================================================================
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob

def load_concentration_file(filepath):
    """Load a concentration CSV file and return as numpy arrays."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 4:
                i, j, k, c = int(parts[0]), int(parts[1]), int(parts[2]), float(parts[3])
                data.append([i, j, k, c])
    
    if not data:
        return None, None, None, None
    
    data = np.array(data)
    return data[:, 0].astype(int), data[:, 1].astype(int), data[:, 2].astype(int), data[:, 3]

def get_grid_dimensions(i_arr, j_arr, k_arr):
    """Get the grid dimensions from index arrays."""
    return int(i_arr.max()) + 1, int(j_arr.max()) + 1, int(k_arr.max()) + 1

def reshape_to_3d(i_arr, j_arr, k_arr, values, nx, ny, nz):
    """Reshape flat data to 3D grid."""
    grid = np.zeros((nx, ny, nz))
    for idx in range(len(values)):
        grid[int(i_arr[idx]), int(j_arr[idx]), int(k_arr[idx])] = values[idx]
    return grid

def analyze_single_timestep(folder, timestep):
    """Analyze concentration and integral concentration for a single timestep."""
    
    conc_file = os.path.join(folder, f"concentration-{timestep}.csv")
    integral_file = os.path.join(folder, f"integralConcentration-{timestep}.csv")
    
    if not os.path.exists(conc_file):
        print(f"Error: Concentration file not found: {conc_file}")
        return None
    
    if not os.path.exists(integral_file):
        print(f"Error: Integral concentration file not found: {integral_file}")
        return None
    
    # Load data
    print(f"\nLoading data for timestep {timestep}...")
    i_c, j_c, k_c, conc = load_concentration_file(conc_file)
    i_ic, j_ic, k_ic, integral_conc = load_concentration_file(integral_file)
    
    if conc is None or integral_conc is None:
        print("Error: Could not load data files")
        return None
    
    # Get dimensions
    nx, ny, nz = get_grid_dimensions(i_c, j_c, k_c)
    print(f"Grid dimensions: {nx} x {ny} x {nz}")
    
    # Reshape to 3D grids
    conc_3d = reshape_to_3d(i_c, j_c, k_c, conc, nx, ny, nz)
    integral_3d = reshape_to_3d(i_ic, j_ic, k_ic, integral_conc, nx, ny, nz)
    
    return {
        'timestep': timestep,
        'nx': nx, 'ny': ny, 'nz': nz,
        'concentration': conc_3d,
        'integral_concentration': integral_3d,
        'conc_flat': conc,
        'integral_flat': integral_conc
    }

def print_statistics(data):
    """Print statistical summary of the data."""
    print("\n" + "="*60)
    print(f"  STATISTICS FOR TIMESTEP {data['timestep']}")
    print("="*60)
    
    conc = data['conc_flat']
    integral = data['integral_flat']
    
    print("\n--- Concentration ---")
    print(f"  Min:    {conc.min():.6e}")
    print(f"  Max:    {conc.max():.6e}")
    print(f"  Mean:   {conc.mean():.6e}")
    print(f"  Std:    {conc.std():.6e}")
    print(f"  Sum:    {conc.sum():.6e}")
    print(f"  Non-zero points: {np.count_nonzero(conc)} / {len(conc)}")
    
    print("\n--- Integral Concentration ---")
    print(f"  Min:    {integral.min():.6e}")
    print(f"  Max:    {integral.max():.6e}")
    print(f"  Mean:   {integral.mean():.6e}")
    print(f"  Std:    {integral.std():.6e}")
    print(f"  Sum:    {integral.sum():.6e}")
    print(f"  Non-zero points: {np.count_nonzero(integral)} / {len(integral)}")
    
    # Correlation
    if conc.std() > 0 and integral.std() > 0:
        correlation = np.corrcoef(conc, integral)[0, 1]
        print(f"\n--- Correlation ---")
        print(f"  Pearson correlation: {correlation:.6f}")

def plot_2d_slice(data, output_folder="."):
    """Create 2D slice plots for visualization."""
    
    timestep = data['timestep']
    conc = data['concentration']
    integral = data['integral_concentration']
    nz = data['nz']
    
    # Take middle slice in z direction
    z_slice = nz // 2
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Concentration
    im1 = axes[0].imshow(conc[:, :, z_slice].T, origin='lower', cmap='hot', aspect='auto')
    axes[0].set_title(f'Concentration (z={z_slice})\nTimestep {timestep}')
    axes[0].set_xlabel('X index')
    axes[0].set_ylabel('Y index')
    plt.colorbar(im1, ax=axes[0], label='Concentration')
    
    # Integral Concentration
    im2 = axes[1].imshow(integral[:, :, z_slice].T, origin='lower', cmap='viridis', aspect='auto')
    axes[1].set_title(f'Integral Concentration (z={z_slice})\nTimestep {timestep}')
    axes[1].set_xlabel('X index')
    axes[1].set_ylabel('Y index')
    plt.colorbar(im2, ax=axes[1], label='Integral Concentration')
    
    # Difference/Ratio plot
    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(conc[:, :, z_slice] > 0, 
                        integral[:, :, z_slice] / conc[:, :, z_slice], 
                        0)
    im3 = axes[2].imshow(ratio.T, origin='lower', cmap='coolwarm', aspect='auto')
    axes[2].set_title(f'Ratio (Integral/Concentration)\nTimestep {timestep}')
    axes[2].set_xlabel('X index')
    axes[2].set_ylabel('Y index')
    plt.colorbar(im3, ax=axes[2], label='Ratio')
    
    plt.tight_layout()
    output_file = os.path.join(output_folder, f'concentration_comparison_t{timestep}.png')
    plt.savefig(output_file, dpi=150)
    print(f"\nSaved 2D slice plot to: {output_file}")
    plt.close()

def plot_time_evolution(folder, output_folder="."):
    """Plot how concentration evolves over time."""
    
    # Find all concentration files
    conc_files = glob.glob(os.path.join(folder, "concentration-*.csv"))
    if not conc_files:
        print("No concentration files found for time evolution plot")
        return
    
    # Extract timesteps and sort
    timesteps = []
    for f in conc_files:
        basename = os.path.basename(f)
        ts = int(basename.replace("concentration-", "").replace(".csv", ""))
        timesteps.append(ts)
    timesteps.sort()
    
    # Sample timesteps if too many
    if len(timesteps) > 20:
        sample_indices = np.linspace(0, len(timesteps)-1, 20, dtype=int)
        timesteps = [timesteps[i] for i in sample_indices]
    
    print(f"\nAnalyzing time evolution for {len(timesteps)} timesteps...")
    
    max_conc = []
    mean_conc = []
    sum_conc = []
    max_integral = []
    mean_integral = []
    sum_integral = []
    valid_timesteps = []
    
    for ts in timesteps:
        data = analyze_single_timestep(folder, ts)
        if data is not None:
            valid_timesteps.append(ts)
            max_conc.append(data['conc_flat'].max())
            mean_conc.append(data['conc_flat'].mean())
            sum_conc.append(data['conc_flat'].sum())
            max_integral.append(data['integral_flat'].max())
            mean_integral.append(data['integral_flat'].mean())
            sum_integral.append(data['integral_flat'].sum())
    
    if not valid_timesteps:
        print("No valid data for time evolution")
        return
    
    # Create time evolution plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Concentration plots
    axes[0, 0].plot(valid_timesteps, max_conc, 'b-o', markersize=4)
    axes[0, 0].set_title('Max Concentration vs Time')
    axes[0, 0].set_xlabel('Timestep')
    axes[0, 0].set_ylabel('Max Concentration')
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].plot(valid_timesteps, mean_conc, 'b-o', markersize=4)
    axes[0, 1].set_title('Mean Concentration vs Time')
    axes[0, 1].set_xlabel('Timestep')
    axes[0, 1].set_ylabel('Mean Concentration')
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[0, 2].plot(valid_timesteps, sum_conc, 'b-o', markersize=4)
    axes[0, 2].set_title('Total Concentration vs Time')
    axes[0, 2].set_xlabel('Timestep')
    axes[0, 2].set_ylabel('Sum of Concentration')
    axes[0, 2].grid(True, alpha=0.3)
    
    # Integral concentration plots
    axes[1, 0].plot(valid_timesteps, max_integral, 'r-s', markersize=4)
    axes[1, 0].set_title('Max Integral Concentration vs Time')
    axes[1, 0].set_xlabel('Timestep')
    axes[1, 0].set_ylabel('Max Integral Concentration')
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].plot(valid_timesteps, mean_integral, 'r-s', markersize=4)
    axes[1, 1].set_title('Mean Integral Concentration vs Time')
    axes[1, 1].set_xlabel('Timestep')
    axes[1, 1].set_ylabel('Mean Integral Concentration')
    axes[1, 1].grid(True, alpha=0.3)
    
    axes[1, 2].plot(valid_timesteps, sum_integral, 'r-s', markersize=4)
    axes[1, 2].set_title('Total Integral Concentration vs Time')
    axes[1, 2].set_xlabel('Timestep')
    axes[1, 2].set_ylabel('Sum of Integral Concentration')
    axes[1, 2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_folder, 'concentration_time_evolution.png')
    plt.savefig(output_file, dpi=150)
    print(f"\nSaved time evolution plot to: {output_file}")
    plt.close()
    
    # Comparison plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    axes[0].plot(valid_timesteps, max_conc, 'b-o', label='Concentration', markersize=4)
    axes[0].plot(valid_timesteps, max_integral, 'r-s', label='Integral Concentration', markersize=4)
    axes[0].set_title('Max Values Comparison')
    axes[0].set_xlabel('Timestep')
    axes[0].set_ylabel('Value')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(valid_timesteps, sum_conc, 'b-o', label='Concentration', markersize=4)
    axes[1].plot(valid_timesteps, sum_integral, 'r-s', label='Integral Concentration', markersize=4)
    axes[1].set_title('Total Values Comparison')
    axes[1].set_xlabel('Timestep')
    axes[1].set_ylabel('Sum')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_folder, 'concentration_comparison_evolution.png')
    plt.savefig(output_file, dpi=150)
    print(f"Saved comparison evolution plot to: {output_file}")
    plt.close()

def plot_histogram(data, output_folder="."):
    """Create histogram comparison of concentration values."""
    
    timestep = data['timestep']
    conc = data['conc_flat']
    integral = data['integral_flat']
    
    # Filter out zeros for better visualization
    conc_nonzero = conc[conc > 0]
    integral_nonzero = integral[integral > 0]
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    if len(conc_nonzero) > 0:
        axes[0].hist(conc_nonzero, bins=50, color='blue', alpha=0.7, edgecolor='black')
        axes[0].set_title(f'Concentration Distribution (non-zero)\nTimestep {timestep}')
        axes[0].set_xlabel('Concentration Value')
        axes[0].set_ylabel('Frequency')
        axes[0].set_yscale('log')
    else:
        axes[0].text(0.5, 0.5, 'No non-zero values', ha='center', va='center', transform=axes[0].transAxes)
        axes[0].set_title(f'Concentration Distribution\nTimestep {timestep}')
    
    if len(integral_nonzero) > 0:
        axes[1].hist(integral_nonzero, bins=50, color='red', alpha=0.7, edgecolor='black')
        axes[1].set_title(f'Integral Concentration Distribution (non-zero)\nTimestep {timestep}')
        axes[1].set_xlabel('Integral Concentration Value')
        axes[1].set_ylabel('Frequency')
        axes[1].set_yscale('log')
    else:
        axes[1].text(0.5, 0.5, 'No non-zero values', ha='center', va='center', transform=axes[1].transAxes)
        axes[1].set_title(f'Integral Concentration Distribution\nTimestep {timestep}')
    
    plt.tight_layout()
    output_file = os.path.join(output_folder, f'concentration_histogram_t{timestep}.png')
    plt.savefig(output_file, dpi=150)
    print(f"Saved histogram plot to: {output_file}")
    plt.close()

def main():
    """Main analysis function."""
    
    print("="*60)
    print("  Concentration Analysis Tool")
    print("  RADCell Diffusion-Reaction Simulation")
    print("="*60)
    
    # Parse arguments
    folder = "./concentration/"
    timestep = None
    
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    if len(sys.argv) > 2:
        timestep = int(sys.argv[2])
    
    # Check folder exists
    if not os.path.exists(folder):
        print(f"Error: Folder not found: {folder}")
        sys.exit(1)
    
    print(f"\nAnalyzing data in: {folder}")
    
    # Find available timesteps
    conc_files = glob.glob(os.path.join(folder, "concentration-*.csv"))
    if not conc_files:
        print("Error: No concentration files found!")
        sys.exit(1)
    
    available_timesteps = []
    for f in conc_files:
        basename = os.path.basename(f)
        ts = int(basename.replace("concentration-", "").replace(".csv", ""))
        available_timesteps.append(ts)
    available_timesteps.sort()
    
    print(f"Found {len(available_timesteps)} timesteps: {available_timesteps[0]} to {available_timesteps[-1]}")
    
    # If no timestep specified, find the last one that has both files
    if timestep is None:
        # Find timesteps that have both concentration and integralConcentration files
        for ts in reversed(available_timesteps):
            integral_file = os.path.join(folder, f"integralConcentration-{ts}.csv")
            if os.path.exists(integral_file):
                timestep = ts
                break
        if timestep is None:
            print("Warning: No integralConcentration files found. Using last concentration file.")
            timestep = available_timesteps[-1]
    
    if timestep not in available_timesteps:
        print(f"Error: Timestep {timestep} not found. Available: {available_timesteps[:5]}...")
        sys.exit(1)
    
    # Analyze single timestep
    data = analyze_single_timestep(folder, timestep)
    if data is None:
        sys.exit(1)
    
    # Print statistics
    print_statistics(data)
    
    # Create output folder for plots
    output_folder = os.path.dirname(folder.rstrip('/'))
    if not output_folder:
        output_folder = "."
    
    # Generate plots
    print("\nGenerating visualization plots...")
    plot_2d_slice(data, output_folder)
    plot_histogram(data, output_folder)
    
    # Time evolution analysis
    print("\nGenerating time evolution analysis...")
    plot_time_evolution(folder, output_folder)
    
    print("\n" + "="*60)
    print("  Analysis Complete!")
    print("="*60)

if __name__ == "__main__":
    main()

