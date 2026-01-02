#!/usr/bin/env python3
"""
============================================================================
Cell Phase Transition Analysis Script
============================================================================

Purpose:
    This script analyzes and compares cell phase transition simulation results
    between different versions of the RADCell code. It also compares simulation
    results with theoretical phase ratios based on mean phase durations.

Theoretical Phase Ratios:
    For an asynchronous exponentially growing cell population, the fraction of
    cells in each phase is proportional to the phase duration:
    
    f_G1 = T_G1 / T_total
    f_S  = T_S  / T_total
    f_G2 = T_G2 / T_total
    f_M  = T_M  / T_total
    
    where T_total = T_G1 + T_S + T_G2 + T_M

Author: Ruirui Liu
        Department of Nuclear Engineering and Radiation Health Physics
        Oregon State University

Usage:
    python3 analyze_phase_comparison.py

Output:
    - phase_comparison.png/svg: Comparison plots
    - Console output with statistics

============================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# ============================================================================
# Configuration - Modify these parameters as needed
# ============================================================================

# Cell cycle parameters (from test_phaseTransition.cc)
T_G1 = 9.0      # Mean G1 duration (hours)
T_S = 11.0      # Mean S duration (hours)
T_G2 = 1.0      # Mean G2 duration (hours)
T_M = 1.0       # Mean M duration (hours)

# File paths
NEW_VERSION_FILE = './phase_output/cellPhase.csv'
OLD_VERSION_FILE = '/home/user/Geant4Projects/CellStateTransitionPaper/VersionInPhDThesis/RAD-build-new/cellPhase.csv'
OUTPUT_DIR = './phase_output/'

# ============================================================================
# Calculate Theoretical Phase Ratios
# ============================================================================

def calculate_theoretical_ratios(tG1, tS, tG2, tM):
    """
    Calculate theoretical phase ratios based on mean phase durations.
    
    For an asynchronous exponentially growing population, the fraction of cells
    in each phase is proportional to the phase duration.
    
    Args:
        tG1: Mean G1 phase duration (hours)
        tS: Mean S phase duration (hours)
        tG2: Mean G2 phase duration (hours)
        tM: Mean M phase duration (hours)
    
    Returns:
        Dictionary with theoretical ratios for each phase
    """
    T_total = tG1 + tS + tG2 + tM
    
    return {
        'G1': tG1 / T_total * 100,
        'S': tS / T_total * 100,
        'G2': tG2 / T_total * 100,
        'M': tM / T_total * 100,
        'G0': 0.0,  # No G0 in theoretical (assumes all cells cycling)
        'T_total': T_total
    }

# ============================================================================
# Load Data Functions
# ============================================================================

def load_new_version_data(filepath):
    """
    Load data from new version output file.
    Format: Step,G0,G1,S,G2,M,TotalCells (with header)
    """
    data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
    return {
        'time': data[:, 0],
        'G0': data[:, 1],
        'G1': data[:, 2],
        'S': data[:, 3],
        'G2': data[:, 4],
        'M': data[:, 5],
        'total': data[:, 6]
    }

def load_old_version_data(filepath):
    """
    Load data from old version output file.
    Format: Step,G1,S,G2,M,TotalCells (no header, no G0 column)
    """
    data = np.genfromtxt(filepath, delimiter=',')
    total = data[:, 5]
    G1 = data[:, 1]
    S = data[:, 2]
    G2 = data[:, 3]
    M = data[:, 4]
    G0 = total - (G1 + S + G2 + M)  # Calculate G0 as remainder
    
    return {
        'time': data[:, 0],
        'G0': G0,
        'G1': G1,
        'S': S,
        'G2': G2,
        'M': M,
        'total': total
    }

def convert_time_to_hours(time_steps, dt_seconds=60):
    """Convert time steps to hours."""
    return time_steps * dt_seconds / 3600

def calculate_ratios(data):
    """Calculate phase ratios as percentages."""
    total = data['total']
    return {
        'G0': data['G0'] / total * 100,
        'G1': data['G1'] / total * 100,
        'S': data['S'] / total * 100,
        'G2': data['G2'] / total * 100,
        'M': data['M'] / total * 100
    }

# ============================================================================
# Plotting Functions
# ============================================================================

def plot_phase_comparison_with_theory(new_data, old_data, theoretical, output_dir, title_suffix=""):
    """
    Create comprehensive comparison plots including theoretical ratios.
    """
    # Convert time to hours
    new_time = convert_time_to_hours(new_data['time'])
    old_time = convert_time_to_hours(old_data['time'])
    
    # Calculate ratios
    new_ratios = calculate_ratios(new_data)
    old_ratios = calculate_ratios(old_data)
    
    # Create figure with 4 subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Cell Phase Comparison: Simulation vs Theory{title_suffix}\n'
                 f'(T_G1={T_G1}h, T_S={T_S}h, T_G2={T_G2}h, T_M={T_M}h, T_total={theoretical["T_total"]}h)', 
                 fontsize=14, fontweight='bold')
    
    # Color scheme
    colors = {'old': 'blue', 'new': 'red', 'theory': 'green'}
    
    # Plot 1: Population Growth
    ax1 = axes[0, 0]
    ax1.plot(old_time, old_data['total'], 'b-', label='Old Version', linewidth=2, alpha=0.7)
    ax1.plot(new_time, new_data['total'], 'r--', label='New Version', linewidth=2, alpha=0.7)
    ax1.set_xlabel('Time (hours)', fontsize=11)
    ax1.set_ylabel('Total Cell Count', fontsize=11)
    ax1.set_title('Population Growth', fontsize=12, fontweight='bold')
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: G1 Phase with Theory
    ax2 = axes[0, 1]
    ax2.plot(old_time, old_ratios['G1'], 'b-', label='Old Version', linewidth=1.5, alpha=0.7)
    ax2.plot(new_time, new_ratios['G1'], 'r--', label='New Version', linewidth=1.5, alpha=0.7)
    ax2.axhline(y=theoretical['G1'], color='green', linestyle=':', linewidth=2, 
                label=f'Theory ({theoretical["G1"]:.1f}%)')
    ax2.set_xlabel('Time (hours)', fontsize=11)
    ax2.set_ylabel('G1 Phase (%)', fontsize=11)
    ax2.set_title('G1 Phase Ratio', fontsize=12, fontweight='bold')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: S Phase with Theory
    ax3 = axes[1, 0]
    ax3.plot(old_time, old_ratios['S'], 'b-', label='Old Version', linewidth=1.5, alpha=0.7)
    ax3.plot(new_time, new_ratios['S'], 'r--', label='New Version', linewidth=1.5, alpha=0.7)
    ax3.axhline(y=theoretical['S'], color='green', linestyle=':', linewidth=2, 
                label=f'Theory ({theoretical["S"]:.1f}%)')
    ax3.set_xlabel('Time (hours)', fontsize=11)
    ax3.set_ylabel('S Phase (%)', fontsize=11)
    ax3.set_title('S Phase Ratio', fontsize=12, fontweight='bold')
    ax3.legend(loc='best')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: G2+M Phase with Theory
    ax4 = axes[1, 1]
    old_G2M = old_ratios['G2'] + old_ratios['M']
    new_G2M = new_ratios['G2'] + new_ratios['M']
    theory_G2M = theoretical['G2'] + theoretical['M']
    ax4.plot(old_time, old_G2M, 'b-', label='Old Version', linewidth=1.5, alpha=0.7)
    ax4.plot(new_time, new_G2M, 'r--', label='New Version', linewidth=1.5, alpha=0.7)
    ax4.axhline(y=theory_G2M, color='green', linestyle=':', linewidth=2, 
                label=f'Theory ({theory_G2M:.1f}%)')
    ax4.set_xlabel('Time (hours)', fontsize=11)
    ax4.set_ylabel('G2+M Phase (%)', fontsize=11)
    ax4.set_title('G2+M Phase Ratio', fontsize=12, fontweight='bold')
    ax4.legend(loc='best')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plots
    plt.savefig(os.path.join(output_dir, 'phase_comparison_with_theory.png'), dpi=150)
    plt.savefig(os.path.join(output_dir, 'phase_comparison_with_theory.svg'))
    print(f"\nPlots saved to {output_dir}")
    
    return fig

def plot_all_phases_comparison(new_data, old_data, theoretical, output_dir, title_suffix=""):
    """
    Create a comprehensive 6-panel plot showing all phases.
    """
    # Convert time to hours
    new_time = convert_time_to_hours(new_data['time'])
    old_time = convert_time_to_hours(old_data['time'])
    
    # Calculate ratios
    new_ratios = calculate_ratios(new_data)
    old_ratios = calculate_ratios(old_data)
    
    # Create figure
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle(f'Cell Phase Comparison: Simulation vs Theory{title_suffix}\n'
                 f'(T_G1={T_G1}h, T_S={T_S}h, T_G2={T_G2}h, T_M={T_M}h)', 
                 fontsize=14, fontweight='bold')
    
    # Plot 1: Population Growth
    ax = axes[0, 0]
    ax.plot(old_time, old_data['total'], 'b-', label='Old Version', linewidth=2, alpha=0.7)
    ax.plot(new_time, new_data['total'], 'r--', label='New Version', linewidth=2, alpha=0.7)
    ax.set_xlabel('Time (hours)', fontsize=11)
    ax.set_ylabel('Total Cell Count', fontsize=11)
    ax.set_title('Population Growth', fontsize=12, fontweight='bold')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)
    
    # Plot 2: G1 Phase
    ax = axes[0, 1]
    ax.plot(old_time, old_ratios['G1'], 'b-', label='Old Version', linewidth=1.5, alpha=0.7)
    ax.plot(new_time, new_ratios['G1'], 'r--', label='New Version', linewidth=1.5, alpha=0.7)
    ax.axhline(y=theoretical['G1'], color='green', linestyle=':', linewidth=2, 
               label=f'Theory ({theoretical["G1"]:.1f}%)')
    ax.set_xlabel('Time (hours)', fontsize=11)
    ax.set_ylabel('G1 Phase (%)', fontsize=11)
    ax.set_title('G1 Phase Ratio', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Plot 3: S Phase
    ax = axes[1, 0]
    ax.plot(old_time, old_ratios['S'], 'b-', label='Old Version', linewidth=1.5, alpha=0.7)
    ax.plot(new_time, new_ratios['S'], 'r--', label='New Version', linewidth=1.5, alpha=0.7)
    ax.axhline(y=theoretical['S'], color='green', linestyle=':', linewidth=2, 
               label=f'Theory ({theoretical["S"]:.1f}%)')
    ax.set_xlabel('Time (hours)', fontsize=11)
    ax.set_ylabel('S Phase (%)', fontsize=11)
    ax.set_title('S Phase Ratio', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Plot 4: G2 Phase
    ax = axes[1, 1]
    ax.plot(old_time, old_ratios['G2'], 'b-', label='Old Version', linewidth=1.5, alpha=0.7)
    ax.plot(new_time, new_ratios['G2'], 'r--', label='New Version', linewidth=1.5, alpha=0.7)
    ax.axhline(y=theoretical['G2'], color='green', linestyle=':', linewidth=2, 
               label=f'Theory ({theoretical["G2"]:.1f}%)')
    ax.set_xlabel('Time (hours)', fontsize=11)
    ax.set_ylabel('G2 Phase (%)', fontsize=11)
    ax.set_title('G2 Phase Ratio', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Plot 5: M Phase
    ax = axes[2, 0]
    ax.plot(old_time, old_ratios['M'], 'b-', label='Old Version', linewidth=1.5, alpha=0.7)
    ax.plot(new_time, new_ratios['M'], 'r--', label='New Version', linewidth=1.5, alpha=0.7)
    ax.axhline(y=theoretical['M'], color='green', linestyle=':', linewidth=2, 
               label=f'Theory ({theoretical["M"]:.1f}%)')
    ax.set_xlabel('Time (hours)', fontsize=11)
    ax.set_ylabel('M Phase (%)', fontsize=11)
    ax.set_title('M Phase (Mitosis) Ratio', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Plot 6: G0 Phase
    ax = axes[2, 1]
    ax.plot(old_time, old_ratios['G0'], 'b-', label='Old Version', linewidth=1.5, alpha=0.7)
    ax.plot(new_time, new_ratios['G0'], 'r--', label='New Version', linewidth=1.5, alpha=0.7)
    ax.axhline(y=theoretical['G0'], color='green', linestyle=':', linewidth=2, 
               label=f'Theory ({theoretical["G0"]:.1f}%)')
    ax.set_xlabel('Time (hours)', fontsize=11)
    ax.set_ylabel('G0 Phase (%)', fontsize=11)
    ax.set_title('G0 Phase (Quiescent) Ratio', fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plots
    plt.savefig(os.path.join(output_dir, 'all_phases_comparison.png'), dpi=150)
    plt.savefig(os.path.join(output_dir, 'all_phases_comparison.svg'))
    
    return fig

def print_summary(new_data, old_data, theoretical):
    """Print summary statistics."""
    new_ratios = calculate_ratios(new_data)
    old_ratios = calculate_ratios(old_data)
    
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    
    print(f"\nTheoretical Phase Ratios (based on mean durations):")
    print(f"  T_G1 = {T_G1} hours, T_S = {T_S} hours, T_G2 = {T_G2} hours, T_M = {T_M} hours")
    print(f"  T_total = {theoretical['T_total']} hours")
    print(f"  f_G1 = {theoretical['G1']:.2f}%")
    print(f"  f_S  = {theoretical['S']:.2f}%")
    print(f"  f_G2 = {theoretical['G2']:.2f}%")
    print(f"  f_M  = {theoretical['M']:.2f}%")
    
    print("\n" + "-" * 80)
    print(f"{'Phase':<10} {'Theory':>12} {'Old (final)':>15} {'New (final)':>15} {'Old-Theory':>12} {'New-Theory':>12}")
    print("-" * 80)
    
    phases = ['G0', 'G1', 'S', 'G2', 'M']
    for phase in phases:
        theory = theoretical[phase]
        old_final = old_ratios[phase][-1]
        new_final = new_ratios[phase][-1]
        old_diff = old_final - theory
        new_diff = new_final - theory
        print(f"{phase:<10} {theory:>11.2f}% {old_final:>14.2f}% {new_final:>14.2f}% {old_diff:>+11.2f}% {new_diff:>+11.2f}%")
    
    print("\n" + "-" * 80)
    print("Final Cell Counts:")
    print(f"  Old Version: {int(old_data['total'][-1])} cells")
    print(f"  New Version: {int(new_data['total'][-1])} cells")
    
    # Calculate correlation
    min_len = min(len(old_data['total']), len(new_data['total']))
    corr = np.corrcoef(old_data['total'][:min_len], new_data['total'][:min_len])[0, 1]
    print(f"\nPopulation Growth Correlation: {corr:.6f}")
    
    # Check match
    if int(old_data['total'][-1]) == int(new_data['total'][-1]):
        print(f"\n✓ MATCH! Both versions end with {int(new_data['total'][-1])} cells")
    else:
        print(f"\n✗ MISMATCH: Old={int(old_data['total'][-1])}, New={int(new_data['total'][-1])}")

# ============================================================================
# Main Function
# ============================================================================

def main():
    """Main function to run the analysis."""
    print("=" * 80)
    print("Cell Phase Transition Analysis")
    print("=" * 80)
    
    # Calculate theoretical ratios
    theoretical = calculate_theoretical_ratios(T_G1, T_S, T_G2, T_M)
    print(f"\nTheoretical phase ratios calculated:")
    print(f"  G1: {theoretical['G1']:.2f}%, S: {theoretical['S']:.2f}%, "
          f"G2: {theoretical['G2']:.2f}%, M: {theoretical['M']:.2f}%")
    
    # Load data
    print(f"\nLoading data from:")
    print(f"  New version: {NEW_VERSION_FILE}")
    print(f"  Old version: {OLD_VERSION_FILE}")
    
    try:
        new_data = load_new_version_data(NEW_VERSION_FILE)
        print(f"  New version: {len(new_data['time'])} time points loaded")
    except Exception as e:
        print(f"  Error loading new version data: {e}")
        return
    
    try:
        old_data = load_old_version_data(OLD_VERSION_FILE)
        print(f"  Old version: {len(old_data['time'])} time points loaded")
    except Exception as e:
        print(f"  Error loading old version data: {e}")
        return
    
    # Create output directory if needed
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Generate plots
    print("\nGenerating plots...")
    plot_phase_comparison_with_theory(new_data, old_data, theoretical, OUTPUT_DIR)
    plot_all_phases_comparison(new_data, old_data, theoretical, OUTPUT_DIR)
    
    # Print summary
    print_summary(new_data, old_data, theoretical)
    
    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)
    
    plt.show()

if __name__ == "__main__":
    main()

