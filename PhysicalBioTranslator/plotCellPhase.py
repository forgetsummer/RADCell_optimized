#!/usr/bin/env python3
"""
================================================================================
Cell Phase Analysis Script
================================================================================
Purpose:
    Analyze and visualize cell phase distribution data from the 
    CellStateModel simulation.

Input:
    - cellPhase.csv: CSV file with columns (Step, G0, G1, S, G2, M, TotalCells)

Output:
    - cellPhase_analysis.png: PNG plot of phase distribution and population
    - cellPhase_analysis.svg: SVG plot (vector format)
    - Console statistics summary

Usage:
    python3 plotCellPhase.py [path_to_cellPhase.csv]
    
    If no path is provided, it looks for ./phase_output/cellPhase.csv

Author: Generated for PhysicalBioTranslator module
Date: December 2024
================================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def analyze_cell_phase(fileName):
    """
    Analyze cell phase data and generate plots.
    
    Args:
        fileName: Path to the cellPhase.csv file
    """
    # Load the data
    data = np.genfromtxt(fileName, delimiter=',', skip_header=1)

    print("=" * 60)
    print("Cell Phase Analysis")
    print("=" * 60)
    print(f"Input file: {fileName}")
    print(f"Data shape: {data.shape}")
    print("Columns: Step, G0, G1, S, G2, M, TotalCells")

    # Extract columns
    time = data[:,0]  # Step number
    G0Num = data[:,1]
    G1Num = data[:,2]
    SNum = data[:,3]
    G2Num = data[:,4]
    MNum = data[:,5]
    num_total = data[:,6]

    # Calculate ratios
    G0Ratio = G0Num / num_total
    G1Ratio = G1Num / num_total
    SRatio = SNum / num_total
    G2Ratio = G2Num / num_total
    MRatio = MNum / num_total

    # Print statistics
    print("\n--- Statistics ---")
    print(f"Initial cells: {int(num_total[0])}")
    print(f"Final cells: {int(num_total[-1])}")
    print(f"Population growth: {(num_total[-1]/num_total[0] - 1)*100:.1f}%")
    print(f"Total time steps: {len(time)}")

    print("\nInitial phase distribution:")
    print(f"  G0: {G0Ratio[0]*100:.1f}%")
    print(f"  G1: {G1Ratio[0]*100:.1f}%")
    print(f"  S:  {SRatio[0]*100:.1f}%")
    print(f"  G2: {G2Ratio[0]*100:.1f}%")
    print(f"  M:  {MRatio[0]*100:.1f}%")

    print("\nFinal phase distribution:")
    print(f"  G0: {G0Ratio[-1]*100:.1f}%")
    print(f"  G1: {G1Ratio[-1]*100:.1f}%")
    print(f"  S:  {SRatio[-1]*100:.1f}%")
    print(f"  G2: {G2Ratio[-1]*100:.1f}%")
    print(f"  M:  {MRatio[-1]*100:.1f}%")

    # Convert time steps to hours (assuming T=60 seconds per step)
    time_hours = time * 60 / 3600  # 60 seconds per step, convert to hours

    # Save plots in the same directory as the input file
    output_dir = os.path.dirname(fileName)
    if not output_dir:
        output_dir = '.'

    # Figure 1: Phase distribution over time
    plt.figure(figsize=(10, 6))
    plt.plot(time_hours, G0Ratio, 'm-', label='G0 Phase (Quiescent)', linewidth=1)
    plt.plot(time_hours, G1Ratio, 'k-', label='G1 Phase', linewidth=1)
    plt.plot(time_hours, SRatio, 'b-', label='S Phase (DNA Synthesis)', linewidth=1)
    plt.plot(time_hours, G2Ratio, 'g-', label='G2 Phase', linewidth=1)
    plt.plot(time_hours, MRatio, 'r-', label='M Phase (Mitosis)', linewidth=1)

    plt.xlabel('Time (hours)', fontsize=12)
    plt.ylabel('Cell Phase Ratio', fontsize=12)
    plt.title('Cell Phase Distribution Over Time', fontsize=14)
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.ylim(0, 1)
    plt.tight_layout()
    
    phase_png = os.path.join(output_dir, 'cellPhase_distribution.png')
    phase_svg = os.path.join(output_dir, 'cellPhase_distribution.svg')
    plt.savefig(phase_png, dpi=150)
    plt.savefig(phase_svg)

    # Figure 2: Total cell count over time
    plt.figure(figsize=(10, 6))
    plt.plot(time_hours, num_total, 'b-', linewidth=2)
    plt.xlabel('Time (hours)', fontsize=12)
    plt.ylabel('Total Cell Count', fontsize=12)
    plt.title('Cell Population Growth', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    growth_png = os.path.join(output_dir, 'cellPopulation_growth.png')
    growth_svg = os.path.join(output_dir, 'cellPopulation_growth.svg')
    plt.savefig(growth_png, dpi=150)
    plt.savefig(growth_svg)
    
    print("\n--- Output Files ---")
    print(f"  Phase distribution: {phase_png}, {phase_svg}")
    print(f"  Population growth:  {growth_png}, {growth_svg}")
    print("=" * 60)

    plt.show()


if __name__ == "__main__":
    # Default file path
    default_path = "./phase_output/cellPhase.csv"
    
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
    elif os.path.exists(default_path):
        file_path = default_path
    else:
        print(f"Error: No file specified and default file '{default_path}' not found.")
        print("Usage: python3 plotCellPhase.py [path_to_cellPhase.csv]")
        sys.exit(1)
    
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
    
    analyze_cell_phase(file_path)

