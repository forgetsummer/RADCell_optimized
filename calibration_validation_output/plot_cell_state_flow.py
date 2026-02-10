#!/usr/bin/env python3
"""
Generate complete cell state transition flow diagram based on calibrated parameters
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
import numpy as np
import os
import re

def read_calibrated_params(params_file='calibrated_params.txt'):
    """Read calibrated parameters from file"""
    params = {
        'E1': 0.0,
        'E2': 36.23,
        'E3': 69.49,
        'sigma': 10.0,
        'alpha': 0.529,
        'kappa': 40.0,
        'T21': 10.0,
        'T23': 10.0,
        'k_error': 1e-10,
        'alpha_LQ': 0.75,
        'beta_LQ': 0.11
    }
    
    if os.path.exists(params_file):
        with open(params_file, 'r') as f:
            content = f.read()
            
            # Extract values using regex
            patterns = {
                'E1': r'E1\s*=\s*([\d.]+)',
                'E2': r'E2\s*=\s*([\d.]+)',
                'E3': r'E3\s*=\s*([\d.]+)',
                'sigma': r'sigma\s*=\s*([\d.]+)',
                'alpha': r'alpha\s*=\s*([\d.]+)',
                'kappa': r'kappa\s*=\s*([\d.]+)',
                'T21': r'T21\s*=\s*([\d.]+)',
                'T23': r'T23\s*=\s*([\d.]+)',
                'k_error': r'k_error\s*=\s*([\d.e-]+)',
                'alpha_LQ': r'alpha_LQ\s*=\s*([\d.]+)',
                'beta_LQ': r'beta_LQ\s*=\s*([\d.]+)'
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, content)
                if match:
                    try:
                        params[key] = float(match.group(1))
                    except:
                        pass
    
    return params

def calculate_dsb_thresholds(params):
    """Calculate DSB thresholds from energy thresholds"""
    # E = alpha * N_DSB, so N_DSB = E / alpha
    dsb_thresholds = {
        'E1_DSB': params['E1'] / params['alpha'] if params['alpha'] > 0 else 0,
        'E2_DSB': params['E2'] / params['alpha'] if params['alpha'] > 0 else 0,
        'E3_DSB': params['E3'] / params['alpha'] if params['alpha'] > 0 else 0
    }
    
    # Dose thresholds (using kappa = DSB/Gy)
    dose_thresholds = {
        'E2_dose': params['E2'] / (params['alpha'] * params['kappa']) if params['alpha'] * params['kappa'] > 0 else 0,
        'E3_dose': params['E3'] / (params['alpha'] * params['kappa']) if params['alpha'] * params['kappa'] > 0 else 0
    }
    
    return dsb_thresholds, dose_thresholds

def plot_cell_state_flow(output_dir='.', params_file='calibrated_params.txt'):
    """Generate complete cell state transition flow diagram"""
    
    # Read parameters
    params = read_calibrated_params(params_file)
    dsb_thresh, dose_thresh = calculate_dsb_thresholds(params)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # Define state positions (centered layout)
    state_positions = {
        'S1': (2, 7),      # Top left - Healthy
        'S2': (5, 7),      # Top center - Repair
        'S3': (8, 7),      # Top right - Death
        'DSB': (2, 4),     # Bottom left - DSB input
        'Repair': (5, 4),  # Bottom center - Repair process
        'Checkpoint': (8, 4)  # Bottom right - Checkpoint
    }
    
    # Define colors
    colors = {
        'S1': '#4CAF50',      # Green - Healthy
        'S2': '#FF9800',      # Orange - Repair
        'S3': '#F44336',      # Red - Death
        'DSB': '#2196F3',     # Blue - Damage
        'Repair': '#9C27B0',  # Purple - Repair
        'Checkpoint': '#607D8B'  # Gray - Checkpoint
    }
    
    # Draw state boxes
    state_boxes = {}
    box_width = 1.5
    box_height = 1.0
    
    for state, pos in state_positions.items():
        if state.startswith('S'):
            # Main state boxes
            x, y = pos
            box = FancyBboxPatch((x - box_width/2, y - box_height/2), 
                               box_width, box_height,
                               boxstyle="round,pad=0.1", 
                               facecolor=colors[state],
                               edgecolor='black', linewidth=2.5,
                               alpha=0.8)
            ax.add_patch(box)
            state_boxes[state] = box
            
            # State labels
            state_labels = {
                'S1': 'S1\nHealthy\nState',
                'S2': 'S2\nRepair\nState',
                'S3': 'S3\nDeath\nState'
            }
            ax.text(x, y, state_labels[state], 
                   ha='center', va='center', fontsize=14, fontweight='bold',
                   color='white')
    
    # Draw parameter boxes
    param_y = 1.5
    
    # E1, E2, E3 thresholds
    ax.text(2, param_y, f"E1 = {params['E1']:.1f}\n(Healthy baseline)", 
           ha='center', va='center', fontsize=10,
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    ax.text(5, param_y, f"E2 = {params['E2']:.1f}\n(Repair threshold)\n{dose_thresh['E2_dose']:.2f} Gy\n{dsb_thresh['E2_DSB']:.1f} DSB", 
           ha='center', va='center', fontsize=10,
           bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))
    
    ax.text(8, param_y, f"E3 = {params['E3']:.1f}\n(Death threshold)\n{dose_thresh['E3_dose']:.2f} Gy\n{dsb_thresh['E3_DSB']:.1f} DSB", 
           ha='center', va='center', fontsize=10,
           bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))
    
    # Draw transition arrows with labels
    arrow_style = dict(arrowstyle='->', lw=2.5, color='black')
    
    # S1 -> S2 (damage-induced transition)
    arrow1 = FancyArrowPatch((state_positions['S1'][0] + box_width/2, state_positions['S1'][1]),
                            (state_positions['S2'][0] - box_width/2, state_positions['S2'][1]),
                            connectionstyle="arc3,rad=0", **arrow_style)
    ax.add_patch(arrow1)
    ax.text(3.5, 7.3, 'Damage\nE > E2', ha='center', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # S1 -> S3 (direct lethal)
    arrow2 = FancyArrowPatch((state_positions['S1'][0] + box_width/2, state_positions['S1'][1] - 0.2),
                            (state_positions['S3'][0] - box_width/2, state_positions['S3'][1] - 0.2),
                            connectionstyle="arc3,rad=-0.3", **arrow_style)
    ax.add_patch(arrow2)
    ax.text(5, 6.2, 'Direct Lethal\nE > E3', ha='center', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # S2 -> S1 (repair)
    arrow3 = FancyArrowPatch((state_positions['S2'][0] - box_width/2, state_positions['S2'][1] - 0.3),
                            (state_positions['S1'][0] + box_width/2, state_positions['S1'][1] - 0.3),
                            connectionstyle="arc3,rad=0.3", **arrow_style)
    ax.add_patch(arrow3)
    ax.text(3.5, 6.5, f'Repair\nT21 = {params["T21"]:.1f} h', ha='center', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    # S2 -> S3 (death progression)
    arrow4 = FancyArrowPatch((state_positions['S2'][0] + box_width/2, state_positions['S2'][1] - 0.3),
                            (state_positions['S3'][0] - box_width/2, state_positions['S3'][1] - 0.3),
                            connectionstyle="arc3,rad=-0.3", **arrow_style)
    ax.add_patch(arrow4)
    ax.text(6.5, 6.5, f'Death\nT23 = {params["T23"]:.1f} h', ha='center', fontsize=9,
           bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
    
    # S3 is absorbing (self-loop)
    circle = Circle((state_positions['S3'][0], state_positions['S3'][1] - 0.6), 0.15,
                   fill=False, edgecolor='black', linewidth=2)
    ax.add_patch(circle)
    ax.text(state_positions['S3'][0], state_positions['S3'][1] - 0.9, 'Absorbing\nState', 
           ha='center', fontsize=9, style='italic')
    
    # Add DSB input box
    dsb_box = FancyBboxPatch((state_positions['DSB'][0] - 0.7, state_positions['DSB'][1] - 0.4),
                           1.4, 0.8,
                           boxstyle="round,pad=0.05",
                           facecolor=colors['DSB'], edgecolor='black', linewidth=2,
                           alpha=0.7)
    ax.add_patch(dsb_box)
    ax.text(state_positions['DSB'][0], state_positions['DSB'][1], 
           f'DSB Input\nκ = {params["kappa"]:.0f} DSB/Gy\nα = {params["alpha"]:.3f}\nE = α × N_DSB',
           ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    # Arrow from DSB to S1
    arrow_dsb = FancyArrowPatch((state_positions['DSB'][0], state_positions['DSB'][1] + 0.4),
                               (state_positions['S1'][0], state_positions['S1'][1] - box_height/2),
                               connectionstyle="arc3,rad=0.2", **arrow_style)
    ax.add_patch(arrow_dsb)
    
    # Add misrepair box
    misrepair_box = FancyBboxPatch((state_positions['Repair'][0] - 0.7, state_positions['Repair'][1] - 0.4),
                                 1.4, 0.8,
                                 boxstyle="round,pad=0.05",
                                 facecolor=colors['Repair'], edgecolor='black', linewidth=2,
                                 alpha=0.7)
    ax.add_patch(misrepair_box)
    k_error_str = f'{params["k_error"]:.2e}' if params['k_error'] > 1e-10 else 'Disabled'
    ax.text(state_positions['Repair'][0], state_positions['Repair'][1],
           f'Misrepair\nk_error = {k_error_str}\nper (DSB×hour)',
           ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    # Arrow from misrepair to S2->S3 transition
    arrow_misrepair = FancyArrowPatch((state_positions['Repair'][0] + 0.7, state_positions['Repair'][1]),
                                     (6.2, 6.7),
                                     connectionstyle="arc3,rad=-0.1", 
                                     arrowstyle='->', lw=1.5, color='purple', linestyle='--')
    ax.add_patch(arrow_misrepair)
    ax.text(6.0, 5.8, 'Misrepair\n→ Death', ha='center', fontsize=8, style='italic',
           color='purple')
    
    # Add checkpoint box
    checkpoint_box = FancyBboxPatch((state_positions['Checkpoint'][0] - 0.7, state_positions['Checkpoint'][1] - 0.4),
                                  1.4, 0.8,
                                  boxstyle="round,pad=0.05",
                                  facecolor=colors['Checkpoint'], edgecolor='black', linewidth=2,
                                  alpha=0.7)
    ax.add_patch(checkpoint_box)
    ax.text(state_positions['Checkpoint'][0], state_positions['Checkpoint'][1],
           'Checkpoint\nG2→M + S1→S2\n(Original Design)',
           ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    # Arrow from checkpoint to S2
    arrow_checkpoint = FancyArrowPatch((state_positions['Checkpoint'][0], state_positions['Checkpoint'][1] + 0.4),
                                     (state_positions['S2'][0], state_positions['S2'][1] - box_height/2),
                                     connectionstyle="arc3,rad=-0.2",
                                     arrowstyle='->', lw=1.5, color='gray', linestyle='--')
    ax.add_patch(arrow_checkpoint)
    
    # Add model parameters box (top center)
    param_box = FancyBboxPatch((4, 9.2), 2, 0.6,
                              boxstyle="round,pad=0.1",
                              facecolor='lightblue', edgecolor='black', linewidth=2,
                              alpha=0.9)
    ax.add_patch(param_box)
    ax.text(5, 9.5, 'Repair Mediated Model With Misrepair', 
           ha='center', va='center', fontsize=12, fontweight='bold')
    
    # Add transition probability formulas (right side)
    formula_x = 9.2
    formula_y_start = 8.5
    
    formulas = [
        'Transition Probabilities:',
        f'P(S1→S2) = 2Φ(-|E - E2|/2σ)',
        f'P(S1→S3) = 2Φ(-|E - E3|/2σ)',
        f'P(S2→S1) = 2Φ(-|E - E1|/2σ)',
        f'P(S2→S3) = 2Φ(-|E - E3|/2σ)',
        '',
        'Timescales:',
        f'T21 = {params["T21"]:.1f} h (repair)',
        f'T23 = {params["T23"]:.1f} h (death)',
        '',
        f'σ = {params["sigma"]:.1f}'
    ]
    
    for i, formula in enumerate(formulas):
        ax.text(formula_x, formula_y_start - i*0.25, formula,
               ha='left', va='center', fontsize=8,
               family='monospace')
    
    # Add LQ model reference (bottom)
    lq_text = f'LQ Model: α={params["alpha_LQ"]:.2f} Gy⁻¹, β={params["beta_LQ"]:.2f} Gy⁻²'
    ax.text(5, 0.3, lq_text, ha='center', va='center', fontsize=10,
           bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plots
    png_file = os.path.join(output_dir, 'cell_state_transition_flow.png')
    pdf_file = os.path.join(output_dir, 'cell_state_transition_flow.pdf')
    svg_file = os.path.join(output_dir, 'cell_state_transition_flow.svg')
    
    plt.savefig(png_file, dpi=150, bbox_inches='tight')
    plt.savefig(pdf_file, bbox_inches='tight')
    plt.savefig(svg_file, bbox_inches='tight')
    
    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")
    print(f"Saved: {svg_file}")
    
    return fig

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if not script_dir:
        script_dir = '.'
    
    print("=" * 60)
    print("Generating Cell State Transition Flow Diagram")
    print("=" * 60)
    
    params_file = os.path.join(script_dir, 'calibrated_params.txt')
    plot_cell_state_flow(script_dir, params_file)
    
    print("=" * 60)
    print("Flow diagram generation complete!")
    print("=" * 60)


