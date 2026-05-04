# RADCell - Radiation-induced Cell Response Simulation

A modular simulation framework for modeling radiation-induced bystander effects (RIBE) and cell state transitions in cell populations.

## Overview

RADCell integrates Monte Carlo radiation transport (via Geant4) with biological response models to simulate:
- **DNA damage** from direct radiation exposure
- **Bystander signaling** through reaction-diffusion of signaling molecules
- **Cell cycle progression** (G0, G1, S, G2, M phases)
- **Cell state transitions** (healthy → stressed → damaged states)
- **Cell proliferation and contact inhibition**

## Architecture

The project is organized into independent modules that can be built and tested separately:

```
RADCellSimulation/
├── CellMaker/                    # Cell layout initialization (Geant4-independent)
│   ├── include/
│   ├── src/
│   └── test_cellLayoutInitializer.cc
│
├── ReactionDiffusionSolver/      # Bystander signal diffusion simulation
│   ├── include/
│   ├── src/
│   └── test_reactionDiffusionKernel.cc
│
├── PhysicalBioTranslator/        # Cell state/phase models (Geant4-independent)
│   ├── include/
│   ├── src/
│   ├── test_multicomponent_validation.cc   # Main validation executable
│   └── test_phaseTransition.cc
│
├── CellStateModelCalibration/    # Analytical calibration library
│   ├── include/
│   ├── src/
│   └── test_repair_mediated_calibration.cc
│
├── PIDE3.4/                      # PIDE experimental database
├── PIDEDataReader/               # PIDE data analysis tools
├── include/                      # Main project headers
├── src/                          # Main project sources
└── python_wrapper/               # Python bindings via SWIG
```

## Modules

### CellMaker
Handles cell population layout initialization:
- 2D monolayer configurations
- 3D volume distributions
- Customizable cell spacing and dimensions
- **Geant4-independent** - can be used standalone

### ReactionDiffusionSolver
Simulates bystander signal propagation:
- 3D reaction-diffusion equations
- Configurable diffusion coefficients and reaction rates
- Signal secretion from damaged cells
- Integral concentration calculation at cell locations

### PhysicalBioTranslator
Models biological cell responses:
- Cell cycle phase progression (G1 → S → G2 → M)
- Cell state transitions (S1/S2/S3) based on damage energy
- Multi-component energy model (Er repairable + Ep persistent)
- Contact inhibition and quiescence (G0)
- Mitotic catastrophe for cells dividing with high damage
- Cell division mechanics
- **Geant4-independent** - can be used standalone

### CellStateModelCalibration
Analytical calibration library for cell state model parameters:
- 5-parameter `RepairMediatedCalibrator` (e2, e3, a, T21, T23)
- 6-parameter `RepairMediatedMisrepairCalibrator` (adds k_error)
- Nelder-Mead optimizer with negative log-likelihood objective
- Survival curve prediction from calibrated parameters

## Requirements

- **Geant4** (v11.x recommended) - for radiation transport
- **CMake** (≥3.16)
- **Python 3** - for analysis scripts and SWIG bindings
- **C++17** compatible compiler

## Building

### Full Project (with Geant4)

```bash
mkdir build && cd build
cmake -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4 ../RADCellSimulation
make -j$(nproc)
```

### Individual Modules (standalone testing)

Each module can be built independently:

```bash
# Build and test CellMaker
cd build/CellMaker
./test_cellLayoutInitializer

# Build and test ReactionDiffusionSolver
cd build/ReactionDiffusionSolver
./test_reactionDiffusion

# Build and test PhysicalBioTranslator
cd build/PhysicalBioTranslator
./test_phaseTransition
```

## Usage Examples

### Cell Layout Initialization
```cpp
#include "CellLayoutInitializer.hh"

CellLayoutInitializer layout;
layout.SetCellDimension(1.0, 1.0, 0.0);  // 2D monolayer (mm)
layout.SetGridSize(0.01);                 // 10 μm grid
layout.SetCellNumber(1000);
layout.InitializeCellPosition();

// Access cell positions
for (int i = 0; i < layout.GetCellNumber(); i++) {
    double x = layout.GetCellPositionX(i);
    double y = layout.GetCellPositionY(i);
    double z = layout.GetCellPositionZ(i);
}
```

### Cell Phase Simulation
```cpp
#include "Cell.hh"
#include "CellStateModel.hh"

CellStateModel cellModel;

// Setup cell parameters
Cell myCell;
myCell.SetCellType("fibroblast");
myCell.SetCellCycleParameter({9.0, 1.0, 11.0, 1.0, 1.0, 1.0, 1.0, 1.0});

cellModel.CellStateModelParameterSetup(myCell);
cellModel.CellPhaseInitializationRandom(cellID);

// Simulate phase progression
cellModel.CellPhaseUpdate(cellID, proliferationAllowed, deltaT, frequency);
```

### Reaction-Diffusion Simulation
```cpp
#include "DiffusionReactionSolver.hh"

DiffusionReactionSolver solver;
solver.SetGridSize(0.01);
solver.SetDiffusionCoefficient(1.0e-9);
solver.SetReactionRate(0.63e-17);

// Run simulation
solver.DiffusionReactionCalculation(deltaT);
solver.WriteConcentrationToFile(timestep, outputPath);
```

## Analysis Scripts

Python scripts for visualization are included:

```bash
# Analyze cell phase distribution
python3 PhysicalBioTranslator/plotCellPhase.py phase_output/cellPhase.csv

# Analyze concentration profiles
python3 ReactionDiffusionSolver/analyze_concentration.py
```

---

## PIDE Data Analysis Tools

The project includes Python tools for analyzing cell survival data from the GSI **PIDE (Particle Irradiation Data Ensemble) Database** version 3.4. These tools allow you to explore experimental survival fraction data for various cell lines irradiated with different ions.

### Prerequisites

Install required Python packages:

```bash
pip3 install pandas numpy matplotlib openpyxl
```

For the interactive GUI (optional):
```bash
pip3 install tkinter  # Usually pre-installed with Python
```

### Setup: Convert PIDE Excel to CSV (One-time)

Before using the analysis tools, you must convert the Excel database to CSV format:

```bash
cd RADCellSimulation/PIDE3.4
python3 convert_pide_xlsx.py
```

This creates `PIDE3.4.csv` from `PIDE3.4.xlsx`.

### Available Scripts

#### 1. `plot_pide_survival.py` - Batch Plotting (Headless)

Generates static survival curve plots from pre-generated CSV data files. Works without a display (headless/SSH sessions).

**Location:** `PIDEDataReader/plot_pide_survival.py`

**Required input files** (generated by C++ `test_pide_reader`):
- `pide_v79_12C_rawdata.csv` - Raw experimental data points
- `pide_v79_LQ_parameters.csv` - Linear-Quadratic model parameters
- `pide_v79_survival_curves.csv` - Computed survival curves

**Usage:**
```bash
# First, run the C++ reader to generate CSV files
cd build/PIDEDataReader
./test_pide_reader

# Then generate plots
python3 plot_pide_survival.py
```

**Output files:**
- `pide_survival_curves.png/pdf` - 2x2 panel with:
  - Raw experimental data (V79 + 12C)
  - LQ model curves for different ions
  - Alpha vs LET relationship
  - Data grouped by LET range
- `pide_v79_12C_survival.png` - Detailed V79+12C survival plot

#### 2. `pide_gui.py` - Interactive GUI Explorer

A graphical interface for exploring the PIDE database interactively. **Requires a display** (X11, Wayland, or similar).

**Location:** `PIDEDataReader/pide_gui.py`

**Features:**
- Select cell line and ion type
- Display raw experimental data points
- Overlay LQ model curves
- Show photon reference data
- View experiment details in a table
- Three selection modes:
  - By Cell Line → Ion
  - By Ion Type → LET
  - By Photon Type

**Usage:**
```bash
# From the PIDE3.4 directory (after running convert_pide_xlsx.py)
python3 ../PIDEDataReader/pide_gui.py PIDE3.4/

# Or with explicit path
python3 pide_gui.py /path/to/PIDE3.4/
```

**Troubleshooting GUI issues:**

If you see a Qt/xcb error like:
```
qt.qpa.plugin: Could not load the Qt platform plugin "xcb"
```

This means you're running in a headless environment. Solutions:
1. Use `plot_pide_survival.py` instead (works headless)
2. Connect via SSH with X11 forwarding: `ssh -X user@host`
3. Use a VNC or remote desktop connection

#### 3. `convert_pide_xlsx.py` - Database Converter

Converts the PIDE Excel database to CSV format for use by other tools.

**Location:** `PIDE3.4/convert_pide_xlsx.py`

**Usage:**
```bash
cd PIDE3.4
python3 convert_pide_xlsx.py
```

**Output:**
- `PIDE3.4.csv` - CSV version of the database
- Prints summary statistics (experiments, cell lines, ions, LET range)

### PIDE Data Structure

The PIDE 3.4 database contains:
- **3000+ experiments** from published literature
- **100+ cell lines** (tumor and normal)
- **Multiple ion types** (protons, helium, carbon, neon, iron, etc.)
- **LET range:** 0.2 - 1000+ keV/μm
- **LQ parameters** (α, β) for each experiment
- **Raw survival data** (dose vs. SF)

### Workflow Example

Complete workflow for analyzing V79 cell survival with carbon ions:

```bash
# 1. Setup (one-time)
cd RADCellSimulation/PIDE3.4
python3 convert_pide_xlsx.py

# 2. Build and run C++ data reader
cd ../build
make test_pide_reader
cd PIDEDataReader
./test_pide_reader

# 3. Generate plots
python3 plot_pide_survival.py

# 4. View results
ls *.png *.pdf
# pide_survival_curves.png
# pide_survival_curves.pdf  
# pide_v79_12C_survival.png
```

---

## Cell State Model Calibration

The project includes tools for calibrating cell state transition model parameters to match Linear-Quadratic (LQ) survival data.

### Location

`CellStateModelCalibration/` - Standalone calibration library
`PhysicalBioTranslator/test_calibration_validation.cc` - End-to-end validation test

### Usage

```bash
cd build/PhysicalBioTranslator

# Run calibration validation (generates survival curves and cell state plots)
./test_calibration_validation --doses 0,1,2,3,4,5,6 --alpha 0.35 --beta 0.035 --replicates 5

# Generate plots
cd calibration_validation_output
python3 plot_calibration_validation.py
```

### Output Files

- `comparison.csv` - LQ target vs. simulated survival fractions
- `calibrated_params.txt` - Fitted model parameters (E2, E3, sigma, k_error)
- `survival_fraction_comparison.pdf` - SF curve comparison plot
- `cell_state_evolution.pdf` - Cell state (S1/S2/S3) over time
- `cell_state_vs_dose.pdf` - Final cell state distribution vs. dose
- `cell_state_transition_logic.txt` - Detailed model logic documentation

---

## Cell State Transition Model with Multi-Component Energy (v1.0)

The cell state model simulates how cells transition between three states after radiation exposure:
- **S1 (Healthy)**: Normal functioning, capable of division
- **S2 (Damaged/Arrested)**: DNA damage detected, repair underway
- **S3 (Dead)**: Irreversibly damaged, removed from population

### Multi-Component Energy Model (Er + Ep)

The model distinguishes between two types of DNA damage energy:

- **Er (Repairable energy)**: Decays via bi-exponential kinetics (fast + slow DNA repair)
- **Ep (Persistent energy)**: Decays slowly with rate lambda_p, representing damage that resists normal repair

When N DSBs are induced, the injected energy is split using a nonlinear saturation function:

```
h(N) = N / (N + Nc)              -- persistent fraction
dEr_inj = (1 - h(N)) * alpha * N -- repairable component
dEp_inj = h(N) * alpha * N       -- persistent component
E_eff = Er + omega_p * Ep        -- effective energy for state transitions
```

Three cell-line-specific parameters control the Er/Ep dynamics:

| Parameter | Symbol | Description | Typical Range |
|-----------|--------|-------------|---------------|
| Half-saturation DSB count | Nc | Controls persistent damage fraction | 10-500 |
| Persistent energy weight | omega_p | Weight of Ep in effective energy | 0.3-3.0 |
| Persistent decay rate | lambda_p | Slow repair rate for Ep (1/hr) | 0.005-0.2 |

### Standardized Calibration Pipeline

The v1.0 pipeline consists of three steps:

**Step 1: Analytical Calibration (6-parameter)**

Uses `RepairMediatedMisrepairCalibrator` to fit (e2, e3, a, k_error, T21, T23) to experimental LQ survival data via Nelder-Mead optimization. The 6-param calibrator is used even though k_error is disabled in simulation, because it helps find better e2/e3/a values.

**Step 2: Er/Ep Parameter Optimization (Latin Hypercube Sampling)**

Per-cell-line optimization of (Nc, omega_p, lambda_p):
1. Generate 15 LHS samples + 1 baseline (ErEp disabled)
2. Run all simulations in parallel (5 workers)
3. Select the parameter combination with lowest error per cell line

Completes in ~5.6 minutes for 5 cell lines (9.3x faster than Optuna).

**Step 3: Fine Validation**

Production-quality simulation: 1000 cells, 5 replicates, 10 dose points (0-6 Gy).

### Running the Standard Pipeline

```bash
# Build
cd build
cmake -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4 ../RADCellSimulation
make -j$(nproc) test_multicomponent_validation

# Run LHS parameter search (from build/PhysicalBioTranslator/)
python3 lhs_surrogate_search.py

# Run fine validation with optimized parameters
./test_multicomponent_validation \
  --pide-dir /path/to/PIDE3.4 \
  --cells 1000 --replicates 5 --max-cell-types 5 \
  --use-fitted-timescales \
  --no-checkpoint --no-misrepair --force-6param-calib \
  --params-file lhs_surrogate_results/best_params_percell.csv \
  --output-dir ./validation_output
```

| Flag | Meaning |
|------|---------|
| `--use-fitted-timescales` | Use T21=T23=10h from calibration |
| `--no-checkpoint` | Disable damage diversion and S2 gating |
| `--no-misrepair` | Disable misrepair channel (k_error=0 in simulation) |
| `--force-6param-calib` | Use 6-param calibrator despite --no-misrepair |
| `--params-file` | Load per-cell-line Nc, omega_p, lambda_p from CSV |

### Validated Results (5 Cell Lines, 60Co)

| Cell Line | alpha | beta | Nc | omega_p | lambda_p | Avg Error |
|-----------|-------|------|----|---------|----------|-----------|
| CHO-10B | 0.235 | 0.031 | 20.6 | 0.380 | 0.045 | 0.037 |
| HS-23 | 0.340 | 0.028 | 208.5 | 0.845 | 0.014 | 0.047 |
| C3H10T1/2 | 0.391 | 0.021 | 483.3 | 2.914 | 0.090 | 0.086 |
| V79 | 0.189 | 0.013 | 20.6 | 0.380 | 0.045 | 0.023 |
| AG1522 | 0.337 | 0.109 | 67.2 | 0.643 | 0.046 | 0.171 |
| **Overall** | | | | | | **0.073** |

Error metric: avg |log10(SF_sim) - log10(SF_exp)| across dose points.

### Extending to New Cell Lines

1. Obtain experimental alpha and beta (LQ parameters) from PIDE database or literature
2. Add cell line to `CELL_LINES` in `lhs_surrogate_search.py` and run (~5 min per cell line)
3. Fine validation with best parameters (~3 min per cell line)
4. Target: AvgLog10Error < 0.1 (good) or < 0.15 (acceptable)

For full technical details, see `ProjectSummary_CellStateModel_v1.md`.

---

## Key Features

- **Modular Design**: Components can be developed and tested independently
- **Geant4 Integration**: Full Monte Carlo radiation transport
- **Multi-threaded**: Supports Geant4 multi-threading and OpenMP parallelization
- **Multi-Component Energy Model**: Er/Ep split with cell-line-specific persistent damage
- **Automated Calibration**: 6-param analytical calibration + LHS parameter optimization
- **Fast Pipeline**: ~8 minutes to calibrate and validate a new cell line
- **Python Bindings**: SWIG-generated Python interface
- **PIDE Integration**: Direct access to experimental survival data from the PIDE 3.4 database

## References

This simulation framework is based on research in radiation-induced bystander effects and cell state transition modeling. Key concepts include:

1. Radiation-induced bystander effect (RIBE) modeling
2. Cell cycle kinetics and phase transitions
3. Reaction-diffusion equations for signal propagation
4. DNA damage quantification and cellular response
5. Multi-component damage models (repairable vs persistent)
6. Cell state transition theory for radiation survival prediction

## Author

Rui Liu  
Email: liuruirui.nova@gmail.com

## License

This project is provided for academic and research purposes.

---

*RADCell_optimized - May 2026 (v1.0)*

