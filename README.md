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
│   └── test_phaseTransition.cc
│
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
- Cell state transitions based on damage/stress
- Contact inhibition and quiescence (G0)
- Cell division mechanics
- **Geant4-independent** - can be used standalone

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

## Key Features

- **Modular Design**: Components can be developed and tested independently
- **Geant4 Integration**: Full Monte Carlo radiation transport
- **Multi-threaded**: Supports Geant4 multi-threading
- **Python Bindings**: SWIG-generated Python interface
- **Optimized**: Pre-computed constants and efficient algorithms

## References

This simulation framework is based on research in radiation-induced bystander effects. Key concepts include:

1. Radiation-induced bystander effect (RIBE) modeling
2. Cell cycle kinetics and phase transitions
3. Reaction-diffusion equations for signal propagation
4. DNA damage quantification and cellular response

## Author

Rui Liu  
Email: liuruirui.nova@gmail.com

## License

This project is provided for academic and research purposes.

---

*RADCell_optimized - December 2024*

