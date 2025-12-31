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

