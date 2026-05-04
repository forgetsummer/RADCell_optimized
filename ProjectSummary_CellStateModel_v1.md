# Cell State Transition Model -- Project Summary (v1.0)

## Overview

This document summarizes the current state of the Cell State Transition Model for radiation cell survival prediction. The model simulates how cells transition between healthy (S1), damaged (S2), and dead (S3) states after radiation exposure, incorporating a multi-component energy model (Er/Ep) that distinguishes between repairable and persistent DNA damage. A standardized calibration and simulation pipeline has been established and validated against experimental data for five 60Co-irradiated cell lines.

---

## Model Architecture

### Cell State Transition Framework

Cells exist in one of three states:
- **S1 (Healthy)**: Normal functioning, capable of division
- **S2 (Damaged/Arrested)**: DNA damage detected, repair underway
- **S3 (Dead)**: Irreversibly damaged, removed from population

State transitions are governed by overlap scores between the cell's current energy and the characteristic energies E1=0, E2, E3 of each state, computed as:

```
w(E, Ej) = 2 * Phi(-|E - Ej| / (2*sigma))
```

where Phi is the standard normal CDF. Transition probabilities are normalized: P_j = w_j / (w1 + w2 + w3).

### Multi-Component Energy Model (Er + Ep)

The total cell state energy is split into two components:

- **Er (Repairable energy)**: Decays via bi-exponential kinetics using DNA repair rates (fast and slow repair fractions f1, f2 with rates lambda1, lambda2)
- **Ep (Persistent energy)**: Decays via mono-exponential kinetics with a slow rate lambda_p

When new DSBs are induced (N DSBs from dose D), the injected energy is split using a nonlinear saturation function:

```
h(N) = N / (N + Nc)          -- persistent fraction
dEr_inj = (1 - h(N)) * alpha * N
dEp_inj = h(N) * alpha * N
```

where Nc is the half-saturation DSB count (cell-line specific).

The effective energy used for S2 transition rates and mitotic catastrophe checks is:

```
E_eff = Er + omega_p * Ep
```

where omega_p is the persistent energy weight (cell-line specific).

### New Parameters (Er/Ep)

| Parameter | Symbol | Description | Range |
|-----------|--------|-------------|-------|
| Half-saturation DSB count | Nc | Controls how much damage becomes persistent | 10-500 |
| Persistent energy weight | omega_p | Weight of Ep in effective energy | 0.3-3.0 |
| Persistent decay rate | lambda_p | Slow repair rate for Ep (1/hr) | 0.005-0.2 |

These three parameters are **cell-line specific** and determined via Latin Hypercube Sampling optimization (see below).

---

## Standardized Simulation Pipeline (v1.0)

### Step 1: Analytical Calibration (6-parameter)

The `RepairMediatedMisrepairCalibrator` fits six parameters to experimental LQ survival data:

| Parameter | Description |
|-----------|-------------|
| e2 | Reduced S2 energy boundary (E2/sigma) |
| e3 | Reduced S3 energy boundary (E3/sigma) |
| a | Reduced dose-response coefficient (alpha/sigma) |
| k_error | Misrepair hazard coefficient |
| T21 | S2->S1 transition timescale |
| T23 | S2->S3 transition timescale |

The calibrator uses Nelder-Mead optimization to minimize negative log-likelihood against synthetic LQ data points generated from experimental alpha/beta values.

**Configuration:**
- kappa = 40.0 DSB/Gy
- dose_max_fit = 6.0 Gy
- T_assay = 27.8 hours
- T21 = T23 = T_cellCycle = 10.0 hours (fixed)
- sigma = 10.0 (chosen for numerical stability)

**Key design decision:** The 6-parameter calibrator (including k_error) is used even though the simulation runs with k_error=0. This is because k_error provides an additional degree of freedom that helps the optimizer find e2/e3/a values that translate better into the full simulation. Without k_error, the optimizer distorts e2/e3/a to compensate, producing parameters that diverge badly in simulation.

### Step 2: Er/Ep Parameter Optimization (Latin Hypercube Sampling)

Per-cell-line optimization of (Nc, omega_p, lambda_p) using LHS + best-observed-point selection.

**Method:**
1. Generate 15 Latin Hypercube samples in (log10(Nc), omega_p, log10(lambda_p)) space
2. Add 1 baseline point with ErEp disabled (Nc=1e9, omega_p=1, lambda_p=0)
3. Run all 80 simulations (16 per cell line x 5 cell lines) in parallel
4. For each cell line, select the parameter combination with lowest AvgLog10Error

**Search ranges:**
- Nc: 10 to 500 (log scale)
- omega_p: 0.3 to 3.0 (linear)
- lambda_p: 0.005 to 0.2 (log scale)

**Performance:** 80 simulations completed in **5.6 minutes** (with 5 parallel workers, OMP_NUM_THREADS=22 each), compared to ~52 minutes for Optuna (30 trials x 5 cell lines, sequential). This is a **9.3x speedup**.

**Script:** `lhs_surrogate_search.py`

### Step 3: Fine Validation

Production-quality simulation with best parameters:
- 1000 cells per simulation
- 5 replicates per dose point
- Dose points: 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0 Gy

**Executable:** `test_multicomponent_validation`

### Simulation Flags (v1.0 Standard)

```bash
./test_multicomponent_validation \
  --pide-dir <PIDE3.4_path> \
  --cells 1000 --replicates 5 --max-cell-types 5 \
  --use-fitted-timescales \
  --no-checkpoint --no-misrepair --force-6param-calib \
  --params-file <best_params_percell.csv> \
  --output-dir <output_directory>
```

| Flag | Meaning |
|------|---------|
| `--use-fitted-timescales` | Use T21=T23=10h from calibration (not effective repair time) |
| `--no-checkpoint` | Disable damage diversion and S2 gating |
| `--no-misrepair` | Disable misrepair channel (k_error=0 in simulation) |
| `--force-6param-calib` | Use 6-param calibrator despite --no-misrepair |
| `--params-file` | Load per-cell-line Nc, omega_p, lambda_p from CSV |

### Active Mechanisms in v1.0 Simulation

| Mechanism | Status | Description |
|-----------|--------|-------------|
| Overlap-based state assignment | ON | S1/S2/S3 via Gaussian overlaps |
| S2 competing risks (recovery/death) | ON | lambda_21, lambda_23 from overlap scores |
| Er bi-exponential decay | ON | Fast+slow DNA repair kinetics |
| Ep mono-exponential decay | ON | Persistent damage with rate lambda_p |
| E_eff = Er + omega_p * Ep | ON | Weighted effective energy for transitions |
| h(N) = N/(N+Nc) | ON | Nonlinear persistent fraction |
| Mitotic catastrophe | ON | Cells in M-phase with high energy die |
| Contact inhibition | ON | Cells stop dividing when surrounded |
| Checkpoint (damage diversion) | OFF | Not used in v1.0 |
| Checkpoint (S2 gating) | OFF | Not used in v1.0 |
| Misrepair (k_error death channel) | OFF | Not used in v1.0 |

---

## Validated Results

### Calibrated Parameters (6-parameter analytical)

| Cell Line | e2 | e3 | a | k_error | T21 | T23 |
|-----------|------|------|--------|---------|-----|-----|
| CHO-10B | 3.004 | 5.766 | 0.0197 | 2.06e-4 | 10 | 10 |
| HS-23 | 2.510 | 5.721 | 0.0182 | 3.26e-3 | 10 | 10 |
| C3H10T1/2 | 2.467 | 6.101 | 0.0182 | 9.66e-3 | 10 | 10 |
| V79 | 2.888 | 6.164 | 0.0165 | 1.06e-4 | 10 | 10 |
| AG1522 | 3.086 | 6.557 | 0.0305 | 1.13e-3 | 10 | 10 |

Note: k_error values are calibrated but NOT used in the simulation (set to 0).

### Optimized Er/Ep Parameters (LHS)

| Cell Line | Nc | omega_p | lambda_p | Ep half-life (h) |
|-----------|------|---------|----------|-------------------|
| CHO-10B | 20.6 | 0.380 | 0.045 | 15.4 |
| HS-23 | 208.5 | 0.845 | 0.014 | 49.5 |
| C3H10T1/2 | 483.3 | 2.914 | 0.090 | 7.7 |
| V79 | 20.6 | 0.380 | 0.045 | 15.4 |
| AG1522 | 67.2 | 0.643 | 0.046 | 15.1 |

### Simulation Performance (avg |log10(SF)| error)

| Cell Line | alpha | beta | alpha/beta | Baseline (no ErEp) | Optuna (ckpt+misrep ON) | LHS v1.0 |
|-----------|-------|------|------------|-----------|--------|----------|
| CHO-10B | 0.235 | 0.031 | 7.7 | 0.125 | 0.057 | **0.037** |
| HS-23 | 0.340 | 0.028 | 12.0 | 0.218 | 0.055 | **0.047** |
| C3H10T1/2 | 0.391 | 0.021 | 18.8 | 0.271 | 0.040 | **0.086** |
| V79 | 0.189 | 0.013 | 14.6 | 0.059 | 0.059 | **0.023** |
| AG1522 | 0.337 | 0.109 | 3.1 | 0.510 | 0.129 | **0.171** |
| **Overall** | | | | **0.236** | **0.068** | **0.073** |

- **Baseline**: 6-param calibration, checkpoints ON, misrepair ON, no ErEp
- **Optuna**: 6-param calibration, checkpoints ON, misrepair ON, per-cell-line ErEp (full model)
- **LHS v1.0**: 6-param calibration, checkpoints OFF, misrepair OFF, per-cell-line ErEp (simplified model)

The v1.0 simplified model achieves overall error 0.073, comparable to the full-model Optuna result (0.068), with a much cleaner model architecture. CHO-10B and V79 (previously the most problematic cell lines) show the best results (0.037 and 0.023).

---

## Computational Performance

### LHS Parameter Search

| Metric | Value |
|--------|-------|
| Total simulations | 80 (15 LHS + 1 baseline per cell line x 5 lines) |
| Wall time | 5.6 minutes |
| Per-simulation time | 20.4s avg (min 19.8s, max 21.2s) |
| Parallelism | 5 workers, OMP_NUM_THREADS=22 each |
| Speedup vs Optuna | 9.3x |

### Fine Validation Run

| Metric | Value |
|--------|-------|
| Cells per simulation | 1000 |
| Replicates per dose | 5 |
| Dose points | 10 (0 to 6 Gy) |
| Cell lines | 5 |
| Wall time | 859 seconds (~14.3 minutes) |

### Total Pipeline Time (for a new cell line)

| Step | Time |
|------|------|
| Analytical calibration | < 1 second |
| LHS search (16 sims, 500 cells, 1 rep) | ~5 minutes |
| Fine validation (1000 cells, 5 reps) | ~3 minutes per cell line |
| **Total** | **~8-9 minutes per new cell line** |

---

## File Locations

### Source Code

| File | Description |
|------|-------------|
| `PhysicalBioTranslator/src/CellStateModel.cc` | Core simulation: energy update, state transitions, Er/Ep split |
| `PhysicalBioTranslator/include/CellStateModel.hh` | Model header with StateInfo (Er, Ep), CellStateParaInfo (Nc, omega_p, lambda_p) |
| `PhysicalBioTranslator/test_multicomponent_validation.cc` | Main test executable with CLI flags |
| `CellStateModelCalibration/src/RepairMediatedModel_with_k_error.cc` | 6-param analytical survival model |
| `CellStateModelCalibration/src/RepairMediatedCalibrator_with_k_error.cc` | 6-param Nelder-Mead calibrator |

### Scripts

| File | Description |
|------|-------------|
| `RAD-build/PhysicalBioTranslator/lhs_surrogate_search.py` | LHS parameter optimization script |
| `RAD-build/PhysicalBioTranslator/plot_lhs_comparison.py` | Comparison plotting script |

### Results

| Directory | Description |
|-----------|-------------|
| `RAD-build/PhysicalBioTranslator/lhs_surrogate_results/` | LHS search output (80 sims, timing, best params CSV) |
| `RAD-build/PhysicalBioTranslator/lhs_nockpt_fine/` | Fine validation results (1000 cells, 5 reps) |
| `RAD-build/PhysicalBioTranslator/lhs_comparison_plots/` | SF comparison plots |
| `RAD-build/PhysicalBioTranslator/grid_fittedT_baseline_fine/` | Previous baseline results (no ErEp) |
| `RAD-build/PhysicalBioTranslator/optuna_percell_fine/` | Previous Optuna results (full model) |

### Configuration Files

| File | Description |
|------|-------------|
| `lhs_surrogate_results/best_params_percell.csv` | Per-cell-line (Nc, omega_p, lambda_p) for v1.0 |

---

## Standard Procedure for Extending to New Cell Lines

To add a new cell line to the model:

1. **Obtain experimental data**: alpha and beta from LQ fit to 60Co survival data (from PIDE database or literature).

2. **Run LHS parameter search**:
   ```bash
   # Edit lhs_surrogate_search.py to add the new cell line to CELL_LINES list
   # Or run test_multicomponent_validation directly:
   ./test_multicomponent_validation \
     --pide-dir <PIDE3.4_path> \
     --cell-line "<CellLineName>" \
     --cells 500 --replicates 1 \
     --use-fitted-timescales \
     --no-checkpoint --no-misrepair --force-6param-calib \
     --nc <Nc_value> --omega-p <wp_value> --lambda-p <lp_value> \
     --output-dir <trial_output>
   ```
   Run for each of 15+ LHS sample points and select the best.

3. **Fine validation**:
   ```bash
   ./test_multicomponent_validation \
     --pide-dir <PIDE3.4_path> \
     --cell-line "<CellLineName>" \
     --cells 1000 --replicates 5 \
     --use-fitted-timescales \
     --no-checkpoint --no-misrepair --force-6param-calib \
     --nc <best_Nc> --omega-p <best_wp> --lambda-p <best_lp> \
     --output-dir <fine_output>
   ```

4. **Verify**: Check that AvgLog10Error < 0.1 (good) or < 0.15 (acceptable).

---

## Known Limitations and Future Work

1. **AG1522 remains challenging** (error 0.171): This cell line has an unusually high beta (0.109), creating steep high-dose curvature that the model struggles to reproduce. May require model refinements or a dedicated calibration approach.

2. **Analytical-simulation gap**: The analytical calibration model is structurally simpler than the simulation (no mitotic catastrophe, no contact inhibition, no continuous energy re-evaluation). The 6-param calibrator's k_error parameter partially compensates for these unmodeled mechanisms.

3. **Surrogate model extrapolation**: The quadratic response surface tends to extrapolate to unrealistic (negative) errors at parameter boundaries. The current approach falls back to best-observed LHS point in these cases. A Gaussian Process surrogate or adaptive sampling could improve this.

4. **Checkpoint and misrepair mechanisms**: These are implemented but disabled in v1.0. They could be re-enabled for specific applications where their biological effects are important.

5. **Cell cycle phase effects**: The simulation includes phase-specific energy scaling (fS, fG2, fM) that the analytical model does not account for. This contributes to the analytical-simulation gap.

---

*Summary prepared: May 2026*
*Build directory: /home/user/Geant4Projects/CellStateTransitionPaper/RAD-build/*
*Source directory: /home/user/Geant4Projects/CellStateTransitionPaper/RADCellSimulation/*
