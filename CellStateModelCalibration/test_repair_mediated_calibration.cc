/**
 * @file test_repair_mediated_calibration.cc
 * @brief Test program for Repair Mediated model calibration (S1->S2->S3)
 * 
 * This test demonstrates fitting the 3-parameter Repair Mediated model
 * to experimental cell survival data. The model includes:
 * - S1 (healthy): Undamaged or lightly damaged cells
 * - S2 (repair): Cells with moderate damage undergoing repair
 * - S3 (dead): Cells with lethal damage
 * 
 * Fitted Parameters (from survival data):
 * - e2 = E2/sigma: Repair state threshold
 * - e3 = E3/sigma: Death state threshold
 * - a = alpha/sigma: Damage per DSB
 * 
 * Fixed Parameters (set before fitting):
 * - T21: Repair timescale (S2->S1) - set to T_cellCycle
 * - T23: Death timescale (S2->S3) - set to T_cellCycle
 */

#include "RepairMediatedCalibrator.hh"
#include "DataTypes.hh"
#include "MathUtilities.hh"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

using namespace CellStateCalibration;

int main() {
    std::cout << "========================================\n";
    std::cout << "Repair Mediated Model Calibration Test\n";
    std::cout << "(S1 -> S2 -> S3 transitions)\n";
    std::cout << "========================================\n\n";
    
    // =====================================================
    // STEP 1: Define biological timescales (FIXED INPUTS)
    // =====================================================
    
    // Cell cycle time - this is the key biological timescale
    // Typical values: 10-24 hours depending on cell type
    const double T_cellCycle = 10.0;  // hours
    
    std::cout << "=== FIXED Biological Timescales ===\n";
    std::cout << "  T_cellCycle = " << T_cellCycle << " hours\n";
    std::cout << "  T21 (repair) = T_cellCycle = " << T_cellCycle << " hours\n";
    std::cout << "  T23 (death)  = T_cellCycle = " << T_cellCycle << " hours\n\n";
    
    // =====================================================
    // STEP 2: Configure the calibration
    // =====================================================
    
    FitConfigS2 cfg;
    cfg.kappa = 40.0;           // DSB/Gy
    cfg.dose_max_fit = 6.0;     // Fit all data
    cfg.T_assay_h = 24.0 * 10.0; // 10 days for clonogenic assay
    
    // Set T21 and T23 to cell cycle time (FIXED, not fitted!)
    cfg.T21_fixed = T_cellCycle;
    cfg.T23_fixed = T_cellCycle;
    
    std::cout << "=== Configuration ===\n";
    std::cout << "  kappa (DSB/Gy): " << cfg.kappa << "\n";
    std::cout << "  dose_max_fit: " << cfg.dose_max_fit << " Gy\n";
    std::cout << "  T_assay: " << cfg.T_assay_h << " hours (" << cfg.T_assay_h/24.0 << " days)\n";
    std::cout << "  T21_fixed: " << cfg.T21_fixed << " hours (FIXED)\n";
    std::cout << "  T23_fixed: " << cfg.T23_fixed << " hours (FIXED)\n\n";
    
    // =====================================================
    // STEP 3: Load experimental data
    // =====================================================
    
    // Example experimental data (replace with your data)
    std::vector<DataPoint> data = {
        {0.0, 1.0, 0.0},
        {1.0, 0.70, 0.0},
        {2.0, 0.45, 0.0},
        {3.0, 0.25, 0.0},
        {4.0, 0.12, 0.0},
        {6.0, 0.02, 0.0}
    };
    
    std::cout << "=== Experimental Data ===\n";
    std::cout << " Dose (Gy)        SF_obs\n";
    std::cout << "  -------------------------\n";
    for (const auto& dp : data) {
        std::cout << std::setw(10) << dp.D << std::setw(15) << dp.SF_obs << "\n";
    }
    std::cout << "\n";
    
    // =====================================================
    // STEP 4: Create calibrator and fit
    // =====================================================
    
    RepairMediatedCalibrator calibrator(cfg);
    calibrator.setData(data);
    calibrator.setDefaultUncertainties(0.2);
    
    // Initial guess for the 3 fitted parameters (e2, e3, a)
    // T21 and T23 are NOT part of the optimization anymore!
    ParamsS2 x0(1.5, 3.0, 0.02);
    ParamsS2 step(0.5, 1.0, 0.01);
    
    std::cout << "=== Fitting Model ===\n";
    std::cout << "Fitting 3 parameters: e2, e3, a\n";
    std::cout << "T21 and T23 are FIXED at " << T_cellCycle << " hours\n";
    std::cout << "Initial guess: e2=" << x0.e2 << ", e3=" << x0.e3 << ", a=" << x0.a << "\n\n";
    
    // Fit
    RepairMediatedFitResult result = calibrator.fit(x0, step);
    
    // =====================================================
    // STEP 5: Display results
    // =====================================================
    
    std::cout << "=== Fit Results ===\n\n";
    std::cout << "Fitted Parameters (from survival data):\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  e2  = E2/sigma = " << result.params.e2 << " (repair threshold)\n";
    std::cout << "  e3  = E3/sigma = " << result.params.e3 << " (death threshold)\n";
    std::cout << "  a   = alpha/sigma = " << result.params.a << " (damage per DSB)\n\n";
    
    std::cout << "Fixed Parameters (from config):\n";
    std::cout << "  T21 = " << cfg.T21_fixed << " hours (repair timescale)\n";
    std::cout << "  T23 = " << cfg.T23_fixed << " hours (death timescale)\n\n";
    
    std::cout << "Optimization:\n";
    std::cout << "  Negative log-likelihood: " << std::setprecision(6) << result.negLogLikelihood << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n\n";
    
    // Physical interpretation
    std::cout << "Physical Interpretation:\n";
    std::cout << "  - Repair DSB threshold (e2/a): " << std::setprecision(2) 
              << result.params.e2 / result.params.a << " DSBs\n";
    std::cout << "  - Death DSB threshold (e3/a): " 
              << result.params.e3 / result.params.a << " DSBs\n";
    std::cout << "  - Repair dose threshold: " 
              << (result.params.e2 / result.params.a) / cfg.kappa << " Gy\n";
    std::cout << "  - Death dose threshold: " 
              << (result.params.e3 / result.params.a) / cfg.kappa << " Gy\n\n";
    
    // =====================================================
    // STEP 6: Show how to get physical parameters
    // =====================================================
    
    std::cout << "========================================\n";
    std::cout << "CALIBRATED PARAMETERS FOR CELL STATE MODEL\n";
    std::cout << "========================================\n\n";
    
    std::cout << "To use these in your CellStateModel, choose sigma:\n\n";
    
    double sigma_example = 10.0;
    std::cout << ">>> Example with sigma = " << sigma_example << " <<<\n";
    CellStateModelParamsS2 params = CellStateModelParamsS2::fromReduced(
        result.params, sigma_example, cfg.kappa);
    params.print();
    
    std::cout << "\n--- Copy these values to your Cell setup ---\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "E1    = " << params.E1 << "\n";
    std::cout << "E2    = " << params.E2 << "\n";
    std::cout << "E3    = " << params.E3 << "\n";
    std::cout << "sigma = " << params.sigma << "\n";
    std::cout << "alpha = " << params.alpha << "\n";
    std::cout << "kappa = " << params.kappa << " DSB/Gy\n";
    std::cout << "T21   = " << params.T21 << " hours (from config, = T_cellCycle)\n";
    std::cout << "T23   = " << params.T23 << " hours (from config, = T_cellCycle)\n\n";
    
    std::cout << "--- For other sigma values, use these formulas ---\n";
    std::cout << "E1    = 0\n";
    std::cout << "E2    = " << result.params.e2 << " * sigma\n";
    std::cout << "E3    = " << result.params.e3 << " * sigma\n";
    std::cout << "alpha = " << result.params.a << " * sigma\n";
    std::cout << "T21   = " << cfg.T21_fixed << " hours (FIXED, independent of sigma)\n";
    std::cout << "T23   = " << cfg.T23_fixed << " hours (FIXED, independent of sigma)\n\n";
    
    // =====================================================
    // STEP 7: Model predictions vs data
    // =====================================================
    
    std::cout << "=== Model Predictions vs Data ===\n";
    std::cout << " Dose (Gy)      SF_obs    SF_model   Residual\n";
    std::cout << "  ------------------------------------------------\n";
    for (const auto& dp : data) {
        double sf_model = calibrator.getModel().survivalFraction(dp.D, result.params);
        double residual = std::log(dp.SF_obs + 1e-12) - std::log(sf_model);
        std::cout << std::setw(10) << std::setprecision(2) << dp.D
                  << std::setw(12) << std::setprecision(4) << dp.SF_obs
                  << std::setw(12) << sf_model
                  << std::setw(11) << residual << "\n";
    }
    std::cout << "\n";
    
    // Generate predicted survival curve
    std::cout << "=== Predicted Survival Curve ===\n";
    std::cout << " Dose (Gy)  SF_predicted\n";
    std::cout << "  -------------------------\n";
    std::vector<double> doses;
    for (double d = 0.0; d <= 8.0; d += 0.5) {
        doses.push_back(d);
    }
    std::vector<double> sf_pred = calibrator.predictSurvival(result.params, doses);
    for (size_t i = 0; i < doses.size(); ++i) {
        std::cout << std::setw(10) << std::setprecision(2) << doses[i]
                  << std::setw(14) << std::setprecision(6) << sf_pred[i] << "\n";
    }
    std::cout << "\n";
    
    // Save results to CSV
    {
        std::ofstream fout("repair_mediated_results.csv");
        fout << "Dose_Gy,SF_obs,SF_model\n";
        for (const auto& dp : data) {
            double sf_model = calibrator.getModel().survivalFraction(dp.D, result.params);
            fout << dp.D << "," << dp.SF_obs << "," << sf_model << "\n";
        }
        fout.close();
    }
    
    {
        std::ofstream fout("repair_mediated_curve.csv");
        fout << "Dose_Gy,SF_predicted\n";
        for (size_t i = 0; i < doses.size(); ++i) {
            fout << doses[i] << "," << sf_pred[i] << "\n";
        }
        fout.close();
    }
    
    std::cout << "Results saved to:\n";
    std::cout << "  - repair_mediated_results.csv\n";
    std::cout << "  - repair_mediated_curve.csv\n\n";
    
    // Summary
    std::cout << "========================================\n";
    std::cout << "SUMMARY\n";
    std::cout << "========================================\n";
    std::cout << "This calibration uses:\n";
    std::cout << "  - FIXED timescales: T21 = T23 = T_cellCycle = " << T_cellCycle << " hours\n";
    std::cout << "  - FITTED parameters: e2, e3, a (3 parameters)\n\n";
    std::cout << "The rationale is that T21 and T23 represent biological\n";
    std::cout << "timescales that can be measured directly (cell cycle time),\n";
    std::cout << "while e2, e3, a are fitted from survival data.\n\n";
    std::cout << "========================================\n";
    std::cout << "Calibration test completed successfully!\n";
    std::cout << "========================================\n";
    
    return 0;
}
