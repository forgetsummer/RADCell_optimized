/**
 * @file test_misrepair_calibration.cc
 * @brief Test program for Repair Mediated model calibration with stochastic misrepair
 *
 * This test demonstrates fitting the Repair Mediated model with stochastic misrepair
 * to experimental cell survival data. The model includes:
 * - S1 (healthy): Undamaged or lightly damaged cells
 * - S2 (repair): Cells with moderate damage undergoing repair
 * - S3 (dead): Cells with lethal damage
 *
 * IMPORTANT: T21 and T23 are now FIXED by user (not optimized)
 * This prevents parameter inflation issues where the optimizer would
 * exploit large T21/T23 values to artificially increase survival.
 *
 * Fixed Parameters (user-defined):
 * - T21: Repair timescale (S2->S1) = 12 hours
 * - T23: Death timescale (S2->S3) = 24 hours
 *
 * Fitted Parameters (4 from survival data):
 * - e2 = E2/sigma: Repair state threshold
 * - e3 = E3/sigma: Death state threshold
 * - a = alpha/sigma: Damage per DSB
 * - k_error: Misrepair rate per (DSB*hour), contributing to S2->S3
 *
 * The stochastic misrepair channel adds a DSB-dependent death hazard:
 *   lambda_mis = k_error * N_dsb
 * This captures the biological reality that cells with more DSBs have
 * higher probability of misrepair events leading to cell death.
 */

#include "RepairMediatedCalibrator_with_k_error.hh"
#include "RepairMediatedModel_with_k_error.hh"
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
    std::cout << "Stochastic Misrepair Model Calibration Test\n";
    std::cout << "(S1 -> S2 -> S3 with k_error parameter)\n";
    std::cout << "T21 and T23 FIXED by user for stable fitting\n";
    std::cout << "========================================\n\n";

    // =====================================================
    // STEP 1: Configure the calibration
    // =====================================================

    FitConfigS2 cfg;
    cfg.kappa = 40.0;           // DSB/Gy
    cfg.dose_max_fit = 6.0;     // Fit all data
    cfg.T_assay_h = 24.0 * 10.0; // 10 days for clonogenic assay

    // FIX T21 and T23 to user-defined values
    // This prevents parameter inflation issues in the optimizer
    cfg.fixTimescales(12.0, 24.0);  // T21 = 12 hours, T23 = 24 hours

    std::cout << "=== Configuration ===\n";
    std::cout << "  kappa (DSB/Gy): " << cfg.kappa << "\n";
    std::cout << "  dose_max_fit: " << cfg.dose_max_fit << " Gy\n";
    std::cout << "  T_assay: " << cfg.T_assay_h << " hours (" << cfg.T_assay_h/24.0 << " days)\n";
    std::cout << "  T21 (FIXED): " << cfg.T21_fixed << " hours\n";
    std::cout << "  T23 (FIXED): " << cfg.T23_fixed << " hours\n\n";

    // =====================================================
    // STEP 2: Load experimental data
    // =====================================================

    // Example experimental data (replace with your data)
    // These data points represent typical cell survival curve
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
    // STEP 3: Create calibrator and fit
    // =====================================================

    RepairMediatedMisrepairCalibrator calibrator(cfg);
    calibrator.setData(data);
    calibrator.setDefaultUncertainties(0.2);

    // Initial guess for 4 fitted parameters (T21/T23 are FIXED):
    // (e2, e3, a, k_error)
    ParamsS2 x0;
    x0.e2 = 1.5;        // Repair threshold
    x0.e3 = 3.0;        // Death threshold
    x0.a = 0.02;        // Damage per DSB
    x0.T21 = cfg.T21_fixed;  // Will use fixed value
    x0.T23 = cfg.T23_fixed;  // Will use fixed value
    x0.k_error = 0.001; // Misrepair rate per (DSB*hour)

    ParamsS2 step;
    step.e2 = 0.5;
    step.e3 = 1.0;
    step.a = 0.01;
    step.T21 = 0.0;     // Not used (T21 is fixed)
    step.T23 = 0.0;     // Not used (T23 is fixed)
    step.k_error = 0.0005;

    std::cout << "=== Fitting Model ===\n";
    std::cout << "Fitting 4 parameters: e2, e3, a, k_error\n";
    std::cout << "T21 = " << cfg.T21_fixed << " hours (FIXED by user)\n";
    std::cout << "T23 = " << cfg.T23_fixed << " hours (FIXED by user)\n";
    std::cout << "Initial guess:\n";
    std::cout << "  e2 = " << x0.e2 << "\n";
    std::cout << "  e3 = " << x0.e3 << "\n";
    std::cout << "  a  = " << x0.a << "\n";
    std::cout << "  k_error = " << x0.k_error << " per (DSB*hour)\n\n";

    // Fit
    RepairMediatedMisrepairFitResult result = calibrator.fit(x0, step);

    // =====================================================
    // STEP 4: Display results
    // =====================================================

    std::cout << "=== Fit Results ===\n\n";
    std::cout << "Fitted Parameters:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  e2      = E2/sigma   = " << result.params.e2 << " (repair threshold) - FITTED\n";
    std::cout << "  e3      = E3/sigma   = " << result.params.e3 << " (death threshold) - FITTED\n";
    std::cout << "  a       = alpha/sigma = " << result.params.a << " (damage per DSB) - FITTED\n";
    std::cout << "  k_error = " << result.params.k_error << " per (DSB*hour) (misrepair rate) - FITTED\n";
    std::cout << "Fixed Parameters:\n";
    std::cout << "  T21     = " << result.params.T21 << " hours (repair timescale) - USER FIXED\n";
    std::cout << "  T23     = " << result.params.T23 << " hours (death timescale) - USER FIXED\n\n";

    std::cout << "Optimization:\n";
    std::cout << "  Negative log-likelihood: " << result.negLogLikelihood << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n\n";

    // Physical interpretation
    std::cout << "Physical Interpretation:\n";
    std::cout << std::setprecision(2);
    if (result.params.a > 1e-9) {
        std::cout << "  - Repair DSB threshold (e2/a): "
                  << result.params.e2 / result.params.a << " DSBs\n";
        std::cout << "  - Death DSB threshold (e3/a): "
                  << result.params.e3 / result.params.a << " DSBs\n";
        std::cout << "  - Repair dose threshold: "
                  << (result.params.e2 / result.params.a) / cfg.kappa << " Gy\n";
        std::cout << "  - Death dose threshold: "
                  << (result.params.e3 / result.params.a) / cfg.kappa << " Gy\n";
    }
    std::cout << "  - Repair half-time: " << std::log(2.0) / (1.0 / result.params.T21) << " hours\n";
    std::cout << "  - Death half-time: " << std::log(2.0) / (1.0 / result.params.T23) << " hours\n";
    std::cout << "\n";

    // Misrepair interpretation
    std::cout << "Misrepair Channel:\n";
    std::cout << "  - k_error = " << std::setprecision(6) << result.params.k_error << " per (DSB*hour)\n";
    std::cout << "  - For a cell with 10 DSBs:\n";
    std::cout << "    lambda_mis = k_error * 10 = " << result.params.k_error * 10.0 << " per hour\n";
    std::cout << "    P(misrepair death in 1h) ~ " << std::setprecision(4)
              << (1.0 - std::exp(-result.params.k_error * 10.0)) << "\n";
    std::cout << "  - For a cell with 50 DSBs:\n";
    std::cout << "    lambda_mis = k_error * 50 = " << std::setprecision(6) << result.params.k_error * 50.0 << " per hour\n";
    std::cout << "    P(misrepair death in 1h) ~ " << std::setprecision(4)
              << (1.0 - std::exp(-result.params.k_error * 50.0)) << "\n\n";

    // =====================================================
    // STEP 5: Show how to get physical parameters
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
    std::cout << "E1      = " << params.E1 << "\n";
    std::cout << "E2      = " << params.E2 << "\n";
    std::cout << "E3      = " << params.E3 << "\n";
    std::cout << "sigma   = " << params.sigma << "\n";
    std::cout << "alpha   = " << params.alpha << "\n";
    std::cout << "kappa   = " << params.kappa << " DSB/Gy\n";
    std::cout << "T21     = " << params.T21 << " hours (from optimization)\n";
    std::cout << "T23     = " << params.T23 << " hours (from optimization)\n";
    std::cout << "k_error = " << result.params.k_error << " per (DSB*hour)\n\n";

    std::cout << "--- For other sigma values, use these formulas ---\n";
    std::cout << "E1      = 0\n";
    std::cout << "E2      = " << result.params.e2 << " * sigma\n";
    std::cout << "E3      = " << result.params.e3 << " * sigma\n";
    std::cout << "alpha   = " << result.params.a << " * sigma\n";
    std::cout << "T21     = " << result.params.T21 << " hours (USER FIXED, independent of sigma)\n";
    std::cout << "T23     = " << result.params.T23 << " hours (USER FIXED, independent of sigma)\n";
    std::cout << "k_error = " << result.params.k_error << " per (DSB*hour) (independent of sigma)\n\n";

    // =====================================================
    // STEP 6: Model predictions vs data
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
        std::ofstream fout("misrepair_calibration_results.csv");
        fout << "Dose_Gy,SF_obs,SF_model\n";
        for (const auto& dp : data) {
            double sf_model = calibrator.getModel().survivalFraction(dp.D, result.params);
            fout << dp.D << "," << dp.SF_obs << "," << sf_model << "\n";
        }
        fout.close();
    }

    {
        std::ofstream fout("misrepair_calibration_curve.csv");
        fout << "Dose_Gy,SF_predicted\n";
        for (size_t i = 0; i < doses.size(); ++i) {
            fout << doses[i] << "," << sf_pred[i] << "\n";
        }
        fout.close();
    }

    // Save fitted parameters
    {
        std::ofstream fout("misrepair_calibration_params.csv");
        fout << "Parameter,Value,Description\n";
        fout << "e2," << result.params.e2 << ",E2/sigma (repair threshold)\n";
        fout << "e3," << result.params.e3 << ",E3/sigma (death threshold)\n";
        fout << "a," << result.params.a << ",alpha/sigma (damage per DSB)\n";
        fout << "T21," << result.params.T21 << ",Repair timescale (hours)\n";
        fout << "T23," << result.params.T23 << ",Death timescale (hours)\n";
        fout << "k_error," << result.params.k_error << ",Misrepair rate per (DSB*hour)\n";
        fout << "kappa," << cfg.kappa << ",DSB per Gy\n";
        fout << "T_assay," << cfg.T_assay_h << ",Assay time (hours)\n";
        fout << "NLL," << result.negLogLikelihood << ",Negative log-likelihood\n";
        fout << "converged," << (result.converged ? 1 : 0) << ",Convergence flag\n";
        fout.close();
    }

    std::cout << "Results saved to:\n";
    std::cout << "  - misrepair_calibration_results.csv\n";
    std::cout << "  - misrepair_calibration_curve.csv\n";
    std::cout << "  - misrepair_calibration_params.csv\n\n";

    // =====================================================
    // STEP 7: Compare with model without misrepair (k_error = 0)
    // =====================================================

    std::cout << "========================================\n";
    std::cout << "COMPARISON: With vs Without Misrepair\n";
    std::cout << "========================================\n\n";

    // Create parameters with k_error = 0
    ParamsS2 params_no_misrepair = result.params;
    params_no_misrepair.k_error = 0.0;

    std::cout << " Dose (Gy)  SF_with_k   SF_without  Difference\n";
    std::cout << "  ------------------------------------------------\n";
    for (size_t i = 0; i < doses.size(); ++i) {
        double sf_with = sf_pred[i];
        double sf_without = calibrator.getModel().survivalFraction(doses[i], params_no_misrepair);
        double diff = sf_without - sf_with;
        std::cout << std::setw(10) << std::setprecision(2) << doses[i]
                  << std::setw(12) << std::setprecision(6) << sf_with
                  << std::setw(12) << sf_without
                  << std::setw(12) << diff << "\n";
    }
    std::cout << "\n";
    std::cout << "Note: Positive difference means misrepair increases cell killing\n\n";

    // Summary
    std::cout << "========================================\n";
    std::cout << "SUMMARY\n";
    std::cout << "========================================\n";
    std::cout << "This calibration uses the stochastic misrepair model:\n";
    std::cout << "  - 4 OPTIMIZED parameters: e2, e3, a, k_error\n";
    std::cout << "  - 2 USER FIXED parameters: T21, T23\n";
    std::cout << "  - Misrepair hazard: lambda_mis = k_error * N_dsb\n\n";
    std::cout << "The misrepair channel captures:\n";
    std::cout << "  - DSB-dependent probability of error-prone repair\n";
    std::cout << "  - Additional death pathway beyond energy threshold\n";
    std::cout << "  - More realistic dose-response at higher doses\n\n";
    std::cout << "Key values:\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  e2      = " << result.params.e2 << " (FITTED)\n";
    std::cout << "  e3      = " << result.params.e3 << " (FITTED)\n";
    std::cout << "  a       = " << result.params.a << " (FITTED)\n";
    std::cout << "  k_error = " << result.params.k_error << " per (DSB*hour) (FITTED)\n";
    std::cout << "  T21     = " << result.params.T21 << " hours (USER FIXED)\n";
    std::cout << "  T23     = " << result.params.T23 << " hours (USER FIXED)\n\n";
    std::cout << "Note: Fixing T21 and T23 prevents parameter inflation\n";
    std::cout << "      and produces more stable, identifiable fits.\n";
    std::cout << "========================================\n";
    std::cout << "Calibration test completed successfully!\n";
    std::cout << "========================================\n";

    return 0;
}
