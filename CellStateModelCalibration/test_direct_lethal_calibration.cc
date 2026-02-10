/**
 * @file test_calibration.cc
 * @brief Test program for CellStateModelCalibration module
 * 
 * Demonstrates:
 * 1. Setting up experimental data
 * 2. Fitting model parameters
 * 3. Uncertainty analysis via Hessian
 * 4. Profile likelihood computation
 * 5. Comparing model predictions with data
 */

#include "DirectLethalCalibrator.hh"
#include "DirectLethalModel.hh"
#include "DataTypes.hh"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

using namespace CellStateCalibration;

int main() {
    std::cout << "========================================\n";
    std::cout << "Direct Lethal Model Calibration Test\n";
    std::cout << "(S1 -> S3 transition only)\n";
    std::cout << "========================================\n\n";
    
    //------------------------------------------------------------------
    // 1. Setup configuration
    //------------------------------------------------------------------
    FitConfig config;
    config.kappa = 40.0;         // DSB/Gy
    config.dose_max_fit = 4.0;   // Fit data up to 4 Gy
    config.eps_sf = 1e-12;       // Floor for log calculations
    
    std::cout << "Configuration:\n";
    std::cout << "  kappa (DSB/Gy): " << config.kappa << "\n";
    std::cout << "  dose_max_fit: " << config.dose_max_fit << " Gy\n\n";
    
    //------------------------------------------------------------------
    // 2. Create calibrator and add experimental data
    //------------------------------------------------------------------
    DirectLethalCalibrator calibrator(config);
    
    // Example data (replace with your digitized survival curve)
    // Format: dose (Gy), survival fraction, uncertainty (log-space)
    std::vector<DataPoint> data = {
        {0.0, 1.0, 0.0},
        {1.0, 0.70, 0.0},
        {2.0, 0.45, 0.0},
        {3.0, 0.25, 0.0},
        {4.0, 0.12, 0.0},
        {6.0, 0.02, 0.0},   // Will be excluded due to dose_max_fit
    };
    
    calibrator.setData(data);
    calibrator.setDefaultUncertainties(0.2);  // 20% relative error
    
    std::cout << "Experimental data:\n";
    std::cout << std::setw(10) << "Dose (Gy)" << std::setw(15) << "SF_obs\n";
    std::cout << "  -------------------------\n";
    for (const auto& dp : calibrator.getData()) {
        std::cout << std::setw(10) << dp.D << std::setw(15) << dp.SF_obs << "\n";
    }
    std::cout << "\n";
    
    //------------------------------------------------------------------
    // 3. Fit model parameters
    //------------------------------------------------------------------
    std::cout << "Fitting model parameters...\n\n";
    
    ParamsReduced initial_guess(20.0, 1.0);  // e3, a
    ParamsReduced step(5.0, 0.5);
    
    DirectLethalFitResult result = calibrator.fit(initial_guess, step);
    
    std::cout << "=== Fit Results ===\n\n";
    
    // Reduced (identifiable) parameters
    std::cout << "Reduced Parameters (identifiable from data):\n";
    std::cout << "  e3 = E3/sigma = " << std::fixed << std::setprecision(4) << result.params.e3 << "\n";
    std::cout << "  a  = alpha/sigma = " << std::fixed << std::setprecision(4) << result.params.a << "\n";
    std::cout << "  Negative log-likelihood: " << std::setprecision(6) << result.negLogLikelihood << "\n";
    std::cout << "  Iterations: " << result.iterations << "\n";
    std::cout << "  Converged: " << (result.converged ? "Yes" : "No") << "\n\n";
    
    // Physical parameters (require choice of sigma)
    // The model has 4 parameters: E1, E3, alpha, sigma
    // But only 2 combinations are identifiable: e3 = E3/sigma, a = alpha/sigma
    // Convention: E1 = 0 (reference point)
    // To get physical parameters, we need to choose sigma
    
    std::cout << "Physical Parameters (for Cell State Model):\n";
    std::cout << "  Note: E1, E3, alpha, sigma are not individually identifiable.\n";
    std::cout << "  Only ratios e3=E3/sigma and a=alpha/sigma can be determined.\n";
    std::cout << "  Convention: E1 = 0 (healthy state reference)\n\n";
    
    // Show physical parameters for different choices of sigma
    std::cout << "  Physical parameters for different sigma choices:\n";
    std::cout << "  " << std::setw(10) << "sigma" 
              << std::setw(12) << "E3" 
              << std::setw(12) << "alpha" 
              << std::setw(15) << "E3 (meaning)\n";
    std::cout << "  " << std::string(49, '-') << "\n";
    
    std::vector<double> sigma_choices = {1.0, 5.0, 10.0, 20.0};
    for (double sigma : sigma_choices) {
        double E3 = result.params.e3 * sigma;
        double alpha = result.params.a * sigma;
        std::cout << "  " << std::setw(10) << std::fixed << std::setprecision(2) << sigma
                  << std::setw(12) << std::setprecision(4) << E3
                  << std::setw(12) << std::setprecision(4) << alpha
                  << std::setw(15) << "death threshold\n";
    }
    std::cout << "\n";
    
    // Interpretation
    std::cout << "Parameter Interpretation:\n";
    std::cout << "  - kappa = " << config.kappa << " DSB/Gy (input, from literature)\n";
    std::cout << "  - E1 = 0 (healthy state, by convention)\n";
    std::cout << "  - E3 = e3 * sigma = death threshold (damage level causing cell death)\n";
    std::cout << "  - alpha = a * sigma = damage per DSB (how much each DSB contributes)\n";
    std::cout << "  - sigma = Gaussian width (uncertainty in state transitions)\n\n";
    
    std::cout << "Key derived quantities:\n";
    double n_threshold = result.params.criticalDSBCount();
    std::cout << "  - Critical DSB count (E3/alpha = e3/a): " << std::setprecision(2) << n_threshold << " DSBs\n";
    std::cout << "    (This is the number of DSBs needed to reach death threshold on average)\n";
    double dose_at_threshold = n_threshold / config.kappa;
    std::cout << "  - Corresponding dose: " << std::setprecision(3) << dose_at_threshold << " Gy\n\n";
    
    //------------------------------------------------------------------
    // Get physical parameters for your Cell State Model
    //------------------------------------------------------------------
    std::cout << "========================================\n";
    std::cout << "CALIBRATED PARAMETERS FOR CELL STATE MODEL\n";
    std::cout << "========================================\n\n";
    
    std::cout << "To use these in your CellStateModel, you need to choose sigma.\n";
    std::cout << "Common approaches:\n";
    std::cout << "  1. Use sigma from literature/prior knowledge\n";
    std::cout << "  2. Set sigma = 1 (simplest, parameters become e3 and a directly)\n";
    std::cout << "  3. Choose sigma based on physical reasoning\n\n";
    
    // Example: Get parameters with sigma = 1.0 (simplest case)
    std::cout << "--- Example with sigma = 1.0 ---\n";
    CellStateModelParams params_s1 = CellStateModelParams::fromReduced(
        result.params, 1.0, config.kappa);
    params_s1.print();
    
    // Example: Get parameters with sigma = 10.0
    std::cout << "\n--- Example with sigma = 10.0 ---\n";
    CellStateModelParams params_s10 = CellStateModelParams::fromReduced(
        result.params, 10.0, config.kappa);
    params_s10.print();
    
    //------------------------------------------------------------------
    // FINAL: Choose sigma and get calibrated parameters
    //------------------------------------------------------------------
    std::cout << "\n========================================\n";
    std::cout << "CALIBRATED PARAMETERS FOR CELL STATE MODEL\n";
    std::cout << "========================================\n\n";
    
    std::cout << "Step 1: Choose your sigma value\n";
    std::cout << "Step 2: Calculate E3 = e3 * sigma = " << result.params.e3 << " * sigma\n";
    std::cout << "Step 3: Calculate alpha = a * sigma = " << result.params.a << " * sigma\n";
    std::cout << "Step 4: E1 = 0 (by convention)\n\n";
    
    // Let user choose sigma - here we demonstrate with sigma = 1.0
    double chosen_sigma = 1.0;  // <-- USER CHOOSES THIS VALUE
    
    std::cout << ">>> Using sigma = " << chosen_sigma << " <<<\n\n";
    
    CellStateModelParams final_params = CellStateModelParams::fromReduced(
        result.params, chosen_sigma, config.kappa);
    
    std::cout << "--- CALIBRATED PARAMETERS ---\n";
    final_params.print();
    
    std::cout << "\n--- Copy these values to your Cell setup ---\n";
    std::cout << std::fixed;
    std::cout << "E1    = " << std::setprecision(6) << final_params.E1 << "\n";
    std::cout << "E3    = " << std::setprecision(6) << final_params.E3 << "\n";
    std::cout << "sigma = " << std::setprecision(6) << final_params.sigma << "\n";
    std::cout << "alpha = " << std::setprecision(6) << final_params.alpha << "\n";
    std::cout << "kappa = " << std::setprecision(6) << final_params.kappa << " DSB/Gy\n\n";
    
    // Show formula for other sigma choices
    std::cout << "--- For other sigma values, use these formulas ---\n";
    std::cout << "E1    = 0\n";
    std::cout << "E3    = " << result.params.e3 << " * sigma\n";
    std::cout << "alpha = " << result.params.a << " * sigma\n";
    std::cout << "kappa = " << config.kappa << " DSB/Gy (fixed)\n\n";
    
    // Show how to use in code
    std::cout << "--- How to use in your code ---\n";
    std::cout << "double sigma = " << chosen_sigma << ";  // Choose your sigma value\n";
    std::cout << "CellStateModelParams params = CellStateModelParams::fromReduced(\n";
    std::cout << "    result.params, sigma, " << config.kappa << ");\n\n";
    
    //------------------------------------------------------------------
    // 4. Uncertainty analysis via Hessian
    //------------------------------------------------------------------
    std::cout << "=== Hessian Matrix at Optimum ===\n";
    std::cout << "  H11 (d²NLL/de3²): " << std::scientific << result.hessian.h11 << "\n";
    std::cout << "  H12 (d²NLL/de3da): " << result.hessian.h12 << "\n";
    std::cout << "  H22 (d²NLL/da²): " << result.hessian.h22 << "\n";
    std::cout << "  Determinant: " << result.hessian.determinant() << "\n";
    std::cout << "  Positive definite: " << (result.hessian.isPositiveDefinite() ? "Yes" : "No") << "\n\n";
    
    // Approximate standard errors from inverse Hessian
    if (result.hessian.isPositiveDefinite()) {
        double det = result.hessian.determinant();
        double se_e3 = std::sqrt(result.hessian.h22 / det);
        double se_a = std::sqrt(result.hessian.h11 / det);
        std::cout << "Approximate standard errors:\n";
        std::cout << "  SE(e3) ≈ " << std::fixed << std::setprecision(4) << se_e3 << "\n";
        std::cout << "  SE(a)  ≈ " << se_a << "\n\n";
    }
    
    //------------------------------------------------------------------
    // 5. Profile likelihood
    //------------------------------------------------------------------
    std::cout << "=== Profile Likelihood for e3 ===\n";
    double e3_min = std::max(0.1, result.params.e3 * 0.5);
    double e3_max = result.params.e3 * 1.5;
    
    auto profile_e3 = calibrator.profileLikelihood_e3(result.params, e3_min, e3_max, 21);
    
    std::cout << std::setw(12) << "e3" << std::setw(15) << "Profile NLL\n";
    std::cout << "  ---------------------------\n";
    for (const auto& pt : profile_e3) {
        std::cout << std::setw(12) << std::fixed << std::setprecision(3) << pt.first 
                  << std::setw(15) << std::setprecision(6) << pt.second << "\n";
    }
    std::cout << "\n";
    
    //------------------------------------------------------------------
    // 6. Model predictions vs data
    //------------------------------------------------------------------
    std::cout << "=== Model Predictions vs Data ===\n";
    std::cout << std::setw(10) << "Dose (Gy)" 
              << std::setw(12) << "SF_obs" 
              << std::setw(12) << "SF_model"
              << std::setw(12) << "Residual\n";
    std::cout << "  ------------------------------------------------\n";
    
    for (const auto& dp : calibrator.getData()) {
        double sf_model = calibrator.getModel().survivalFraction(dp.D, result.params);
        double residual = std::log(dp.SF_obs) - std::log(sf_model);
        
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << dp.D
                  << std::setw(12) << std::setprecision(4) << dp.SF_obs
                  << std::setw(12) << sf_model
                  << std::setw(12) << std::setprecision(4) << residual << "\n";
    }
    std::cout << "\n";
    
    //------------------------------------------------------------------
    // 7. Generate smooth survival curve
    //------------------------------------------------------------------
    std::cout << "=== Predicted Survival Curve ===\n";
    std::vector<double> doses;
    for (double d = 0.0; d <= 8.0; d += 0.5) {
        doses.push_back(d);
    }
    
    auto sf_curve = calibrator.predictSurvival(result.params, doses);
    
    std::cout << std::setw(10) << "Dose (Gy)" << std::setw(15) << "SF_predicted\n";
    std::cout << "  -------------------------\n";
    for (size_t i = 0; i < doses.size(); ++i) {
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << doses[i]
                  << std::setw(15) << std::setprecision(6) << sf_curve[i] << "\n";
    }
    
    //------------------------------------------------------------------
    // 8. Save results to CSV
    //------------------------------------------------------------------
    std::ofstream csv_file("calibration_results.csv");
    csv_file << "Dose_Gy,SF_observed,SF_model\n";
    for (const auto& dp : calibrator.getData()) {
        double sf_model = calibrator.getModel().survivalFraction(dp.D, result.params);
        csv_file << dp.D << "," << dp.SF_obs << "," << sf_model << "\n";
    }
    csv_file.close();
    
    std::ofstream curve_file("survival_curve_fit.csv");
    curve_file << "Dose_Gy,SF_model\n";
    for (size_t i = 0; i < doses.size(); ++i) {
        curve_file << doses[i] << "," << sf_curve[i] << "\n";
    }
    curve_file.close();
    
    std::cout << "\nResults saved to:\n";
    std::cout << "  - calibration_results.csv\n";
    std::cout << "  - survival_curve_fit.csv\n";
    
    std::cout << "\n========================================\n";
    std::cout << "Calibration test completed successfully!\n";
    std::cout << "========================================\n";
    
    return 0;
}
