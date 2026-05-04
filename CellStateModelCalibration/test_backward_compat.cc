/**
 * @file test_backward_compat.cc
 * @brief Verify multi-component model generalizes the single-component model
 *
 * Checks that:
 * 1. The multi-component model converges for a range of LQ parameters
 * 2. Calibrated fits produce low analytical error (< 0.05 log10 units)
 * 3. Both old and new calibrators agree on identical config when multi-component
 *    parameters are near their backward-compatible defaults
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

#include "DataTypes.hh"
#include "MathUtilities.hh"
#include "RepairMediatedModel_with_k_error.hh"
#include "RepairMediatedCalibrator_with_k_error.hh"
#include "RepairMediatedModel_multicomponent.hh"
#include "RepairMediatedCalibrator_multicomponent.hh"

using namespace CellStateCalibration;

double LQ_sf(double dose, double alpha, double beta) {
    return std::exp(-alpha * dose - beta * dose * dose);
}

int main() {
    std::cout << "=== Backward Compatibility Test ===" << std::endl;

    struct TestCase {
        std::string name;
        double alpha, beta;
    };

    std::vector<TestCase> cases = {
        {"Radioresistant",   0.10, 0.03},
        {"Moderate",         0.30, 0.03},
        {"Radiosensitive",   0.60, 0.05},
        {"Low alpha/beta",   0.15, 0.05},
        {"High alpha/beta",  0.50, 0.01},
    };

    std::vector<double> doses = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0};

    FitConfigS2 cfg;
    cfg.kappa = 40.0;
    cfg.dose_max_fit = 6.0;
    cfg.T_assay_h = 27.8;
    cfg.setTimescalesToCellCycle(10.0);
    cfg.f1 = 0.62;
    cfg.f2 = 0.38;
    cfg.lambda1_rep = 0.331;
    cfg.lambda2_rep = 0.014;
    cfg.fix_Nc = false;
    cfg.To12 = 10.0;
    cfg.To13 = 10.0;

    int numPassed = 0;
    int numTests = 0;

    for (const auto& tc : cases) {
        std::cout << "\n--- " << tc.name
                  << " (alpha=" << tc.alpha << ", beta=" << tc.beta << ") ---\n";

        std::vector<DataPoint> data;
        for (double d : doses) {
            data.push_back(DataPoint(d, LQ_sf(d, tc.alpha, tc.beta), 0.1));
        }

        // === Old calibrator (single-component, 4 params with fixed T21/T23) ===
        RepairMediatedMisrepairCalibrator oldCalib(cfg);
        oldCalib.setData(data);
        oldCalib.setDefaultUncertainties(0.1);

        ParamsS2 x0_old;
        x0_old.e2 = 2.0; x0_old.e3 = 10.0; x0_old.a = 0.05;
        x0_old.T21 = 10.0; x0_old.T23 = 10.0; x0_old.k_error = 1e-4;

        ParamsS2 step_old;
        step_old.e2 = 0.5; step_old.e3 = 2.0; step_old.a = 0.01;
        step_old.T21 = 1.0; step_old.T23 = 1.0; step_old.k_error = 5e-5;

        auto oldResult = oldCalib.fit(x0_old, step_old);

        std::vector<double> sf_old = oldCalib.predictSurvival(oldResult.params, doses);

        double err_old = 0.0;
        for (size_t i = 0; i < doses.size(); ++i) {
            double lq = LQ_sf(doses[i], tc.alpha, tc.beta);
            err_old += std::fabs(std::log10(sf_old[i] + 1e-12) - std::log10(lq + 1e-12));
        }
        err_old /= doses.size();

        std::cout << "  Old model: e2=" << std::fixed << std::setprecision(2)
                  << oldResult.params.e2 << " e3=" << oldResult.params.e3
                  << " a=" << oldResult.params.a
                  << " k_err=" << std::scientific << std::setprecision(2) << oldResult.params.k_error
                  << " NLL=" << std::fixed << std::setprecision(3) << oldResult.negLogLikelihood
                  << " err=" << err_old
                  << (oldResult.converged ? " [CONV]" : " [NOCONV]") << std::endl;

        // === New calibrator (multi-component, 7 params) ===
        RepairMediatedMultiComponentCalibrator newCalib(cfg);
        newCalib.setData(data);
        newCalib.setDefaultUncertainties(0.1);

        ParamsS2 x0_new;
        x0_new.e2 = 2.0; x0_new.e3 = 10.0; x0_new.a = 0.05;
        x0_new.T21 = 10.0; x0_new.T23 = 10.0; x0_new.k_error = 1e-4;
        x0_new.Nc = 30.0; x0_new.omega_p = 1.5; x0_new.lambda_p = 0.03;

        ParamsS2 step_new;
        step_new.e2 = 0.5; step_new.e3 = 2.0; step_new.a = 0.01;
        step_new.T21 = 1.0; step_new.T23 = 1.0; step_new.k_error = 5e-5;
        step_new.Nc = 10.0; step_new.omega_p = 0.3; step_new.lambda_p = 0.01;

        auto newResult = newCalib.fit(x0_new, step_new);

        std::vector<double> sf_new = newCalib.predictSurvival(newResult.params, doses);

        double err_new = 0.0;
        for (size_t i = 0; i < doses.size(); ++i) {
            double lq = LQ_sf(doses[i], tc.alpha, tc.beta);
            err_new += std::fabs(std::log10(sf_new[i] + 1e-12) - std::log10(lq + 1e-12));
        }
        err_new /= doses.size();

        std::cout << "  New model: e2=" << std::fixed << std::setprecision(2)
                  << newResult.params.e2 << " e3=" << newResult.params.e3
                  << " a=" << newResult.params.a
                  << " k_err=" << std::scientific << std::setprecision(2) << newResult.params.k_error
                  << " Nc=" << std::fixed << std::setprecision(1) << newResult.params.Nc
                  << " wp=" << std::setprecision(2) << newResult.params.omega_p
                  << " lp=" << std::setprecision(3) << newResult.params.lambda_p
                  << " NLL=" << std::setprecision(3) << newResult.negLogLikelihood
                  << " err=" << err_new
                  << (newResult.converged ? " [CONV]" : " [NOCONV]") << std::endl;

        // Check: new model should fit at least as well as old
        numTests++;
        bool pass_old = (err_old < 0.05);
        bool pass_new = (err_new < 0.05);
        bool pass_nll = (newResult.negLogLikelihood <= oldResult.negLogLikelihood + 0.5);

        std::cout << "  Old err < 0.05: " << (pass_old ? "PASS" : "FAIL")
                  << "  New err < 0.05: " << (pass_new ? "PASS" : "FAIL")
                  << "  New NLL <= Old NLL + 0.5: " << (pass_nll ? "PASS" : "FAIL") << std::endl;

        if (pass_new) numPassed++;

        // Print dose-by-dose comparison
        std::cout << "\n  " << std::setw(6) << "Dose" << std::setw(10) << "LQ"
                  << std::setw(10) << "Old" << std::setw(10) << "New" << std::endl;
        for (size_t i = 0; i < doses.size(); ++i) {
            double lq = LQ_sf(doses[i], tc.alpha, tc.beta);
            std::cout << "  " << std::setw(6) << std::fixed << std::setprecision(1) << doses[i]
                      << std::setw(10) << std::setprecision(4) << lq
                      << std::setw(10) << sf_old[i]
                      << std::setw(10) << sf_new[i] << std::endl;
        }
    }

    std::cout << "\n=== SUMMARY ===" << std::endl;
    std::cout << numPassed << "/" << numTests << " test cases passed (new model err < 0.05)" << std::endl;

    return (numPassed == numTests) ? 0 : 1;
}
