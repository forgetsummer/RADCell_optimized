#ifndef REPAIR_MEDIATED_CALIBRATOR_WITH_K_ERROR_HH
#define REPAIR_MEDIATED_CALIBRATOR_WITH_K_ERROR_HH

/**
 * @file RepairMediatedCalibrator_with_k_error.hh
 * @brief Calibrator for Repair Mediated model with stochastic misrepair (S1->S2->S3)
 *
 * Fits the 6-parameter Repair Mediated + misrepair model to experimental survival data
 * using maximum likelihood estimation with Nelder-Mead optimization.
 *
 * ALL 6 parameters are OPTIMIZED from survival data:
 * - e2: Repair state threshold (E2/sigma)
 * - e3: Death state threshold (E3/sigma)
 * - a: Damage per DSB (alpha/sigma)
 * - T21: Repair timescale (hours) - S2->S1
 * - T23: Death timescale (hours) - S2->S3
 * - k_error: Misrepair rate per (DSB*hour) contributing to S2->S3
 */

#include "DataTypes.hh"
#include "RepairMediatedModel_with_k_error.hh"
#include <vector>
#include <functional>

namespace CellStateCalibration {

/**
 * @brief Result of Repair Mediated model fitting with k_error
 */
struct RepairMediatedMisrepairFitResult {
    ParamsS2 params;              ///< Fitted parameters (including k_error)
    double negLogLikelihood;      ///< Final NLL value
    int iterations;               ///< Number of iterations
    bool converged;               ///< Whether optimization converged

    RepairMediatedMisrepairFitResult()
        : negLogLikelihood(std::numeric_limits<double>::infinity())
        , iterations(0)
        , converged(false)
    {}
};

/**
 * @brief Calibrator for Repair Mediated survival model with stochastic misrepair
 *
 * This version optimizes 6 parameters including k_error for the misrepair channel.
 */
class RepairMediatedMisrepairCalibrator {
public:
    /**
     * @brief Constructor with optional configuration
     */
    explicit RepairMediatedMisrepairCalibrator(const FitConfigS2& config = FitConfigS2());

    /**
     * @brief Set configuration parameters
     */
    void setConfig(const FitConfigS2& config);

    /**
     * @brief Set experimental data
     */
    void setData(const std::vector<DataPoint>& data);

    /**
     * @brief Add single data point
     */
    void addDataPoint(const DataPoint& point);

    /**
     * @brief Clear all data
     */
    void clearData();

    /**
     * @brief Set default uncertainties for data points without errors
     * @param rel_err Relative error (default 0.2 = 20%)
     */
    void setDefaultUncertainties(double rel_err = 0.2);

    /**
     * @brief Compute negative log-likelihood for given parameters
     */
    double negativeLogLikelihood(const ParamsS2& params) const;

    /**
     * @brief Fit model to data using Nelder-Mead optimization
     *
     * Fits ALL 6 parameters: e2, e3, a, T21, T23, k_error
     *
     * @param initial_guess Starting point for optimization
     * @param step Initial step sizes for each parameter
     * @return Fit result including parameters and diagnostics
     */
    RepairMediatedMisrepairFitResult fit(
        const ParamsS2& initial_guess = ParamsS2(10.0, 25.0, 1.0, 12.0, 48.0, 0.001),
        const ParamsS2& step = ParamsS2(2.0, 5.0, 0.3, 4.0, 12.0, 0.0005));

    /**
     * @brief Predict survival for given parameters and doses
     */
    std::vector<double> predictSurvival(const ParamsS2& params,
                                         const std::vector<double>& doses) const;

    /**
     * @brief Get the experimental data
     */
    const std::vector<DataPoint>& getData() const;

    /**
     * @brief Get the model
     */
    const RepairMediatedMisrepairModel& getModel() const;

private:
    FitConfigS2 m_config;
    std::vector<DataPoint> m_data;
    RepairMediatedMisrepairModel m_model;

    /**
     * @brief N-dimensional Nelder-Mead optimizer
     */
    std::vector<double> nelderMeadND(
        std::function<double(const std::vector<double>&)> f,
        const std::vector<double>& x0,
        const std::vector<double>& step,
        const std::vector<double>& lo,
        const std::vector<double>& hi,
        double& best_f_out,
        int& iterations_out,
        bool& converged_out) const;
};

} // namespace CellStateCalibration

#endif // REPAIR_MEDIATED_CALIBRATOR_WITH_K_ERROR_HH
