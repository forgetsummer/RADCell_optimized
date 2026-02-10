#ifndef DIRECT_LETHAL_CALIBRATOR_HH
#define DIRECT_LETHAL_CALIBRATOR_HH

/**
 * @file DirectLethalCalibrator.hh
 * @brief Calibrator for Direct Lethal Model (S1->S3 only)
 * 
 * This calibrator fits the Direct Lethal Model to survival curve data.
 * It only considers the direct S1->S3 transition (no repair state S2).
 * 
 * Provides:
 * - Maximum likelihood parameter estimation
 * - Uncertainty analysis via Hessian
 * - Profile likelihood computation
 * 
 * For models including repair (S1->S2->S3), see RepairMediatedCalibrator.
 */

#include "DataTypes.hh"
#include "DirectLethalModel.hh"
#include "NelderMeadOptimizer.hh"
#include <vector>
#include <utility>

namespace CellStateCalibration {

/**
 * @struct DirectLethalFitResult
 * @brief Results from Direct Lethal Model calibration
 */
struct DirectLethalFitResult {
    ParamsReduced params;      ///< Optimal parameters (e3 = E3/σ, a = α/σ)
    double negLogLikelihood;   ///< Negative log-likelihood at optimum
    Hessian2 hessian;          ///< Hessian matrix at optimum
    int iterations;            ///< Number of optimizer iterations
    bool converged;            ///< Whether optimization converged
    
    DirectLethalFitResult() 
        : negLogLikelihood(std::numeric_limits<double>::infinity())
        , iterations(0)
        , converged(false) 
    {}
    
    /**
     * @brief Get physical parameters given a sigma value
     * @param sigma Chosen sigma value
     * @param kappa DSB/Gy value
     * @return Full physical parameters for CellStateModel
     */
    CellStateModelParams getPhysicalParams(double sigma, double kappa = 40.0) const {
        return CellStateModelParams::fromReduced(params, sigma, kappa);
    }
};

/**
 * @class DirectLethalCalibrator
 * @brief Calibrator for Direct Lethal Model (S1->S3 transition only)
 * 
 * Fits the model:
 *   SF(D) = 1 - Σ_n P(N=n | κD) * p13(n)
 * 
 * where:
 *   N ~ Poisson(κD) is DSB count
 *   p13(n) = 2Φ(-|an - e3|/2) is S1->S3 transition probability
 * 
 * Optimizes parameters:
 *   e3 = E3/σ (dimensionless death threshold)
 *   a = α/σ (dimensionless damage sensitivity)
 */
class DirectLethalCalibrator {
public:
    /**
     * @brief Constructor
     * @param config Fitting configuration
     */
    explicit DirectLethalCalibrator(const FitConfig& config = FitConfig());
    
    /**
     * @brief Set fitting configuration
     * @param config New configuration
     */
    void setConfig(const FitConfig& config);
    
    /**
     * @brief Set experimental data
     * @param data Vector of data points
     */
    void setData(const std::vector<DataPoint>& data);
    
    /**
     * @brief Add single data point
     * @param point Data point to add
     */
    void addDataPoint(const DataPoint& point);
    
    /**
     * @brief Clear all data
     */
    void clearData();
    
    /**
     * @brief Set default uncertainties for data points
     * 
     * If no error bars are provided, assumes constant relative
     * error which translates to constant std dev in log-space.
     * 
     * @param rel_err Relative error (e.g., 0.2 for 20%)
     */
    void setDefaultUncertainties(double rel_err = 0.2);
    
    /**
     * @brief Compute negative log-likelihood
     * @param params Model parameters
     * @return Negative log-likelihood value
     */
    double negativeLogLikelihood(const ParamsReduced& params) const;
    
    /**
     * @brief Fit model to data
     * @param initial_guess Initial parameter values
     * @param step Initial step sizes for optimizer
     * @return Fit results including optimal parameters and diagnostics
     */
    DirectLethalFitResult fit(const ParamsReduced& initial_guess = ParamsReduced(20.0, 1.0),
                               const ParamsReduced& step = ParamsReduced(5.0, 0.5));
    
    /**
     * @brief Compute Hessian at given parameters
     * @param params Parameters at which to evaluate Hessian
     * @param h_e3 Step size for e3 finite difference
     * @param h_a Step size for a finite difference
     * @return 2x2 Hessian matrix
     */
    Hessian2 computeHessian(const ParamsReduced& params,
                            double h_e3 = -1.0,  // -1 means auto
                            double h_a = -1.0) const;
    
    /**
     * @brief Compute profile likelihood for e3
     * @param params_opt Optimal parameters (center of grid)
     * @param e3_min Minimum e3 value
     * @param e3_max Maximum e3 value
     * @param n_grid Number of grid points
     * @return Vector of (e3, profile_nll) pairs
     */
    std::vector<std::pair<double, double>> profileLikelihood_e3(
        const ParamsReduced& params_opt,
        double e3_min,
        double e3_max,
        int n_grid = 41) const;
    
    /**
     * @brief Compute profile likelihood for a
     * @param params_opt Optimal parameters
     * @param a_min Minimum a value
     * @param a_max Maximum a value
     * @param n_grid Number of grid points
     * @return Vector of (a, profile_nll) pairs
     */
    std::vector<std::pair<double, double>> profileLikelihood_a(
        const ParamsReduced& params_opt,
        double a_min,
        double a_max,
        int n_grid = 41) const;
    
    /**
     * @brief Get predicted survival curve
     * @param params Model parameters
     * @param doses Dose values
     * @return Predicted survival fractions
     */
    std::vector<double> predictSurvival(const ParamsReduced& params,
                                         const std::vector<double>& doses) const;
    
    /**
     * @brief Get current data
     * @return Reference to data vector
     */
    const std::vector<DataPoint>& getData() const;
    
    /**
     * @brief Get the Direct Lethal Model
     * @return Reference to model
     */
    const DirectLethalModel& getModel() const;

private:
    FitConfig m_config;
    std::vector<DataPoint> m_data;
    DirectLethalModel m_model;
    mutable NelderMeadOptimizer m_optimizer;
};

} // namespace CellStateCalibration

#endif // DIRECT_LETHAL_CALIBRATOR_HH
