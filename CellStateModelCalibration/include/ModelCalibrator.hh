#ifndef MODEL_CALIBRATOR_HH
#define MODEL_CALIBRATOR_HH

/**
 * @file ModelCalibrator.hh
 * @brief Main calibration class for cell state model
 * 
 * Provides high-level interface for:
 * - Fitting model parameters to experimental data
 * - Computing negative log-likelihood
 * - Uncertainty analysis via Hessian
 * - Profile likelihood computation
 */

#include "DataTypes.hh"
#include "SurvivalModel.hh"
#include "NelderMeadOptimizer.hh"
#include <vector>
#include <utility>

namespace CellStateCalibration {

/**
 * @struct FitResult
 * @brief Results from model calibration
 */
struct FitResult {
    ParamsReduced params;      ///< Optimal parameters
    double negLogLikelihood;   ///< Negative log-likelihood at optimum
    Hessian2 hessian;          ///< Hessian matrix at optimum
    int iterations;            ///< Number of optimizer iterations
    bool converged;            ///< Whether optimization converged
    
    FitResult() 
        : negLogLikelihood(std::numeric_limits<double>::infinity())
        , iterations(0)
        , converged(false) 
    {}
};

/**
 * @class ModelCalibrator
 * @brief High-level interface for model calibration
 */
class ModelCalibrator {
public:
    /**
     * @brief Constructor
     * @param config Fitting configuration
     */
    explicit ModelCalibrator(const FitConfig& config = FitConfig());
    
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
    FitResult fit(const ParamsReduced& initial_guess = ParamsReduced(20.0, 1.0),
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
     * 
     * Fixes e3 at grid points and optimizes over a to get
     * the profile likelihood curve.
     * 
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
     * @brief Get survival model
     * @return Reference to survival model
     */
    const SurvivalModel& getModel() const;

private:
    FitConfig m_config;
    std::vector<DataPoint> m_data;
    SurvivalModel m_model;
    mutable NelderMeadOptimizer m_optimizer;
};

} // namespace CellStateCalibration

#endif // MODEL_CALIBRATOR_HH
