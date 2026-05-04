#ifndef SIMULATION_BASED_CALIBRATOR_HH
#define SIMULATION_BASED_CALIBRATOR_HH

/**
 * @file SimulationBasedCalibrator.hh
 * @brief Calibrator that uses an external simulation function instead of an analytical model
 *
 * Accepts a std::function<double(double dose, const ParamsS2& params)> that returns
 * the simulated survival fraction for a given dose and parameters. This allows
 * calibrating model parameters directly against the simulation's own output,
 * bypassing the analytical model entirely.
 *
 * Optimizes all 6 parameters: e2, e3, a, k_error, T21, T23
 * using Nelder-Mead with the same NLL objective as other calibrators.
 */

#include "DataTypes.hh"
#include <vector>
#include <functional>
#include <limits>

namespace CellStateCalibration {

/**
 * @brief Result of simulation-based fitting
 */
struct SimulationBasedFitResult {
    ParamsS2 params;
    double negLogLikelihood;
    int iterations;
    bool converged;

    SimulationBasedFitResult()
        : negLogLikelihood(std::numeric_limits<double>::infinity())
        , iterations(0)
        , converged(false)
    {}
};

/**
 * @brief Calibrator that optimizes parameters against a user-supplied simulation function
 *
 * The survival function is injected via setSurvivalFunction(), keeping this class
 * decoupled from any specific simulation implementation (e.g. CellStateModel).
 */
class SimulationBasedCalibrator {
public:
    using SurvivalFunction = std::function<double(double dose, const ParamsS2& params)>;
    using BatchSurvivalFunction = std::function<std::vector<double>(
        const std::vector<double>& doses, const ParamsS2& params)>;

    explicit SimulationBasedCalibrator(const FitConfigS2& config = FitConfigS2());

    void setConfig(const FitConfigS2& config);
    void setSurvivalFunction(SurvivalFunction fn);
    void setBatchSurvivalFunction(BatchSurvivalFunction fn);
    void setData(const std::vector<DataPoint>& data);
    void addDataPoint(const DataPoint& point);
    void clearData();
    void setDefaultUncertainties(double rel_err = 0.2);

    double negativeLogLikelihood(const ParamsS2& params) const;

    /**
     * @brief Fit by running simulation-in-the-loop Nelder-Mead
     * @param initial_guess Starting parameters (typically from analytical calibration)
     * @param step Initial simplex step sizes
     * @param max_iter Maximum Nelder-Mead iterations
     * @return Fit result with optimized parameters
     */
    SimulationBasedFitResult fit(
        const ParamsS2& initial_guess = ParamsS2(10.0, 25.0, 1.0, 12.0, 48.0, 0.001),
        const ParamsS2& step = ParamsS2(2.0, 5.0, 0.3, 4.0, 12.0, 0.0005),
        int max_iter = 500);

    const std::vector<DataPoint>& getData() const;

private:
    FitConfigS2 m_config;
    std::vector<DataPoint> m_data;
    SurvivalFunction m_survivalFn;
    BatchSurvivalFunction m_batchSurvivalFn;

    std::vector<double> nelderMeadND(
        std::function<double(const std::vector<double>&)> f,
        const std::vector<double>& x0,
        const std::vector<double>& step,
        const std::vector<double>& lo,
        const std::vector<double>& hi,
        double& best_f_out,
        int& iterations_out,
        bool& converged_out,
        int max_iter) const;
};

} // namespace CellStateCalibration

#endif // SIMULATION_BASED_CALIBRATOR_HH
