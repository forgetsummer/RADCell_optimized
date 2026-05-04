#ifndef REPAIR_MEDIATED_CALIBRATOR_MULTICOMPONENT_HH
#define REPAIR_MEDIATED_CALIBRATOR_MULTICOMPONENT_HH

/**
 * @file RepairMediatedCalibrator_multicomponent.hh
 * @brief Calibrator for the multi-component energy model (Er + Ep)
 *
 * Jointly optimizes up to 7 parameters from survival data:
 *   e2, e3, a, k_error, Nc, omega_p, lambda_p
 * (with T21/T23 fixed or optimized via FitConfigS2 flags).
 *
 * This eliminates the need for a separate grid scan over (Nc, omega_p, lambda_p).
 */

#include "DataTypes.hh"
#include "RepairMediatedModel_multicomponent.hh"
#include <vector>
#include <functional>

namespace CellStateCalibration {

struct MultiComponentFitResult {
    ParamsS2 params;
    double negLogLikelihood;
    int iterations;
    bool converged;

    MultiComponentFitResult()
        : negLogLikelihood(std::numeric_limits<double>::infinity())
        , iterations(0)
        , converged(false)
    {}
};

class RepairMediatedMultiComponentCalibrator {
public:
    explicit RepairMediatedMultiComponentCalibrator(const FitConfigS2& config = FitConfigS2());

    void setConfig(const FitConfigS2& config);
    void setData(const std::vector<DataPoint>& data);
    void addDataPoint(const DataPoint& point);
    void clearData();
    void setDefaultUncertainties(double rel_err = 0.2);

    double negativeLogLikelihood(const ParamsS2& params) const;

    /**
     * @brief Fit model to data using Nelder-Mead optimization
     *
     * With T21/T23 fixed (default), optimizes 7 parameters:
     *   e2, e3, a, k_error, Nc, omega_p, lambda_p
     *
     * With T21/T23 free, optimizes up to 9 parameters.
     */
    MultiComponentFitResult fit(
        const ParamsS2& initial_guess,
        const ParamsS2& step);

    std::vector<double> predictSurvival(const ParamsS2& params,
                                         const std::vector<double>& doses) const;

    const std::vector<DataPoint>& getData() const;
    const RepairMediatedMultiComponentModel& getModel() const;

private:
    FitConfigS2 m_config;
    std::vector<DataPoint> m_data;
    RepairMediatedMultiComponentModel m_model;

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

#endif // REPAIR_MEDIATED_CALIBRATOR_MULTICOMPONENT_HH
