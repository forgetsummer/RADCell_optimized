#ifndef REPAIR_MEDIATED_MODEL_MULTICOMPONENT_HH
#define REPAIR_MEDIATED_MODEL_MULTICOMPONENT_HH

/**
 * @file RepairMediatedModel_multicomponent.hh
 * @brief Analytical survival model with multi-component energy (Er + Ep)
 *
 * Extends the single-component misrepair model to include:
 * - Repairable energy Er(t): decays via fast bi-exponential (f1*exp(-l1*t) + f2*exp(-l2*t))
 * - Persistent energy Ep(t): decays via slow mono-exponential exp(-lambda_p * t)
 * - Effective energy for S2 transitions: E_eff(t) = Er(t) + omega_p * Ep(t)
 *
 * The persistent fraction of initial damage depends on DSB count:
 *   h(N) = N / (N + Nc)
 *
 * This makes S2 transition rates TIME-DEPENDENT, requiring numerical
 * integration over [0, T_assay] instead of the closed-form competing
 * exponential used in the single-component model.
 *
 * Backward compatible: when Nc >> N (default 1e9), h(N) ~ 0, Ep ~ 0,
 * and the model reduces to the single-component case.
 */

#include "DataTypes.hh"
#include "MathUtilities.hh"
#include <vector>

namespace CellStateCalibration {

class RepairMediatedMultiComponentModel {
public:
    explicit RepairMediatedMultiComponentModel(const FitConfigS2& config = FitConfigS2());

    void setConfig(const FitConfigS2& config);
    const FitConfigS2& getConfig() const;

    double overlapScore(double e, double ej) const;

    void assignmentProbabilities(int n, const ParamsS2& params,
                                  double& P1, double& P2, double& P3) const;

    double cumulativeProbToRate(double p_cum, double T0) const;

    /**
     * @brief Compute survival fraction with time-dependent multi-component energy
     *
     * SF(D) = Sum_n P(N=n|kappa*D) * [P1(n) + P2(n) * R21(n)]
     *
     * R21(n) is computed by numerically integrating competing time-dependent
     * rates lambda21(t) and lambda23(t) over [0, T_assay].
     */
    double survivalFraction(double dose, const ParamsS2& params) const;

    std::vector<double> survivalCurve(const std::vector<double>& doses,
                                       const ParamsS2& params) const;

private:
    FitConfigS2 m_config;

    /**
     * @brief Compute Er(t) for acute irradiation with n DSBs
     *
     * Er(t) = alpha*n*(1-h) * [f1*exp(-l1*t) + f2*exp(-l2*t)]
     */
    double Er(double t, int n, double alpha_n_1mh) const;

    /**
     * @brief Compute Ep(t) for acute irradiation with n DSBs
     *
     * Ep(t) = alpha*n*h * exp(-lambda_p*t)
     */
    double Ep(double t, double alpha_n_h, double lambda_p) const;

    /**
     * @brief Compute repair probability R21 for an S2 cell with n DSBs
     *
     * Numerically integrates:
     *   R21 = integral_0^T lambda21(t) * S2(t) dt
     * where S2(t) = exp(-integral_0^t [lambda21(s) + lambda23(s)] ds)
     * and the rates are derived from time-dependent E_eff(t).
     */
    double computeS2RepairProbability(int n, const ParamsS2& params) const;

    static const int NUM_INTEGRATION_STEPS = 200;
};

} // namespace CellStateCalibration

#endif // REPAIR_MEDIATED_MODEL_MULTICOMPONENT_HH
