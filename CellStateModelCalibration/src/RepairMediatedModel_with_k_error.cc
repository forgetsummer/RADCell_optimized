/**
 * @file RepairMediatedModel_with_k_error.cc
 * @brief Implementation of Repair Mediated survival model with stochastic misrepair
 */

#include "RepairMediatedModel_with_k_error.hh"
#include <cmath>
#include <algorithm>

namespace CellStateCalibration {

RepairMediatedMisrepairModel::RepairMediatedMisrepairModel(const FitConfigS2& config)
    : m_config(config)
{}

void RepairMediatedMisrepairModel::setConfig(const FitConfigS2& config) {
    m_config = config;
}

const FitConfigS2& RepairMediatedMisrepairModel::getConfig() const {
    return m_config;
}

double RepairMediatedMisrepairModel::overlapScore(double e, double ej) const {
    // w(e, ej) = 2 * Phi(-|e - ej| / 2)
    double z = -std::abs(e - ej) / 2.0;
    double val = 2.0 * normal_cdf(z);
    return clamp(val, 0.0, 1.0);
}

void RepairMediatedMisrepairModel::assignmentProbabilities(int n, const ParamsS2& params,
                                                   double& P1, double& P2, double& P3) const {
    // Reduced energy after n DSBs
    const double e = params.a * static_cast<double>(n);

    // Overlap scores with each state
    const double w1 = overlapScore(e, 0.0);       // S1 at E1=0
    const double w2 = overlapScore(e, params.e2); // S2 at E2
    const double w3 = overlapScore(e, params.e3); // S3 at E3

    const double sum = w1 + w2 + w3;

    if (sum <= 1e-15) {
        // Extremely separated states - assign to closest
        double d1 = std::abs(e - 0.0);
        double d2 = std::abs(e - params.e2);
        double d3 = std::abs(e - params.e3);

        P1 = (d1 <= d2 && d1 <= d3) ? 1.0 : 0.0;
        P2 = (d2 < d1 && d2 <= d3) ? 1.0 : 0.0;
        P3 = (d3 < d1 && d3 < d2) ? 1.0 : 0.0;
        return;
    }

    // Normalize
    P1 = w1 / sum;
    P2 = w2 / sum;
    P3 = w3 / sum;
}

double RepairMediatedMisrepairModel::cumulativeProbToRate(double p_cum, double T0) const {
    // Map cumulative probability to exponential rate
    // P_cum = 1 - exp(-lambda * T0) => lambda = -log(1 - P_cum) / T0
    p_cum = clamp(p_cum, 0.0, 1.0 - 1e-15);
    T0 = std::max(T0, 1e-9);
    return -std::log(1.0 - p_cum) / T0;
}

double RepairMediatedMisrepairModel::competingRiskDeathProb(double lambda21, double lambda23, double T) const {
    // Competing risks: probability of S2->S3 by time T
    // R23 = (lambda23 / (lambda21 + lambda23)) * (1 - exp(-(lambda21 + lambda23) * T))
    lambda21 = std::max(lambda21, 0.0);
    lambda23 = std::max(lambda23, 0.0);

    const double lambda_sum = lambda21 + lambda23;

    if (lambda_sum <= 1e-15) {
        return 0.0;  // No transitions occur
    }

    const double frac = lambda23 / lambda_sum;
    return frac * (1.0 - std::exp(-lambda_sum * T));
}

double RepairMediatedMisrepairModel::survivalFraction(double dose, const ParamsS2& params) const {
    const double lambda = m_config.kappa * dose;  // Expected DSB count
    const int nmax = poisson_nmax(lambda);

    // Compute kinetic rates from state overlaps (constant for given params)
    // These determine how fast S2 cells transition to S1 (repair) or S3 (death)
    const double p21_cum = overlapScore(params.e2, 0.0);       // S2->S1 overlap
    const double p23_cum = overlapScore(params.e2, params.e3); // baseline S2->S3 overlap (energy-threshold channel)

    // Map cumulative probability to exponential rates over reference windows
    const double lambda21 = cumulativeProbToRate(p21_cum, params.T21);
    const double lambda23_energy = cumulativeProbToRate(p23_cum, params.T23);

    double sf = 0.0;  // Total survival fraction (S1 only)

    for (int n = 0; n <= nmax; ++n) {
        // P(N = n | lambda)
        const double pn = std::exp(poisson_log_pmf(n, lambda));

        // State assignment probabilities given n DSBs
        double P1, P2, P3;
        assignmentProbabilities(n, params, P1, P2, P3);

        // --- Stochastic misrepair channel ---
        // Misrepair hazard: lambda_mis = k_error * n
        // Total S2->S3 hazard = energy-threshold hazard + misrepair hazard
        const double lambda23_tot = lambda23_energy + std::max(0.0, params.k_error) * static_cast<double>(n);

        // Probability that an S2 cell repairs to S1 by assay time window
        const double lambda_sum = std::max(0.0, lambda21) + std::max(0.0, lambda23_tot);
        const double R21_n = (lambda_sum > 1e-15)
            ? (lambda21 / lambda_sum) * (1.0 - std::exp(-lambda_sum * m_config.T_assay_h))
            : 0.0;

        // Survival counts ONLY S1 at assay time
        const double survival_given_n = clamp(P1 + P2 * R21_n, 0.0, 1.0);

        sf += pn * survival_given_n;
    }

    return clamp(sf, m_config.eps_sf, 1.0);
}

std::vector<double> RepairMediatedMisrepairModel::survivalCurve(const std::vector<double>& doses,
                                                        const ParamsS2& params) const {
    std::vector<double> sf_values;
    sf_values.reserve(doses.size());

    for (double dose : doses) {
        sf_values.push_back(survivalFraction(dose, params));
    }

    return sf_values;
}

} // namespace CellStateCalibration
