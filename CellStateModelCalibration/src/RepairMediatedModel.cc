/**
 * @file RepairMediatedModel.cc
 * @brief Implementation of Repair Mediated survival model
 */

#include "RepairMediatedModel.hh"
#include <cmath>
#include <algorithm>

namespace CellStateCalibration {

RepairMediatedModel::RepairMediatedModel(const FitConfigS2& config)
    : m_config(config)
{}

void RepairMediatedModel::setConfig(const FitConfigS2& config) {
    m_config = config;
}

const FitConfigS2& RepairMediatedModel::getConfig() const {
    return m_config;
}

double RepairMediatedModel::overlapScore(double e, double ej) const {
    // w(e, ej) = 2 * Φ(-|e - ej| / 2)
    double z = -std::abs(e - ej) / 2.0;
    double val = 2.0 * normal_cdf(z);
    return clamp(val, 0.0, 1.0);
}

void RepairMediatedModel::assignmentProbabilities(int n, const ParamsS2& params,
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

double RepairMediatedModel::cumulativeProbToRate(double p_cum, double T0) const {
    // Map cumulative probability to exponential rate
    // P_cum = 1 - exp(-λ * T0) => λ = -log(1 - P_cum) / T0
    p_cum = clamp(p_cum, 0.0, 1.0 - 1e-15);
    T0 = std::max(T0, 1e-9);
    return -std::log(1.0 - p_cum) / T0;
}

double RepairMediatedModel::competingRiskDeathProb(double lambda21, double lambda23, double T) const {
    // Competing risks: probability of S2->S3 by time T
    // R23 = (λ23 / (λ21 + λ23)) * (1 - exp(-(λ21 + λ23) * T))
    lambda21 = std::max(lambda21, 0.0);
    lambda23 = std::max(lambda23, 0.0);
    
    const double lambda_sum = lambda21 + lambda23;
    
    if (lambda_sum <= 1e-15) {
        return 0.0;  // No transitions occur
    }
    
    const double frac = lambda23 / lambda_sum;
    return frac * (1.0 - std::exp(-lambda_sum * T));
}

double RepairMediatedModel::survivalFraction(double dose, const ParamsS2& params) const {
    const double lambda = m_config.kappa * dose;  // Expected DSB count
    const int nmax = poisson_nmax(lambda);
    
    // Compute kinetic rates from state overlaps (constant for given params)
    // These determine how fast S2 cells transition to S1 (repair) or S3 (death)
    const double p21_cum = overlapScore(params.e2, 0.0);       // S2->S1 overlap
    const double p23_cum = overlapScore(params.e2, params.e3); // S2->S3 overlap
    
    // Use T21 and T23 from params (OPTIMIZED values)
    const double lambda21 = cumulativeProbToRate(p21_cum, params.T21);
    const double lambda23 = cumulativeProbToRate(p23_cum, params.T23);
    
    // Probability of S2->S3 over assay time window
    const double R23 = competingRiskDeathProb(lambda21, lambda23, m_config.T_assay_h);
    
    double death_prob = 0.0;  // Total death probability
    
    for (int n = 0; n <= nmax; ++n) {
        // P(N = n | λ)
        const double pn = std::exp(poisson_log_pmf(n, lambda));
        
        // State assignment probabilities given n DSBs
        double P1, P2, P3;
        assignmentProbabilities(n, params, P1, P2, P3);
        
        // Total death probability given n:
        // - Immediate death (P3)
        // - Death from repair state (P2 * R23)
        const double death_given_n = clamp(P3 + P2 * R23, 0.0, 1.0);
        
        death_prob += pn * death_given_n;
    }
    
    double sf = 1.0 - death_prob;
    return clamp(sf, m_config.eps_sf, 1.0);
}

std::vector<double> RepairMediatedModel::survivalCurve(const std::vector<double>& doses,
                                                        const ParamsS2& params) const {
    std::vector<double> sf_values;
    sf_values.reserve(doses.size());
    
    for (double dose : doses) {
        sf_values.push_back(survivalFraction(dose, params));
    }
    
    return sf_values;
}

} // namespace CellStateCalibration
