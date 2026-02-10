/**
 * @file DirectLethalModel.cc
 * @brief Implementation of Direct Lethal Model (S1->S3 only)
 */

#include "DirectLethalModel.hh"
#include <cmath>

namespace CellStateCalibration {

DirectLethalModel::DirectLethalModel(const FitConfig& config)
    : m_config(config)
{
}

void DirectLethalModel::setConfig(const FitConfig& config) {
    m_config = config;
}

const FitConfig& DirectLethalModel::getConfig() const {
    return m_config;
}

double DirectLethalModel::transitionProbability_S1_S3(int n, const ParamsReduced& params) const {
    // p13(n) = 2 * Φ(-|a*n - e3| / 2)
    // This is the probability of S1->S3 transition given n DSBs
    double z = -std::abs(params.a * static_cast<double>(n) - params.e3) / 2.0;
    double val = 2.0 * normal_cdf(z);
    return clamp(val, 0.0, 1.0);
}

double DirectLethalModel::survivalFraction(double dose, const ParamsReduced& params) const {
    // SF(D) = 1 - Σ_n P(N=n | κD) * p13(n)
    // where N ~ Poisson(κD) is the DSB count
    
    const double lambda = m_config.kappa * dose;
    const int nmax = poisson_nmax(lambda);
    
    double death_prob = 0.0;
    
    // Sum over Poisson distribution
    for (int n = 0; n <= nmax; ++n) {
        double log_pn = poisson_log_pmf(n, lambda);
        double pn = std::exp(log_pn);
        death_prob += pn * transitionProbability_S1_S3(n, params);
    }
    
    double sf = 1.0 - death_prob;
    return clamp(sf, m_config.eps_sf, 1.0);
}

std::vector<double> DirectLethalModel::survivalCurve(const std::vector<double>& doses,
                                                      const ParamsReduced& params) const {
    std::vector<double> sf_values;
    sf_values.reserve(doses.size());
    
    for (double dose : doses) {
        sf_values.push_back(survivalFraction(dose, params));
    }
    
    return sf_values;
}

} // namespace CellStateCalibration
