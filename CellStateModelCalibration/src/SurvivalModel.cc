/**
 * @file SurvivalModel.cc
 * @brief Implementation of cell survival model
 */

#include "SurvivalModel.hh"
#include <cmath>

namespace CellStateCalibration {

SurvivalModel::SurvivalModel(const FitConfig& config)
    : m_config(config)
{
}

void SurvivalModel::setConfig(const FitConfig& config) {
    m_config = config;
}

const FitConfig& SurvivalModel::getConfig() const {
    return m_config;
}

double SurvivalModel::transitionProbability(int n, const ParamsReduced& params) const {
    // p13(n) = 2 * Φ(-|a*n - e3| / 2)
    // This is the probability of S1->S3 transition given n DSBs
    double z = -std::abs(params.a * static_cast<double>(n) - params.e3) / 2.0;
    double val = 2.0 * normal_cdf(z);
    return clamp(val, 0.0, 1.0);
}

double SurvivalModel::survivalFraction(double dose, const ParamsReduced& params) const {
    // SF(D) = 1 - Σ_n P(N=n | κD) * p13(n)
    // where N ~ Poisson(κD) is the DSB count
    
    const double lambda = m_config.kappa * dose;
    const int nmax = poisson_nmax(lambda);
    
    double death_prob = 0.0;
    
    // Sum over Poisson distribution
    for (int n = 0; n <= nmax; ++n) {
        double log_pn = poisson_log_pmf(n, lambda);
        double pn = std::exp(log_pn);
        death_prob += pn * transitionProbability(n, params);
    }
    
    double sf = 1.0 - death_prob;
    return clamp(sf, m_config.eps_sf, 1.0);
}

std::vector<double> SurvivalModel::survivalCurve(const std::vector<double>& doses,
                                                  const ParamsReduced& params) const {
    std::vector<double> sf_values;
    sf_values.reserve(doses.size());
    
    for (double dose : doses) {
        sf_values.push_back(survivalFraction(dose, params));
    }
    
    return sf_values;
}

} // namespace CellStateCalibration
