#ifndef SURVIVAL_MODEL_HH
#define SURVIVAL_MODEL_HH

/**
 * @file SurvivalModel.hh
 * @brief Cell survival model based on state transitions
 * 
 * Implements the forward model for cell survival:
 * - S1 (healthy) -> S3 (dead) transition based on DSB damage
 * - Survival fraction as function of dose
 * - Transition probability based on Gaussian overlap
 */

#include "DataTypes.hh"
#include "MathUtilities.hh"

namespace CellStateCalibration {

/**
 * @class SurvivalModel
 * @brief Forward model for cell survival prediction
 * 
 * The model assumes:
 * 1. DSB damage follows Poisson distribution: N ~ Poisson(κD)
 * 2. Damage state X = α*N (linear in DSB count)
 * 3. Cell dies (S1->S3) if X exceeds threshold E3
 * 4. Transition probability uses Gaussian overlap formula
 */
class SurvivalModel {
public:
    /**
     * @brief Constructor with configuration
     * @param config Fitting configuration parameters
     */
    explicit SurvivalModel(const FitConfig& config = FitConfig());
    
    /**
     * @brief Set configuration
     * @param config New configuration
     */
    void setConfig(const FitConfig& config);
    
    /**
     * @brief Get current configuration
     * @return Current configuration
     */
    const FitConfig& getConfig() const;
    
    /**
     * @brief Calculate S1->S3 transition probability for given DSB count
     * 
     * Uses the Gaussian overlap formula:
     * p13(n) = 2 * Φ(-|a*n - e3| / 2)
     * 
     * @param n Number of DSBs
     * @param params Model parameters
     * @return Transition probability in [0, 1]
     */
    double transitionProbability(int n, const ParamsReduced& params) const;
    
    /**
     * @brief Calculate survival fraction at given dose
     * 
     * SF(D) = 1 - Σ_n P(N=n | κD) * p13(n)
     * 
     * @param dose Dose in Gy
     * @param params Model parameters
     * @return Survival fraction in [eps_sf, 1]
     */
    double survivalFraction(double dose, const ParamsReduced& params) const;
    
    /**
     * @brief Calculate survival curve for multiple doses
     * @param doses Vector of doses
     * @param params Model parameters
     * @return Vector of survival fractions
     */
    std::vector<double> survivalCurve(const std::vector<double>& doses, 
                                       const ParamsReduced& params) const;

private:
    FitConfig m_config;
};

} // namespace CellStateCalibration

#endif // SURVIVAL_MODEL_HH
