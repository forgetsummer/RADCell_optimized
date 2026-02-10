#ifndef DIRECT_LETHAL_MODEL_HH
#define DIRECT_LETHAL_MODEL_HH

/**
 * @file DirectLethalModel.hh
 * @brief Direct Lethal Model: S1 -> S3 transition only
 * 
 * This model considers only the direct path from healthy (S1) to dead (S3):
 * - S1 (healthy) -> S3 (dead) based on DSB damage exceeding threshold
 * - No intermediate repair state (S2) is considered
 * 
 * Use this for:
 * - Simple survival curve fitting
 * - When repair kinetics are not important
 * - Initial parameter estimation
 * 
 * For models including repair (S1->S2->S3), see RepairMediatedModel.
 */

#include "DataTypes.hh"
#include "MathUtilities.hh"

namespace CellStateCalibration {

/**
 * @class DirectLethalModel
 * @brief Forward model for direct lethal (S1->S3) cell death
 * 
 * The model assumes:
 * 1. DSB damage follows Poisson distribution: N ~ Poisson(κD)
 * 2. Damage state X = α*N (linear in DSB count)
 * 3. Cell dies (S1->S3) if X exceeds threshold E3
 * 4. Transition probability uses Gaussian overlap formula:
 *    p13(n) = 2 * Φ(-|a*n - e3| / 2)
 * 
 * Parameters (reduced/identifiable):
 * - e3 = E3/σ : dimensionless death threshold
 * - a = α/σ   : dimensionless damage sensitivity
 */
class DirectLethalModel {
public:
    /**
     * @brief Constructor with configuration
     * @param config Fitting configuration parameters
     */
    explicit DirectLethalModel(const FitConfig& config = FitConfig());
    
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
     * @param params Model parameters (e3, a)
     * @return Transition probability in [0, 1]
     */
    double transitionProbability_S1_S3(int n, const ParamsReduced& params) const;
    
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

#endif // DIRECT_LETHAL_MODEL_HH
