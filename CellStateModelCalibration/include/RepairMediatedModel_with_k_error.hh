#ifndef REPAIR_MEDIATED_MODEL_WITH_K_ERROR_HH
#define REPAIR_MEDIATED_MODEL_WITH_K_ERROR_HH

/**
 * @file RepairMediatedModel_with_k_error.hh
 * @brief Repair Mediated survival model with stochastic misrepair (S1->S2->S3 transitions)
 *
 * This model includes an intermediate repair state S2 between the healthy
 * state S1 and the dead state S3. After radiation damage:
 *
 * 1. Cells can be assigned to S1, S2, or S3 based on damage level
 * 2. From S2, cells compete between:
 *    - Repair (S2->S1) with timescale T21
 *    - Death (S2->S3) with timescale T23 plus stochastic misrepair channel
 *
 * This captures the biological reality that moderately damaged cells
 * may either repair and survive, or fail to repair and die.
 *
 * The stochastic misrepair channel adds:
 *   lambda_mis = k_error * N_dsb
 * This represents DSB-dependent misrepair events leading to cell death.
 *
 * The model uses:
 * - Poisson distribution for DSB counts
 * - Gaussian overlap for state assignment probabilities
 * - Competing exponential kinetics for S2 fate (including stochastic misrepair term)
 */

#include "DataTypes.hh"
#include "MathUtilities.hh"

namespace CellStateCalibration {

/**
 * @brief Repair Mediated survival model with S2 kinetics and stochastic misrepair
 *
 * Computes survival fraction considering:
 * - Immediate assignment to S1, S2, or S3 based on damage
 * - Time-dependent fate of S2 cells (repair vs death with misrepair channel)
 */
class RepairMediatedMisrepairModel {
public:
    /**
     * @brief Constructor with optional configuration
     */
    explicit RepairMediatedMisrepairModel(const FitConfigS2& config = FitConfigS2());

    /**
     * @brief Set configuration parameters
     */
    void setConfig(const FitConfigS2& config);

    /**
     * @brief Get current configuration
     */
    const FitConfigS2& getConfig() const;

    /**
     * @brief Compute Gaussian overlap score between two states
     *
     * w(e, ej) = 2 * Phi(-|e - ej| / 2)
     *
     * This measures how much the damage energy e overlaps with state j
     * centered at ej. Higher overlap means higher probability of assignment.
     *
     * @param e Reduced damage energy (alpha * n / sigma)
     * @param ej Reduced state energy (Ej / sigma)
     * @return Overlap score in [0, 1]
     */
    double overlapScore(double e, double ej) const;

    /**
     * @brief Compute normalized state assignment probabilities
     *
     * Given n DSBs, computes P(S1), P(S2), P(S3) based on Gaussian overlaps.
     *
     * @param n Number of DSBs
     * @param params Model parameters
     * @param P1 Output: Probability of S1 (healthy)
     * @param P2 Output: Probability of S2 (repair)
     * @param P3 Output: Probability of S3 (dead)
     */
    void assignmentProbabilities(int n, const ParamsS2& params,
                                  double& P1, double& P2, double& P3) const;

    /**
     * @brief Convert cumulative probability to exponential rate
     *
     * Maps a cumulative transition probability over characteristic time T0
     * to an equivalent exponential rate.
     *
     * @param p_cum Cumulative probability (0 to 1)
     * @param T0 Characteristic time
     * @return Exponential rate lambda
     */
    double cumulativeProbToRate(double p_cum, double T0) const;

    /**
     * @brief Compute probability of S2->S3 transition under competing risks
     *
     * Starting in S2, computes probability of reaching S3 (death) by time T,
     * given competing repair (S2->S1) and death (S2->S3) processes.
     *
     * R23(T) = (lambda23 / (lambda21 + lambda23)) * (1 - exp(-(lambda21 + lambda23) * T))
     *
     * @param lambda21 Repair rate
     * @param lambda23 Death rate (including misrepair contribution)
     * @param T Time window
     * @return Probability of death from S2
     */
    double competingRiskDeathProb(double lambda21, double lambda23, double T) const;

    /**
     * @brief Compute survival fraction at given dose
     *
     * SF(D) = Sum_n P(N=n|lambda=kappa*D) * [P1(n) + P2(n) * R21]
     *
     * The S2 death rate includes both energy-threshold and misrepair channels:
     *   lambda23_tot = lambda23_energy + k_error * n
     *
     * R21 is the competing-risk probability that an S2 cell repairs to S1 by assay time.
     *
     * @param dose Dose in Gy
     * @param params Model parameters (including k_error)
     * @return Survival fraction
     */
    double survivalFraction(double dose, const ParamsS2& params) const;

    /**
     * @brief Compute survival curve for multiple doses
     *
     * @param doses Vector of dose values
     * @param params Model parameters
     * @return Vector of survival fractions
     */
    std::vector<double> survivalCurve(const std::vector<double>& doses,
                                       const ParamsS2& params) const;

private:
    FitConfigS2 m_config;
};

} // namespace CellStateCalibration

#endif // REPAIR_MEDIATED_MODEL_WITH_K_ERROR_HH
