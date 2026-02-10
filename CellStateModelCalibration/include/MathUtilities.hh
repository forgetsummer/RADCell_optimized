#ifndef MATH_UTILITIES_HH
#define MATH_UTILITIES_HH

/**
 * @file MathUtilities.hh
 * @brief Mathematical utility functions for cell state model calibration
 * 
 * Provides standard mathematical functions including:
 * - Normal CDF using error function
 * - Clamping utilities
 * - Poisson distribution functions (log-space for numerical stability)
 */

#include <cmath>
#include <algorithm>
#include <limits>

namespace CellStateCalibration {

/**
 * @brief Standard normal CDF Φ(x) using complementary error function
 * @param x Input value
 * @return P(Z <= x) where Z ~ N(0,1)
 */
inline double normal_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

/**
 * @brief Clamp value to specified range
 * @param x Value to clamp
 * @param lo Lower bound
 * @param hi Upper bound
 * @return Clamped value in [lo, hi]
 */
inline double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

/**
 * @brief Stable log factorial using lgamma
 * @param n Non-negative integer
 * @return log(n!)
 */
inline double log_factorial(int n) {
    return std::lgamma(n + 1.0);
}

/**
 * @brief Poisson PMF in log form for numerical stability
 * 
 * Computes log P(N=n) = -λ + n*log(λ) - log(n!)
 * 
 * @param n Number of events
 * @param lambda Poisson rate parameter
 * @return log P(N=n | λ)
 */
inline double poisson_log_pmf(int n, double lambda) {
    if (lambda <= 0.0) {
        return (n == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
    }
    return -lambda + n * std::log(lambda) - log_factorial(n);
}

/**
 * @brief Determine truncation point for Poisson sum
 * 
 * Returns n_max such that P(N > n_max) is negligible.
 * Uses rule: n_max ~ λ + 6*sqrt(λ)
 * 
 * @param lambda Poisson rate parameter
 * @return Maximum n to consider in summations
 */
inline int poisson_nmax(double lambda) {
    if (lambda < 1e-8) return 5;
    int nmax = static_cast<int>(std::ceil(lambda + 6.0 * std::sqrt(lambda)));
    return std::max(10, nmax);
}

} // namespace CellStateCalibration

#endif // MATH_UTILITIES_HH
