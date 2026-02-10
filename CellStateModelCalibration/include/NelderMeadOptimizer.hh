#ifndef NELDER_MEAD_OPTIMIZER_HH
#define NELDER_MEAD_OPTIMIZER_HH

/**
 * @file NelderMeadOptimizer.hh
 * @brief Nelder-Mead simplex optimizer for 2D parameter space
 * 
 * Implements the Nelder-Mead algorithm for minimizing functions
 * of two variables. This is a derivative-free method that is
 * robust and easy to implement.
 */

#include "DataTypes.hh"
#include <functional>

namespace CellStateCalibration {

/**
 * @class NelderMeadOptimizer
 * @brief 2D Nelder-Mead simplex optimizer
 * 
 * The algorithm maintains a simplex of 3 points in 2D and
 * iteratively improves by reflection, expansion, contraction,
 * and shrinkage operations.
 */
class NelderMeadOptimizer {
public:
    /**
     * @brief Constructor with options
     * @param options Optimizer configuration
     */
    explicit NelderMeadOptimizer(const NelderMeadOptions& options = NelderMeadOptions());
    
    /**
     * @brief Set optimizer options
     * @param options New options
     */
    void setOptions(const NelderMeadOptions& options);
    
    /**
     * @brief Get current options
     * @return Current options
     */
    const NelderMeadOptions& getOptions() const;
    
    /**
     * @brief Set parameter bounds
     * @param lower Lower bounds for (e3, a)
     * @param upper Upper bounds for (e3, a)
     */
    void setBounds(const ParamsReduced& lower, const ParamsReduced& upper);
    
    /**
     * @brief Optimize the objective function
     * 
     * @param objective Function to minimize: f(ParamsReduced) -> double
     * @param x0 Initial guess
     * @param step Initial step sizes for simplex construction
     * @param[out] best_value Optimal function value found
     * @return Optimal parameters
     */
    ParamsReduced optimize(
        std::function<double(const ParamsReduced&)> objective,
        const ParamsReduced& x0,
        const ParamsReduced& step,
        double& best_value);
    
    /**
     * @brief Get number of iterations used in last optimization
     * @return Iteration count
     */
    int getIterations() const;

private:
    /**
     * @brief Clamp parameters to bounds
     */
    ParamsReduced clampToBounds(const ParamsReduced& p) const;
    
    /**
     * @brief Compute squared distance between two points
     */
    static double squaredDistance(const ParamsReduced& a, const ParamsReduced& b);
    
    NelderMeadOptions m_options;
    ParamsReduced m_lowerBounds;
    ParamsReduced m_upperBounds;
    int m_iterations;
};

} // namespace CellStateCalibration

#endif // NELDER_MEAD_OPTIMIZER_HH
