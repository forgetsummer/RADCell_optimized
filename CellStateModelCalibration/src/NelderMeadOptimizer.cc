/**
 * @file NelderMeadOptimizer.cc
 * @brief Implementation of Nelder-Mead simplex optimizer
 */

#include "NelderMeadOptimizer.hh"
#include "MathUtilities.hh"
#include <algorithm>
#include <cmath>

namespace CellStateCalibration {

NelderMeadOptimizer::NelderMeadOptimizer(const NelderMeadOptions& options)
    : m_options(options)
    , m_lowerBounds(0.0, 0.0)
    , m_upperBounds(200.0, 200.0)
    , m_iterations(0)
{
}

void NelderMeadOptimizer::setOptions(const NelderMeadOptions& options) {
    m_options = options;
}

const NelderMeadOptions& NelderMeadOptimizer::getOptions() const {
    return m_options;
}

void NelderMeadOptimizer::setBounds(const ParamsReduced& lower, const ParamsReduced& upper) {
    m_lowerBounds = lower;
    m_upperBounds = upper;
}

int NelderMeadOptimizer::getIterations() const {
    return m_iterations;
}

ParamsReduced NelderMeadOptimizer::clampToBounds(const ParamsReduced& p) const {
    return ParamsReduced(
        clamp(p.e3, m_lowerBounds.e3, m_upperBounds.e3),
        clamp(p.a, m_lowerBounds.a, m_upperBounds.a)
    );
}

double NelderMeadOptimizer::squaredDistance(const ParamsReduced& a, const ParamsReduced& b) {
    double de3 = a.e3 - b.e3;
    double da = a.a - b.a;
    return de3 * de3 + da * da;
}

ParamsReduced NelderMeadOptimizer::optimize(
    std::function<double(const ParamsReduced&)> objective,
    const ParamsReduced& x0,
    const ParamsReduced& step,
    double& best_value)
{
    m_iterations = 0;
    
    // Initialize simplex: x0, x0 + step_e3, x0 + step_a
    std::vector<ParamsReduced> x(3);
    x[0] = clampToBounds(x0);
    x[1] = clampToBounds(ParamsReduced(x0.e3 + step.e3, x0.a));
    x[2] = clampToBounds(ParamsReduced(x0.e3, x0.a + step.a));
    
    // Evaluate objective at simplex vertices
    std::vector<double> fx(3);
    for (int i = 0; i < 3; ++i) {
        fx[i] = objective(x[i]);
    }
    
    // Lambda to sort simplex by function value
    auto sort_simplex = [&]() {
        std::vector<int> idx = {0, 1, 2};
        std::sort(idx.begin(), idx.end(), [&](int i, int j) {
            return fx[i] < fx[j];
        });
        std::vector<ParamsReduced> x2(3);
        std::vector<double> fx2(3);
        for (int k = 0; k < 3; ++k) {
            x2[k] = x[idx[k]];
            fx2[k] = fx[idx[k]];
        }
        x.swap(x2);
        fx.swap(fx2);
    };
    
    sort_simplex();
    
    // Main optimization loop
    for (int iter = 0; iter < m_options.max_iter; ++iter) {
        m_iterations = iter + 1;
        
        // Check convergence
        double fspan = std::abs(fx[2] - fx[0]);
        double xspan = std::sqrt(std::max({
            squaredDistance(x[2], x[0]),
            squaredDistance(x[1], x[0]),
            squaredDistance(x[2], x[1])
        }));
        
        if (fspan < m_options.tol_f && xspan < m_options.tol_x) {
            break;
        }
        
        // Compute centroid of best two points
        ParamsReduced xc(
            0.5 * (x[0].e3 + x[1].e3),
            0.5 * (x[0].a + x[1].a)
        );
        
        // Reflection
        ParamsReduced xr = clampToBounds(ParamsReduced(
            xc.e3 + m_options.alpha * (xc.e3 - x[2].e3),
            xc.a + m_options.alpha * (xc.a - x[2].a)
        ));
        double fr = objective(xr);
        
        if (fr < fx[0]) {
            // Try expansion
            ParamsReduced xe = clampToBounds(ParamsReduced(
                xc.e3 + m_options.gamma * (xr.e3 - xc.e3),
                xc.a + m_options.gamma * (xr.a - xc.a)
            ));
            double fe = objective(xe);
            
            if (fe < fr) {
                x[2] = xe;
                fx[2] = fe;
            } else {
                x[2] = xr;
                fx[2] = fr;
            }
        } else if (fr < fx[1]) {
            // Accept reflection
            x[2] = xr;
            fx[2] = fr;
        } else {
            // Contraction
            ParamsReduced xcand;
            if (fr < fx[2]) {
                // Outside contraction
                xcand = clampToBounds(ParamsReduced(
                    xc.e3 + m_options.rho * (xr.e3 - xc.e3),
                    xc.a + m_options.rho * (xr.a - xc.a)
                ));
            } else {
                // Inside contraction
                xcand = clampToBounds(ParamsReduced(
                    xc.e3 - m_options.rho * (xc.e3 - x[2].e3),
                    xc.a - m_options.rho * (xc.a - x[2].a)
                ));
            }
            double fc = objective(xcand);
            
            if (fc < fx[2]) {
                x[2] = xcand;
                fx[2] = fc;
            } else {
                // Shrink toward best point
                for (int i = 1; i < 3; ++i) {
                    x[i] = clampToBounds(ParamsReduced(
                        x[0].e3 + m_options.sigma * (x[i].e3 - x[0].e3),
                        x[0].a + m_options.sigma * (x[i].a - x[0].a)
                    ));
                    fx[i] = objective(x[i]);
                }
            }
        }
        
        sort_simplex();
    }
    
    best_value = fx[0];
    return x[0];
}

} // namespace CellStateCalibration
