/**
 * @file RepairMediatedCalibrator.cc
 * @brief Implementation of Repair Mediated model calibrator
 */

#include "RepairMediatedCalibrator.hh"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>

namespace CellStateCalibration {

RepairMediatedCalibrator::RepairMediatedCalibrator(const FitConfigS2& config)
    : m_config(config)
    , m_model(config)
{}

void RepairMediatedCalibrator::setConfig(const FitConfigS2& config) {
    m_config = config;
    m_model.setConfig(config);
}

void RepairMediatedCalibrator::setData(const std::vector<DataPoint>& data) {
    m_data = data;
}

void RepairMediatedCalibrator::addDataPoint(const DataPoint& point) {
    m_data.push_back(point);
}

void RepairMediatedCalibrator::clearData() {
    m_data.clear();
}

void RepairMediatedCalibrator::setDefaultUncertainties(double rel_err) {
    for (auto& dp : m_data) {
        dp.si = rel_err;
    }
}

double RepairMediatedCalibrator::negativeLogLikelihood(const ParamsS2& params) const {
    double nll = 0.0;
    
    for (const auto& dp : m_data) {
        // Skip data points beyond fit range
        if (dp.D > m_config.dose_max_fit) continue;
        
        // Model prediction
        double sf_model = m_model.survivalFraction(dp.D, params);
        
        // Log-space likelihood (normal on ln SF)
        double y = std::log(clamp(dp.SF_obs, m_config.eps_sf, 1.0));
        double mu = std::log(sf_model);
        double s = std::max(dp.si, 1e-6);
        
        double residual = (y - mu) / s;
        nll += 0.5 * residual * residual + std::log(s);
    }
    
    return nll;
}

// Helper: squared distance between vectors
static double vec_dist2(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return s;
}

std::vector<double> RepairMediatedCalibrator::nelderMeadND(
    std::function<double(const std::vector<double>&)> f,
    const std::vector<double>& x0,
    const std::vector<double>& step,
    const std::vector<double>& lo,
    const std::vector<double>& hi,
    double& best_f_out,
    int& iterations_out,
    bool& converged_out) const
{
    // Nelder-Mead options
    const int max_iter = 8000;
    const double tol_f = 1e-10;
    const double tol_x = 1e-8;
    const double alpha = 1.0;   // Reflection
    const double gamma = 2.0;   // Expansion
    const double rho = 0.5;     // Contraction
    const double sigma = 0.5;   // Shrink
    
    const int n = static_cast<int>(x0.size());
    const int m = n + 1;  // Simplex size
    
    // Initialize simplex
    std::vector<std::vector<double>> x(m, x0);
    for (int i = 1; i < m; ++i) {
        x[i][i-1] += step[i-1];
        x[i][i-1] = clamp(x[i][i-1], lo[i-1], hi[i-1]);
    }
    
    // Evaluate function at all vertices
    std::vector<double> fx(m);
    for (int i = 0; i < m; ++i) {
        fx[i] = f(x[i]);
    }
    
    // Sort simplex by function value
    auto sort_simplex = [&]() {
        std::vector<int> idx(m);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return fx[a] < fx[b]; });
        
        std::vector<std::vector<double>> x2(m);
        std::vector<double> fx2(m);
        for (int k = 0; k < m; ++k) {
            x2[k] = x[idx[k]];
            fx2[k] = fx[idx[k]];
        }
        x.swap(x2);
        fx.swap(fx2);
    };
    
    sort_simplex();
    
    converged_out = false;
    iterations_out = 0;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        iterations_out = iter + 1;
        
        // Check convergence
        double fspan = std::abs(fx[m-1] - fx[0]);
        double xspan = 0.0;
        for (int i = 1; i < m; ++i) {
            xspan = std::max(xspan, std::sqrt(vec_dist2(x[i], x[0])));
        }
        
        if (fspan < tol_f && xspan < tol_x) {
            converged_out = true;
            break;
        }
        
        // Compute centroid of best n points (exclude worst)
        std::vector<double> xc(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int d = 0; d < n; ++d) {
                xc[d] += x[i][d];
            }
        }
        for (int d = 0; d < n; ++d) {
            xc[d] /= static_cast<double>(n);
        }
        
        // Bound helper
        auto bounded = [&](std::vector<double> v) {
            for (int d = 0; d < n; ++d) {
                v[d] = clamp(v[d], lo[d], hi[d]);
            }
            return v;
        };
        
        // Reflection
        std::vector<double> xr(n);
        for (int d = 0; d < n; ++d) {
            xr[d] = xc[d] + alpha * (xc[d] - x[m-1][d]);
        }
        xr = bounded(xr);
        double fr = f(xr);
        
        if (fr < fx[0]) {
            // Expansion
            std::vector<double> xe(n);
            for (int d = 0; d < n; ++d) {
                xe[d] = xc[d] + gamma * (xr[d] - xc[d]);
            }
            xe = bounded(xe);
            double fe = f(xe);
            
            if (fe < fr) {
                x[m-1] = xe;
                fx[m-1] = fe;
            } else {
                x[m-1] = xr;
                fx[m-1] = fr;
            }
        } else if (fr < fx[n-1]) {
            // Accept reflection
            x[m-1] = xr;
            fx[m-1] = fr;
        } else {
            // Contraction
            std::vector<double> xk(n);
            if (fr < fx[m-1]) {
                // Outside contraction
                for (int d = 0; d < n; ++d) {
                    xk[d] = xc[d] + rho * (xr[d] - xc[d]);
                }
            } else {
                // Inside contraction
                for (int d = 0; d < n; ++d) {
                    xk[d] = xc[d] - rho * (xc[d] - x[m-1][d]);
                }
            }
            xk = bounded(xk);
            double fk = f(xk);
            
            if (fk < fx[m-1]) {
                x[m-1] = xk;
                fx[m-1] = fk;
            } else {
                // Shrink
                for (int i = 1; i < m; ++i) {
                    for (int d = 0; d < n; ++d) {
                        x[i][d] = x[0][d] + sigma * (x[i][d] - x[0][d]);
                    }
                    x[i] = bounded(x[i]);
                    fx[i] = f(x[i]);
                }
            }
        }
        
        sort_simplex();
    }
    
    best_f_out = fx[0];
    return x[0];
}

RepairMediatedFitResult RepairMediatedCalibrator::fit(
    const ParamsS2& initial_guess,
    const ParamsS2& step)
{
    RepairMediatedFitResult result;
    
    // Pack parameters into vector - ALL 5 parameters (e2, e3, a, T21, T23)
    std::vector<double> x0 = {initial_guess.e2, initial_guess.e3, initial_guess.a,
                              initial_guess.T21, initial_guess.T23};
    // Only optimize e2, e3, a (3 parameters)
    // T21 and T23 are FIXED from config
    std::vector<double> x0_3 = {initial_guess.e2, initial_guess.e3, initial_guess.a};
    std::vector<double> step_vec = {step.e2, step.e3, step.a};
    
    // Parameter bounds for e2, e3, a
    std::vector<double> lo = {0.1, 0.2, 0.001};
    std::vector<double> hi = {80.0, 200.0, 50.0};
    
    // Get fixed T21 and T23 from config
    const double T21_fixed = m_config.T21_fixed;
    const double T23_fixed = m_config.T23_fixed;
    
    // Objective function - only 3 parameters optimized, T21/T23 are fixed
    auto objective = [this, T21_fixed, T23_fixed](const std::vector<double>& v) -> double {
        ParamsS2 p;
        p.e2 = v[0];
        p.e3 = v[1];
        p.a = v[2];
        p.T21 = T21_fixed;  // FIXED from config
        p.T23 = T23_fixed;  // FIXED from config
        
        // Hard constraints
        if (!p.isValid()) {
            return 1e30;
        }
        
        return negativeLogLikelihood(p);
    };
    
    // Run optimization (3 parameters)
    double best_nll = 0.0;
    int iterations = 0;
    bool converged = false;
    
    std::vector<double> xhat = nelderMeadND(objective, x0_3, step_vec, lo, hi,
                                             best_nll, iterations, converged);
    
    // Unpack result - e2, e3, a optimized; T21, T23 fixed
    result.params.e2 = xhat[0];
    result.params.e3 = xhat[1];
    result.params.a = xhat[2];
    result.params.T21 = T21_fixed;  // From config
    result.params.T23 = T23_fixed;  // From config
    result.negLogLikelihood = best_nll;
    result.iterations = iterations;
    result.converged = converged;
    
    return result;
}

std::vector<double> RepairMediatedCalibrator::predictSurvival(
    const ParamsS2& params,
    const std::vector<double>& doses) const
{
    return m_model.survivalCurve(doses, params);
}

const std::vector<DataPoint>& RepairMediatedCalibrator::getData() const {
    return m_data;
}

const RepairMediatedModel& RepairMediatedCalibrator::getModel() const {
    return m_model;
}

} // namespace CellStateCalibration
