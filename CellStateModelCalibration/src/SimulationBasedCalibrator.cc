/**
 * @file SimulationBasedCalibrator.cc
 * @brief Implementation of simulation-based calibrator
 */

#include "SimulationBasedCalibrator.hh"
#include "MathUtilities.hh"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iostream>

namespace CellStateCalibration {

SimulationBasedCalibrator::SimulationBasedCalibrator(const FitConfigS2& config)
    : m_config(config)
{}

void SimulationBasedCalibrator::setConfig(const FitConfigS2& config) {
    m_config = config;
}

void SimulationBasedCalibrator::setSurvivalFunction(SurvivalFunction fn) {
    m_survivalFn = std::move(fn);
}

void SimulationBasedCalibrator::setBatchSurvivalFunction(BatchSurvivalFunction fn) {
    m_batchSurvivalFn = std::move(fn);
}

void SimulationBasedCalibrator::setData(const std::vector<DataPoint>& data) {
    m_data = data;
}

void SimulationBasedCalibrator::addDataPoint(const DataPoint& point) {
    m_data.push_back(point);
}

void SimulationBasedCalibrator::clearData() {
    m_data.clear();
}

void SimulationBasedCalibrator::setDefaultUncertainties(double rel_err) {
    for (auto& dp : m_data) {
        dp.si = rel_err;
    }
}

double SimulationBasedCalibrator::negativeLogLikelihood(const ParamsS2& params) const {
    if (!m_survivalFn && !m_batchSurvivalFn) return 1e30;

    // Collect doses that are within fit range
    std::vector<size_t> indices;
    std::vector<double> doses;
    for (size_t i = 0; i < m_data.size(); i++) {
        if (m_data[i].D <= m_config.dose_max_fit) {
            indices.push_back(i);
            doses.push_back(m_data[i].D);
        }
    }
    if (doses.empty()) return 0.0;

    // Get SF predictions (batch if available for parallelism)
    std::vector<double> sf_models;
    if (m_batchSurvivalFn) {
        sf_models = m_batchSurvivalFn(doses, params);
    } else {
        sf_models.resize(doses.size());
        for (size_t i = 0; i < doses.size(); i++)
            sf_models[i] = m_survivalFn(doses[i], params);
    }

    double nll = 0.0;
    for (size_t k = 0; k < indices.size(); k++) {
        const auto& dp = m_data[indices[k]];
        double sf_model = clamp(sf_models[k], m_config.eps_sf, 1.0);
        double y = std::log(clamp(dp.SF_obs, m_config.eps_sf, 1.0));
        double mu = std::log(sf_model);
        double s = std::max(dp.si, 1e-6);
        double residual = (y - mu) / s;
        nll += 0.5 * residual * residual + std::log(s);
    }
    return nll;
}

static double vec_dist2_simcalib(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return s;
}

std::vector<double> SimulationBasedCalibrator::nelderMeadND(
    std::function<double(const std::vector<double>&)> f,
    const std::vector<double>& x0,
    const std::vector<double>& step,
    const std::vector<double>& lo,
    const std::vector<double>& hi,
    double& best_f_out,
    int& iterations_out,
    bool& converged_out,
    int max_iter) const
{
    const double tol_f = 1e-6;
    const double tol_x = 1e-5;
    const double alpha_nm = 1.0;
    const double gamma = 2.0;
    const double rho = 0.5;
    const double sigma = 0.5;

    const int n = static_cast<int>(x0.size());
    const int m = n + 1;

    std::vector<std::vector<double>> x(m, x0);
    for (int i = 1; i < m; ++i) {
        x[i][i-1] += step[i-1];
        x[i][i-1] = clamp(x[i][i-1], lo[i-1], hi[i-1]);
    }

    std::vector<double> fx(m);
    for (int i = 0; i < m; ++i) {
        fx[i] = f(x[i]);
    }

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

        double fspan = std::abs(fx[m-1] - fx[0]);
        double xspan = 0.0;
        for (int i = 1; i < m; ++i) {
            xspan = std::max(xspan, std::sqrt(vec_dist2_simcalib(x[i], x[0])));
        }

        if (fspan < tol_f && xspan < tol_x) {
            converged_out = true;
            break;
        }

        // Progress logging every 50 iterations
        if (iter % 50 == 0) {
            std::cout << "    NM iter " << iter << ": best NLL = " << fx[0]
                      << ", fspan = " << fspan << std::endl;
        }

        std::vector<double> xc(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int d = 0; d < n; ++d) {
                xc[d] += x[i][d];
            }
        }
        for (int d = 0; d < n; ++d) {
            xc[d] /= static_cast<double>(n);
        }

        auto bounded = [&](std::vector<double> v) {
            for (int d = 0; d < n; ++d) {
                v[d] = clamp(v[d], lo[d], hi[d]);
            }
            return v;
        };

        std::vector<double> xr(n);
        for (int d = 0; d < n; ++d) {
            xr[d] = xc[d] + alpha_nm * (xc[d] - x[m-1][d]);
        }
        xr = bounded(xr);
        double fr = f(xr);

        if (fr < fx[0]) {
            std::vector<double> xe(n);
            for (int d = 0; d < n; ++d) {
                xe[d] = xc[d] + gamma * (xr[d] - xc[d]);
            }
            xe = bounded(xe);
            double fe = f(xe);
            if (fe < fr) { x[m-1] = xe; fx[m-1] = fe; }
            else { x[m-1] = xr; fx[m-1] = fr; }
        } else if (fr < fx[n-1]) {
            x[m-1] = xr; fx[m-1] = fr;
        } else {
            std::vector<double> xk(n);
            if (fr < fx[m-1]) {
                for (int d = 0; d < n; ++d) xk[d] = xc[d] + rho * (xr[d] - xc[d]);
            } else {
                for (int d = 0; d < n; ++d) xk[d] = xc[d] - rho * (xc[d] - x[m-1][d]);
            }
            xk = bounded(xk);
            double fk = f(xk);
            if (fk < fx[m-1]) {
                x[m-1] = xk; fx[m-1] = fk;
            } else {
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

SimulationBasedFitResult SimulationBasedCalibrator::fit(
    const ParamsS2& initial_guess,
    const ParamsS2& step,
    int max_iter)
{
    SimulationBasedFitResult result;

    std::vector<double> x0 = {initial_guess.e2, initial_guess.e3, initial_guess.a,
                              initial_guess.T21, initial_guess.T23, initial_guess.k_error};
    std::vector<double> s  = {step.e2, step.e3, step.a,
                              step.T21, step.T23, step.k_error};

    std::vector<double> lo = {0.1,  0.2,   0.001,  1.0,   1.0,   1e-10};
    std::vector<double> hi = {80.0, 200.0, 50.0,   100.0, 100.0, 1e-2};

    auto objective = [this](const std::vector<double>& v) -> double {
        ParamsS2 p;
        p.e2      = v[0];
        p.e3      = v[1];
        p.a       = v[2];
        p.T21     = v[3];
        p.T23     = v[4];
        p.k_error = v[5];
        if (!p.isValid()) return 1e30;
        return negativeLogLikelihood(p);
    };

    double best_nll = 0.0;
    int iterations = 0;
    bool converged = false;

    std::vector<double> xhat = nelderMeadND(objective, x0, s, lo, hi,
                                             best_nll, iterations, converged, max_iter);

    result.params.e2      = xhat[0];
    result.params.e3      = xhat[1];
    result.params.a       = xhat[2];
    result.params.T21     = xhat[3];
    result.params.T23     = xhat[4];
    result.params.k_error = xhat[5];
    result.negLogLikelihood = best_nll;
    result.iterations = iterations;
    result.converged = converged;

    return result;
}

const std::vector<DataPoint>& SimulationBasedCalibrator::getData() const {
    return m_data;
}

} // namespace CellStateCalibration
