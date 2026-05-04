/**
 * @file RepairMediatedCalibrator_multicomponent.cc
 * @brief Implementation of the multi-component energy calibrator
 */

#include "RepairMediatedCalibrator_multicomponent.hh"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>

namespace CellStateCalibration {

RepairMediatedMultiComponentCalibrator::RepairMediatedMultiComponentCalibrator(
    const FitConfigS2& config)
    : m_config(config)
    , m_model(config)
{}

void RepairMediatedMultiComponentCalibrator::setConfig(const FitConfigS2& config) {
    m_config = config;
    m_model.setConfig(config);
}

void RepairMediatedMultiComponentCalibrator::setData(const std::vector<DataPoint>& data) {
    m_data = data;
}

void RepairMediatedMultiComponentCalibrator::addDataPoint(const DataPoint& point) {
    m_data.push_back(point);
}

void RepairMediatedMultiComponentCalibrator::clearData() {
    m_data.clear();
}

void RepairMediatedMultiComponentCalibrator::setDefaultUncertainties(double rel_err) {
    for (auto& dp : m_data) {
        dp.si = rel_err;
    }
}

double RepairMediatedMultiComponentCalibrator::negativeLogLikelihood(
    const ParamsS2& params) const
{
    double nll = 0.0;

    for (const auto& dp : m_data) {
        if (dp.D > m_config.dose_max_fit) continue;

        double sf_model = m_model.survivalFraction(dp.D, params);

        double y = std::log(clamp(dp.SF_obs, m_config.eps_sf, 1.0));
        double mu = std::log(sf_model);
        double s = std::max(dp.si, 1e-6);

        double residual = (y - mu) / s;
        nll += 0.5 * residual * residual + std::log(s);
    }

    return nll;
}

static double vec_dist2_mc(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return s;
}

std::vector<double> RepairMediatedMultiComponentCalibrator::nelderMeadND(
    std::function<double(const std::vector<double>&)> f,
    const std::vector<double>& x0,
    const std::vector<double>& step,
    const std::vector<double>& lo,
    const std::vector<double>& hi,
    double& best_f_out,
    int& iterations_out,
    bool& converged_out) const
{
    const int max_iter = 10000;
    const double tol_f = 1e-10;
    const double tol_x = 1e-8;
    const double alpha_nm = 1.0;
    const double gamma_nm = 2.0;
    const double rho_nm = 0.5;
    const double sigma_nm = 0.5;

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
            xspan = std::max(xspan, std::sqrt(vec_dist2_mc(x[i], x[0])));
        }

        if (fspan < tol_f && xspan < tol_x) {
            converged_out = true;
            break;
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
                xe[d] = xc[d] + gamma_nm * (xr[d] - xc[d]);
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
            x[m-1] = xr;
            fx[m-1] = fr;
        } else {
            std::vector<double> xk(n);
            if (fr < fx[m-1]) {
                for (int d = 0; d < n; ++d) {
                    xk[d] = xc[d] + rho_nm * (xr[d] - xc[d]);
                }
            } else {
                for (int d = 0; d < n; ++d) {
                    xk[d] = xc[d] - rho_nm * (xc[d] - x[m-1][d]);
                }
            }
            xk = bounded(xk);
            double fk = f(xk);

            if (fk < fx[m-1]) {
                x[m-1] = xk;
                fx[m-1] = fk;
            } else {
                for (int i = 1; i < m; ++i) {
                    for (int d = 0; d < n; ++d) {
                        x[i][d] = x[0][d] + sigma_nm * (x[i][d] - x[0][d]);
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

MultiComponentFitResult RepairMediatedMultiComponentCalibrator::fit(
    const ParamsS2& initial_guess,
    const ParamsS2& step)
{
    MultiComponentFitResult result;

    const bool fix_T21 = m_config.fix_T21;
    const bool fix_T23 = m_config.fix_T23;
    const bool fix_Nc  = m_config.fix_Nc;
    const double T21_fixed = m_config.T21_fixed;
    const double T23_fixed = m_config.T23_fixed;
    const double Nc_fixed  = m_config.Nc_fixed;

    double best_nll = 0.0;
    int iterations = 0;
    bool converged = false;

    // Build the free-parameter vector dynamically based on what's fixed.
    // Base free params: e2, e3, a, k_error (always optimized)
    // Optional: T21, T23 (if not fixed)
    // Multi-component: Nc (if not fixed), omega_p, lambda_p (always optimized)

    // We pack: [e2, e3, a, k_error, <T21?>, <T23?>, <Nc?>, omega_p, lambda_p]
    // and unpack back to ParamsS2 in the objective.

    std::vector<double> x0_vec, step_vec, lo_vec, hi_vec;

    // Indices for unpacking
    enum { IDX_E2 = 0, IDX_E3 = 1, IDX_A = 2, IDX_KERR = 3 };
    x0_vec.push_back(initial_guess.e2);    step_vec.push_back(step.e2);
    lo_vec.push_back(0.1);                 hi_vec.push_back(80.0);

    x0_vec.push_back(initial_guess.e3);    step_vec.push_back(step.e3);
    lo_vec.push_back(0.2);                 hi_vec.push_back(200.0);

    x0_vec.push_back(initial_guess.a);     step_vec.push_back(step.a);
    lo_vec.push_back(0.001);               hi_vec.push_back(50.0);

    x0_vec.push_back(initial_guess.k_error); step_vec.push_back(step.k_error);
    lo_vec.push_back(1e-10);               hi_vec.push_back(1e-2);

    int idx_T21 = -1, idx_T23 = -1, idx_Nc = -1, idx_omega = -1, idx_lambda = -1;

    if (!fix_T21) {
        idx_T21 = static_cast<int>(x0_vec.size());
        x0_vec.push_back(initial_guess.T21);  step_vec.push_back(step.T21);
        lo_vec.push_back(1.0);                 hi_vec.push_back(100.0);
    }
    if (!fix_T23) {
        idx_T23 = static_cast<int>(x0_vec.size());
        x0_vec.push_back(initial_guess.T23);  step_vec.push_back(step.T23);
        lo_vec.push_back(1.0);                 hi_vec.push_back(100.0);
    }
    if (!fix_Nc) {
        idx_Nc = static_cast<int>(x0_vec.size());
        x0_vec.push_back(initial_guess.Nc);   step_vec.push_back(step.Nc);
        lo_vec.push_back(5.0);                 hi_vec.push_back(500.0);
    }

    idx_omega = static_cast<int>(x0_vec.size());
    x0_vec.push_back(initial_guess.omega_p);   step_vec.push_back(step.omega_p);
    lo_vec.push_back(0.5);                     hi_vec.push_back(5.0);

    idx_lambda = static_cast<int>(x0_vec.size());
    x0_vec.push_back(initial_guess.lambda_p);  step_vec.push_back(step.lambda_p);
    lo_vec.push_back(0.005);                   hi_vec.push_back(0.5);

    auto objective = [this, fix_T21, fix_T23, fix_Nc,
                      T21_fixed, T23_fixed, Nc_fixed,
                      idx_T21, idx_T23, idx_Nc, idx_omega, idx_lambda]
        (const std::vector<double>& v) -> double
    {
        ParamsS2 p;
        p.e2      = v[IDX_E2];
        p.e3      = v[IDX_E3];
        p.a       = v[IDX_A];
        p.k_error = v[IDX_KERR];

        p.T21 = fix_T21 ? T21_fixed : v[idx_T21];
        p.T23 = fix_T23 ? T23_fixed : v[idx_T23];
        p.Nc  = fix_Nc  ? Nc_fixed  : v[idx_Nc];

        p.omega_p  = v[idx_omega];
        p.lambda_p = v[idx_lambda];

        if (!p.isValid()) return 1e30;
        return negativeLogLikelihood(p);
    };

    std::vector<double> xhat = nelderMeadND(objective, x0_vec, step_vec, lo_vec, hi_vec,
                                             best_nll, iterations, converged);

    // Unpack result
    result.params.e2      = xhat[IDX_E2];
    result.params.e3      = xhat[IDX_E3];
    result.params.a       = xhat[IDX_A];
    result.params.k_error = xhat[IDX_KERR];

    result.params.T21 = fix_T21 ? T21_fixed : xhat[idx_T21];
    result.params.T23 = fix_T23 ? T23_fixed : xhat[idx_T23];
    result.params.Nc  = fix_Nc  ? Nc_fixed  : xhat[idx_Nc];

    result.params.omega_p  = xhat[idx_omega];
    result.params.lambda_p = xhat[idx_lambda];

    result.negLogLikelihood = best_nll;
    result.iterations = iterations;
    result.converged = converged;

    return result;
}

std::vector<double> RepairMediatedMultiComponentCalibrator::predictSurvival(
    const ParamsS2& params,
    const std::vector<double>& doses) const
{
    return m_model.survivalCurve(doses, params);
}

const std::vector<DataPoint>& RepairMediatedMultiComponentCalibrator::getData() const {
    return m_data;
}

const RepairMediatedMultiComponentModel& RepairMediatedMultiComponentCalibrator::getModel() const {
    return m_model;
}

} // namespace CellStateCalibration
