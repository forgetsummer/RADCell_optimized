/**
 * @file DirectLethalCalibrator.cc
 * @brief Implementation of Direct Lethal Model calibrator
 */

#include "DirectLethalCalibrator.hh"
#include "MathUtilities.hh"
#include <cmath>
#include <algorithm>

namespace CellStateCalibration {

DirectLethalCalibrator::DirectLethalCalibrator(const FitConfig& config)
    : m_config(config)
    , m_model(config)
{
}

void DirectLethalCalibrator::setConfig(const FitConfig& config) {
    m_config = config;
    m_model.setConfig(config);
}

void DirectLethalCalibrator::setData(const std::vector<DataPoint>& data) {
    m_data = data;
}

void DirectLethalCalibrator::addDataPoint(const DataPoint& point) {
    m_data.push_back(point);
}

void DirectLethalCalibrator::clearData() {
    m_data.clear();
}

void DirectLethalCalibrator::setDefaultUncertainties(double rel_err) {
    // If SF has ~rel_err relative error, then std dev of ln(SF) ~ rel_err
    // because Var(ln X) ~ Var(X)/E[X]^2 for small noise
    for (auto& dp : m_data) {
        dp.si = rel_err;
    }
}

double DirectLethalCalibrator::negativeLogLikelihood(const ParamsReduced& params) const {
    // Gaussian likelihood on log(SF)
    // NLL = Σ 0.5 * ((ln SF_obs - ln SF_model) / si)^2 + ln(si)
    
    double nll = 0.0;
    
    for (const auto& dp : m_data) {
        // Skip data points beyond dose_max_fit
        if (dp.D > m_config.dose_max_fit) continue;
        
        double sf_model = m_model.survivalFraction(dp.D, params);
        double y = std::log(clamp(dp.SF_obs, m_config.eps_sf, 1.0));
        double mu = std::log(sf_model);
        double s = std::max(dp.si, 1e-6);
        double r = (y - mu) / s;
        
        nll += 0.5 * r * r + std::log(s);
    }
    
    return nll;
}

DirectLethalFitResult DirectLethalCalibrator::fit(const ParamsReduced& initial_guess,
                                                   const ParamsReduced& step) {
    DirectLethalFitResult result;
    
    // Create objective function
    auto objective = [this](const ParamsReduced& p) -> double {
        // Basic sanity check
        if (!(p.e3 > 0.0 && p.a > 0.0)) return 1e30;
        return negativeLogLikelihood(p);
    };
    
    // Configure optimizer
    NelderMeadOptions opt;
    opt.max_iter = 3000;
    opt.tol_f = 1e-10;
    opt.tol_x = 1e-8;
    m_optimizer.setOptions(opt);
    m_optimizer.setBounds(ParamsReduced(0.0, 0.0), ParamsReduced(200.0, 200.0));
    
    // Run optimization
    result.params = m_optimizer.optimize(objective, initial_guess, step, 
                                          result.negLogLikelihood);
    result.iterations = m_optimizer.getIterations();
    
    // Compute Hessian at optimum
    result.hessian = computeHessian(result.params);
    
    // Check convergence (Hessian should be positive definite)
    result.converged = result.hessian.isPositiveDefinite();
    
    return result;
}

Hessian2 DirectLethalCalibrator::computeHessian(const ParamsReduced& params,
                                                 double h_e3, double h_a) const {
    // Use automatic step sizes if not specified
    if (h_e3 < 0) h_e3 = 1e-3 * std::max(1.0, params.e3);
    if (h_a < 0) h_a = 1e-3 * std::max(1.0, params.a);
    
    // Create evaluation function
    auto eval = [this, &params](double de3, double da) -> double {
        ParamsReduced y(
            clamp(params.e3 + de3, 0.0, 200.0),
            clamp(params.a + da, 0.0, 200.0)
        );
        return negativeLogLikelihood(y);
    };
    
    // Central differences
    double f00 = eval(0, 0);
    double fpp = eval(+h_e3, +h_a);
    double fpm = eval(+h_e3, -h_a);
    double fmp = eval(-h_e3, +h_a);
    double fmm = eval(-h_e3, -h_a);
    double fxp = eval(+h_e3, 0);
    double fxm = eval(-h_e3, 0);
    double fyp = eval(0, +h_a);
    double fym = eval(0, -h_a);
    
    Hessian2 H;
    H.h11 = (fxp - 2*f00 + fxm) / (h_e3 * h_e3);
    H.h22 = (fyp - 2*f00 + fym) / (h_a * h_a);
    H.h12 = (fpp - fpm - fmp + fmm) / (4 * h_e3 * h_a);
    
    return H;
}

std::vector<std::pair<double, double>> DirectLethalCalibrator::profileLikelihood_e3(
    const ParamsReduced& params_opt,
    double e3_min,
    double e3_max,
    int n_grid) const
{
    std::vector<std::pair<double, double>> profile;
    profile.reserve(n_grid);
    
    for (int i = 0; i < n_grid; ++i) {
        double e3 = e3_min + (e3_max - e3_min) * static_cast<double>(i) / (n_grid - 1);
        
        // 1D optimization over 'a' with e3 fixed
        double best_nll = std::numeric_limits<double>::infinity();
        
        // Coarse grid search over 'a'
        double a_min = std::max(0.01, params_opt.a * 0.2);
        double a_max = std::min(200.0, params_opt.a * 3.0);
        
        for (int k = 0; k < 200; ++k) {
            double a = a_min + (a_max - a_min) * static_cast<double>(k) / 199.0;
            ParamsReduced p(e3, a);
            double nll = negativeLogLikelihood(p);
            best_nll = std::min(best_nll, nll);
        }
        
        profile.push_back({e3, best_nll});
    }
    
    return profile;
}

std::vector<std::pair<double, double>> DirectLethalCalibrator::profileLikelihood_a(
    const ParamsReduced& params_opt,
    double a_min,
    double a_max,
    int n_grid) const
{
    std::vector<std::pair<double, double>> profile;
    profile.reserve(n_grid);
    
    for (int i = 0; i < n_grid; ++i) {
        double a = a_min + (a_max - a_min) * static_cast<double>(i) / (n_grid - 1);
        
        // 1D optimization over 'e3' with a fixed
        double best_nll = std::numeric_limits<double>::infinity();
        
        // Coarse grid search over 'e3'
        double e3_min_search = std::max(0.01, params_opt.e3 * 0.2);
        double e3_max_search = std::min(200.0, params_opt.e3 * 3.0);
        
        for (int k = 0; k < 200; ++k) {
            double e3 = e3_min_search + (e3_max_search - e3_min_search) * static_cast<double>(k) / 199.0;
            ParamsReduced p(e3, a);
            double nll = negativeLogLikelihood(p);
            best_nll = std::min(best_nll, nll);
        }
        
        profile.push_back({a, best_nll});
    }
    
    return profile;
}

std::vector<double> DirectLethalCalibrator::predictSurvival(const ParamsReduced& params,
                                                             const std::vector<double>& doses) const {
    return m_model.survivalCurve(doses, params);
}

const std::vector<DataPoint>& DirectLethalCalibrator::getData() const {
    return m_data;
}

const DirectLethalModel& DirectLethalCalibrator::getModel() const {
    return m_model;
}

} // namespace CellStateCalibration
