/**
 * @file RepairMediatedModel_multicomponent.cc
 * @brief Implementation of multi-component energy survival model
 *
 * This model mirrors the MC simulation's step-by-step logic:
 * - Energy decays each time step (bi-exponential for Er, mono-exponential for Ep)
 * - Transition probabilities are re-evaluated every step from current energy
 * - S1 uses total energy E = Er + Ep for overlaps
 * - S2 uses effective energy E_eff = Er + omega_p * Ep for overlaps
 * - Misrepair uses (k_error/alpha) * Er, conditional on recovery attempt
 * - First step (DSB injection): instantaneous probability (overlap score directly)
 * - Subsequent steps: delayed probability 1 - (1-p)^(dt/To)
 */

#include "RepairMediatedModel_multicomponent.hh"
#include <cmath>
#include <algorithm>

namespace CellStateCalibration {

RepairMediatedMultiComponentModel::RepairMediatedMultiComponentModel(const FitConfigS2& config)
    : m_config(config)
{}

void RepairMediatedMultiComponentModel::setConfig(const FitConfigS2& config) {
    m_config = config;
}

const FitConfigS2& RepairMediatedMultiComponentModel::getConfig() const {
    return m_config;
}

double RepairMediatedMultiComponentModel::overlapScore(double e, double ej) const {
    double z = -std::abs(e - ej) / 2.0;
    double val = 2.0 * normal_cdf(z);
    return clamp(val, 0.0, 1.0);
}

void RepairMediatedMultiComponentModel::assignmentProbabilities(
    int n, const ParamsS2& params,
    double& P1, double& P2, double& P3) const
{
    const double e = params.a * static_cast<double>(n);

    const double w1 = overlapScore(e, 0.0);
    const double w2 = overlapScore(e, params.e2);
    const double w3 = overlapScore(e, params.e3);

    const double sum = w1 + w2 + w3;

    if (sum <= 1e-15) {
        double d1 = std::abs(e - 0.0);
        double d2 = std::abs(e - params.e2);
        double d3 = std::abs(e - params.e3);

        P1 = (d1 <= d2 && d1 <= d3) ? 1.0 : 0.0;
        P2 = (d2 < d1 && d2 <= d3) ? 1.0 : 0.0;
        P3 = (d3 < d1 && d3 < d2) ? 1.0 : 0.0;
        return;
    }

    P1 = w1 / sum;
    P2 = w2 / sum;
    P3 = w3 / sum;
}

double RepairMediatedMultiComponentModel::cumulativeProbToRate(double p_cum, double T0) const {
    p_cum = clamp(p_cum, 0.0, 1.0 - 1e-15);
    T0 = std::max(T0, 1e-9);
    return -std::log(1.0 - p_cum) / T0;
}

double RepairMediatedMultiComponentModel::Er(
    double t, int /*n*/, double alpha_n_1mh) const
{
    const double f1 = m_config.f1;
    const double f2 = m_config.f2;
    const double l1 = m_config.lambda1_rep;
    const double l2 = m_config.lambda2_rep;
    return alpha_n_1mh * (f1 * std::exp(-l1 * t) + f2 * std::exp(-l2 * t));
}

double RepairMediatedMultiComponentModel::Ep(
    double t, double alpha_n_h, double lambda_p) const
{
    return alpha_n_h * std::exp(-lambda_p * t);
}

double RepairMediatedMultiComponentModel::computeS2RepairProbability(
    int n, const ParamsS2& params) const
{
    const double nd = static_cast<double>(n);
    const double Nc = std::max(params.Nc, 1e-9);
    const double h_n = nd / (nd + Nc);

    const double alpha_n = params.a * nd;
    const double alpha_n_1mh = alpha_n * (1.0 - h_n);
    const double alpha_n_h   = alpha_n * h_n;

    const double omega_p  = params.omega_p;
    const double lambda_p = params.lambda_p;
    const double T_assay  = m_config.T_assay_h;

    // Misrepair: (k_error / a) * Er(t) — matches simulation's (k_error/alpha) * Er
    // In reduced space, alpha = a*sigma, Er_phys = Er_reduced * sigma
    // so (k_error/alpha) * Er_phys = (k_error/(a*sigma)) * (Er_reduced*sigma) = (k_error/a) * Er_reduced
    // But Er_reduced here is alpha_n_1mh * decay, where alpha_n_1mh = a * n * (1-h)
    // So (k_error/a) * Er_reduced = k_error * n * (1-h) * decay
    const double k_error_over_a = (params.a > 1e-15) ? (params.k_error / params.a) : 0.0;

    const double dt = T_assay / static_cast<double>(NUM_INTEGRATION_STEPS);

    double cumulative_hazard = 0.0;
    double R21 = 0.0;

    double e_eff_prev = Er(0.0, n, alpha_n_1mh) + omega_p * Ep(0.0, alpha_n_h, lambda_p);
    double p21_cum_prev = overlapScore(e_eff_prev, 0.0);
    double p23_cum_prev = overlapScore(e_eff_prev, params.e3);
    double lam21_prev = cumulativeProbToRate(p21_cum_prev, params.T21);

    // Misrepair: lambda_mis = (k_error/a) * Er(t)
    double Er_prev = Er(0.0, n, alpha_n_1mh);
    double lambda_mis_prev = k_error_over_a * std::max(0.0, Er_prev);

    // In simulation, p23_mis = p21 * p_mis, p21_eff = p21 * (1-p_mis)
    // Continuous-time equivalent: the recovery rate is split into
    //   effective_lam21 = lam21 * (1 - p_mis_rate)
    //   lam23_mis = lam21 * p_mis_rate
    // where p_mis_rate represents the fraction of recovery attempts that fail
    // For small dt: p_mis ≈ lambda_mis * dt, so in continuous time the misrepair
    // acts as an additive hazard that diverts recovery into death.
    // Total lambda23 = lambda23_energy + lambda_mis (continuous-time competing risk)
    double lam23_prev = cumulativeProbToRate(p23_cum_prev, params.T23) + lambda_mis_prev;
    double total_prev = lam21_prev + lam23_prev;
    double S2_prev = 1.0;
    double integrand_prev = lam21_prev * S2_prev;

    for (int i = 1; i <= NUM_INTEGRATION_STEPS; ++i) {
        const double t = dt * static_cast<double>(i);

        const double e_r = Er(t, n, alpha_n_1mh);
        const double e_p = Ep(t, alpha_n_h, lambda_p);
        const double e_eff = e_r + omega_p * e_p;

        const double p21_cum = overlapScore(e_eff, 0.0);
        const double p23_cum = overlapScore(e_eff, params.e3);

        const double lam21 = cumulativeProbToRate(p21_cum, params.T21);
        const double lambda_mis = k_error_over_a * std::max(0.0, e_r);
        const double lam23 = cumulativeProbToRate(p23_cum, params.T23) + lambda_mis;
        const double total = lam21 + lam23;

        cumulative_hazard += 0.5 * (total_prev + total) * dt;
        const double S2_cur = std::exp(-cumulative_hazard);

        const double integrand_cur = lam21 * S2_cur;

        R21 += 0.5 * (integrand_prev + integrand_cur) * dt;

        total_prev = total;
        integrand_prev = integrand_cur;
    }

    return clamp(R21, 0.0, 1.0);
}

double RepairMediatedMultiComponentModel::survivalFraction(
    double dose, const ParamsS2& params) const
{
    const double lambda = m_config.kappa * dose;
    const int nmax = poisson_nmax(lambda);

    const double T_assay = m_config.T_assay_h;
    // Use ~6 minute steps: fine enough to capture energy decay dynamics
    // while keeping calibration tractable (200 steps over 20h assay)
    const int nsteps = NUM_INTEGRATION_STEPS;
    const double dt_step = T_assay / static_cast<double>(nsteps);

    const double f1 = m_config.f1;
    const double f2 = m_config.f2;
    const double l1 = m_config.lambda1_rep;
    const double l2 = m_config.lambda2_rep;
    const double lambda_p = params.lambda_p;
    const double omega_p = params.omega_p;

    // Bi-exponential decay factor per step
    const double decay_Er = f1 * std::exp(-l1 * dt_step) + f2 * std::exp(-l2 * dt_step);
    const double decay_Ep = std::exp(-lambda_p * dt_step);

    // S1 observation timescales (from config, typically cell cycle time)
    const double To12 = m_config.To12;
    const double To13 = m_config.To13;

    // Misrepair rate conversion (k_error/a in reduced space, matching (k_error/alpha) in physical)
    const double k_error_over_a = (params.a > 1e-15) ? (params.k_error / params.a) : 0.0;

    double sf = 0.0;

    for (int n = 0; n <= nmax; ++n) {
        const double pn = std::exp(poisson_log_pmf(n, lambda));
        const double nd = static_cast<double>(n);
        const double Nc = std::max(params.Nc, 1e-9);
        const double h_n = nd / (nd + Nc);

        // Initial energy injection (reduced space)
        const double E_inj = params.a * nd;
        const double dEr_inj = E_inj * (1.0 - h_n);
        const double dEp_inj = E_inj * h_n;

        // Track population fractions in three states
        double frac_S1 = 1.0;
        double frac_S2 = 0.0;
        // frac_S3 = 1 - frac_S1 - frac_S2 (absorbing)

        // Track energy (starts at 0, injection happens at step 1)
        double cur_Er = 0.0;
        double cur_Ep = 0.0;

        for (int step = 0; step < nsteps; ++step) {
            // Energy update: decay then inject (injection only on first step)
            cur_Er = cur_Er * decay_Er;
            cur_Ep = cur_Ep * decay_Ep;
            if (step == 0) {
                cur_Er += dEr_inj;
                cur_Ep += dEp_inj;
            }

            const double E_total = cur_Er + cur_Ep;  // S1 uses total energy
            const double E_eff = cur_Er + omega_p * cur_Ep;  // S2 uses effective energy

            // Is this the irradiation step? (determines instantaneous vs delayed)
            const bool is_irradiation_step = (step == 0 && n > 0);

            // === S1 transitions ===
            if (frac_S1 > 1e-15) {
                double p12_overlap = overlapScore(E_total, params.e2);
                double p13_overlap = overlapScore(E_total, params.e3);

                double p12_step, p13_step;
                if (is_irradiation_step) {
                    p12_step = p12_overlap;
                    p13_step = p13_overlap;
                } else {
                    // Delayed: 1 - (1-p)^(dt/To) — S1 uses To12/To13 (cell cycle timescale)
                    p12_step = 1.0 - std::pow(std::max(0.0, 1.0 - p12_overlap), dt_step / To12);
                    p13_step = 1.0 - std::pow(std::max(0.0, 1.0 - p13_overlap), dt_step / To13);
                }

                p12_step = clamp(p12_step, 0.0, 1.0);
                p13_step = clamp(p13_step, 0.0, 1.0);

                // Normalize if sum > 1
                const double pSum12_13 = p12_step + p13_step;
                if (pSum12_13 > 1.0) {
                    p12_step /= pSum12_13;
                    p13_step /= pSum12_13;
                }

                // Population transfers from S1
                const double to_S2 = frac_S1 * p12_step;
                const double to_S3 = frac_S1 * p13_step;
                frac_S1 -= (to_S2 + to_S3);
                frac_S2 += to_S2;
                // to_S3 goes to absorbing state
            }

            // === S2 transitions ===
            if (frac_S2 > 1e-15) {
                double p21_overlap = overlapScore(E_eff, 0.0);
                double p23_overlap = overlapScore(E_eff, params.e3);

                double p21_step, p23_0_step;
                if (is_irradiation_step) {
                    p21_step = p21_overlap;
                    p23_0_step = p23_overlap;
                } else {
                    p21_step = 1.0 - std::pow(std::max(0.0, 1.0 - p21_overlap), dt_step / params.T21);
                    p23_0_step = 1.0 - std::pow(std::max(0.0, 1.0 - p23_overlap), dt_step / params.T23);
                }

                p21_step = clamp(p21_step, 0.0, 1.0);
                p23_0_step = clamp(p23_0_step, 0.0, 1.0);

                // Misrepair: lambda_mis = (k_error/a) * Er, p_mis = 1 - exp(-lambda_mis * dt)
                const double lambda_mis = k_error_over_a * std::max(0.0, cur_Er);
                double p_mis = 1.0 - std::exp(-lambda_mis * dt_step);
                p_mis = clamp(p_mis, 0.0, 1.0);

                // Misrepair diverts recovery attempts into death
                const double p23_mis_step = p21_step * p_mis;
                const double p21_eff_step = p21_step * (1.0 - p_mis);
                double p23_step = p23_0_step + p23_mis_step;
                p23_step = clamp(p23_step, 0.0, 1.0);

                // Normalize if needed
                const double pSum = p21_eff_step + p23_step;
                double p21_final = p21_eff_step;
                double p23_final = p23_step;
                if (pSum > 1.0) {
                    p21_final /= pSum;
                    p23_final /= pSum;
                }

                // Population transfers from S2
                const double to_S1 = frac_S2 * p21_final;
                const double to_S3 = frac_S2 * p23_final;
                frac_S2 -= (to_S1 + to_S3);
                frac_S1 += to_S1;
            }

            // Safety clamp
            frac_S1 = clamp(frac_S1, 0.0, 1.0);
            frac_S2 = clamp(frac_S2, 0.0, 1.0);
        }

        // Survival = S1 population at end of assay
        sf += pn * clamp(frac_S1, 0.0, 1.0);
    }

    return clamp(sf, m_config.eps_sf, 1.0);
}

std::vector<double> RepairMediatedMultiComponentModel::survivalCurve(
    const std::vector<double>& doses, const ParamsS2& params) const
{
    std::vector<double> sf_values;
    sf_values.reserve(doses.size());
    for (double dose : doses) {
        sf_values.push_back(survivalFraction(dose, params));
    }
    return sf_values;
}

} // namespace CellStateCalibration
