#ifndef DATA_TYPES_HH
#define DATA_TYPES_HH

/**
 * @file DataTypes.hh
 * @brief Data structures for cell survival curve fitting
 * 
 * Defines the core data types used in the calibration process:
 * - DataPoint: Individual dose-survival measurement
 * - FitConfig: Configuration parameters for fitting
 * - ParamsReduced: Reduced (identifiable) model parameters
 * - Hessian2: 2x2 Hessian matrix for uncertainty analysis
 */

#include <vector>
#include <iostream>
#include <iomanip>

namespace CellStateCalibration {

/**
 * @brief Single experimental data point
 */
struct DataPoint {
    double D;        ///< Dose in Gy
    double SF_obs;   ///< Observed survival fraction
    double si;       ///< Log-space standard deviation (on ln SF)
    
    DataPoint() : D(0.0), SF_obs(1.0), si(0.2) {}
    DataPoint(double dose, double sf, double sigma = 0.2) 
        : D(dose), SF_obs(sf), si(sigma) {}
};

/**
 * @brief Configuration for the fitting procedure
 */
struct FitConfig {
    double kappa;        ///< DSB/Gy (double-strand breaks per Gray)
    double eps_sf;       ///< Floor value to avoid log(0)
    double dose_max_fit; ///< Maximum dose to include in fit (Gy)
    
    FitConfig() 
        : kappa(40.0)
        , eps_sf(1e-12)
        , dose_max_fit(4.0) 
    {}
};

/**
 * @brief Reduced (identifiable) model parameters
 * 
 * The full model has parameters (E1, E3, alpha, sigma) but only
 * two combinations are identifiable from survival data:
 * - e3 = E3/sigma (dimensionless threshold)
 * - a = alpha/sigma (dimensionless damage sensitivity)
 * 
 * Note: E1 is set to 0 (reference point)
 */
struct ParamsReduced {
    double e3;  ///< Dimensionless threshold E3/sigma
    double a;   ///< Dimensionless damage sensitivity alpha/sigma
    
    ParamsReduced() : e3(20.0), a(1.0) {}
    ParamsReduced(double e3_val, double a_val) : e3(e3_val), a(a_val) {}
    
    /**
     * @brief Get critical DSB count (independent of sigma)
     * @return Number of DSBs at threshold = E3/alpha = e3/a
     */
    double criticalDSBCount() const {
        return (a > 0) ? e3 / a : 0.0;
    }
};

/**
 * @brief Full physical parameters for the cell state model
 * 
 * These are the actual parameters needed by the CellStateModel:
 * - E1: Healthy state energy (reference, typically 0)
 * - E3: Death state threshold
 * - sigma: Gaussian width for state transitions
 * - alpha: Damage per DSB
 */
struct CellStateModelParams {
    double E1;      ///< Healthy state energy (typically 0)
    double E3;      ///< Death state threshold
    double sigma;   ///< Gaussian width for transitions
    double alpha;   ///< Damage per DSB
    double kappa;   ///< DSB per Gy (for reference)
    
    CellStateModelParams() 
        : E1(0.0), E3(0.0), sigma(1.0), alpha(0.0), kappa(40.0) {}
    
    CellStateModelParams(double e1, double e3, double sig, double alph, double kap = 40.0)
        : E1(e1), E3(e3), sigma(sig), alpha(alph), kappa(kap) {}
    
    /**
     * @brief Create from reduced parameters by specifying sigma
     * @param reduced Reduced parameters from calibration
     * @param sigma_value Value of sigma to use
     * @param kappa_value DSB/Gy value
     * @return Full physical parameters
     */
    static CellStateModelParams fromReduced(const ParamsReduced& reduced, 
                                             double sigma_value,
                                             double kappa_value = 40.0) {
        CellStateModelParams params;
        params.E1 = 0.0;  // Convention
        params.E3 = reduced.e3 * sigma_value;
        params.sigma = sigma_value;
        params.alpha = reduced.a * sigma_value;
        params.kappa = kappa_value;
        return params;
    }
    
    /**
     * @brief Create from reduced parameters using the constraint sigma = sqrt(E3)
     * 
     * This matches the constraint used in the existing CellStateModel:
     *   sigma = sqrt(E3)
     *   E3 = 4 * x^2  where x = Phi^-1(p_sp/2)
     * 
     * With this constraint:
     *   e3 = E3/sigma = E3/sqrt(E3) = sqrt(E3)
     *   So: E3 = e3^2 and sigma = e3
     * 
     * @param reduced Reduced parameters from calibration
     * @param kappa_value DSB/Gy value
     * @return Full physical parameters with sigma = sqrt(E3) constraint
     */
    static CellStateModelParams fromReducedWithSigmaConstraint(
        const ParamsReduced& reduced,
        double kappa_value = 40.0) {
        CellStateModelParams params;
        params.E1 = 0.0;
        // With sigma = sqrt(E3), we have e3 = sqrt(E3), so E3 = e3^2
        params.E3 = reduced.e3 * reduced.e3;
        params.sigma = reduced.e3;  // sigma = sqrt(E3) = e3
        // alpha = a * sigma = a * e3
        params.alpha = reduced.a * reduced.e3;
        params.kappa = kappa_value;
        return params;
    }
    
    /**
     * @brief Get critical DSB count
     */
    double criticalDSBCount() const {
        return (alpha > 0) ? E3 / alpha : 0.0;
    }
    
    /**
     * @brief Get dose at which average damage equals threshold
     */
    double criticalDose() const {
        return criticalDSBCount() / kappa;
    }
    
    /**
     * @brief Print parameters
     */
    void print() const;
};

/**
 * @brief 2x2 Hessian matrix for parameter uncertainty
 */
struct Hessian2 {
    double h11;  ///< d²f/de3²
    double h12;  ///< d²f/de3da
    double h22;  ///< d²f/da²
    
    Hessian2() : h11(0.0), h12(0.0), h22(0.0) {}
    
    /**
     * @brief Compute determinant
     */
    double determinant() const {
        return h11 * h22 - h12 * h12;
    }
    
    /**
     * @brief Check if matrix is positive definite
     */
    bool isPositiveDefinite() const {
        return h11 > 0 && determinant() > 0;
    }
};

/**
 * @brief Options for Nelder-Mead optimizer
 */
struct NelderMeadOptions {
    int max_iter;      ///< Maximum iterations
    double tol_f;      ///< Tolerance on function value
    double tol_x;      ///< Tolerance on parameter values
    double alpha;      ///< Reflection coefficient
    double gamma;      ///< Expansion coefficient
    double rho;        ///< Contraction coefficient
    double sigma;      ///< Shrink coefficient
    
    NelderMeadOptions()
        : max_iter(2000)
        , tol_f(1e-8)
        , tol_x(1e-6)
        , alpha(1.0)
        , gamma(2.0)
        , rho(0.5)
        , sigma(0.5)
    {}
};

/**
 * @brief Extended fit configuration for Repair Mediated model
 *
 * Includes additional parameters for the S2 (repair) state kinetics.
 *
 * T21 and T23 can be either:
 * 1. FIXED by user (fix_T21=true, fix_T23=true): Use T21_fixed and T23_fixed values
 * 2. OPTIMIZED during fitting (fix_T21=false, fix_T23=false): Fitted along with other params
 *
 * When fixed, they represent measurable biological timescales:
 * - T21: Repair timescale (S2->S1), typically related to DNA repair half-life
 * - T23: Death timescale (S2->S3), typically related to cell cycle
 */
struct FitConfigS2 {
    double kappa;        ///< DSB/Gy (double-strand breaks per Gray)
    double eps_sf;       ///< Floor value to avoid log(0)
    double dose_max_fit; ///< Maximum dose to include in fit (Gy)
    double T_assay_h;    ///< Assay time window in hours (e.g., colony formation ~240h)
    double T21_fixed;    ///< User-defined repair timescale in hours (S2->S1)
    double T23_fixed;    ///< User-defined death timescale in hours (S2->S3)
    bool fix_T21;        ///< If true, T21 is fixed at T21_fixed; if false, T21 is optimized
    bool fix_T23;        ///< If true, T23 is fixed at T23_fixed; if false, T23 is optimized

    FitConfigS2()
        : kappa(40.0)
        , eps_sf(1e-12)
        , dose_max_fit(6.0)
        , T_assay_h(24.0 * 10.0)  // 10 days default for clonogenic assay
        , T21_fixed(12.0)         // Default: 12 hours
        , T23_fixed(24.0)         // Default: 24 hours
        , fix_T21(false)          // Default: optimize T21
        , fix_T23(false)          // Default: optimize T23
    {}

    /**
     * @brief Set timescales to cell cycle time and fix them
     * @param T_cellCycle Cell cycle time in hours
     */
    void setTimescalesToCellCycle(double T_cellCycle) {
        T21_fixed = T_cellCycle;
        T23_fixed = T_cellCycle;
        fix_T21 = true;
        fix_T23 = true;
    }

    /**
     * @brief Fix T21 and T23 to specific user-defined values
     * @param T21_val Repair timescale in hours
     * @param T23_val Death timescale in hours
     */
    void fixTimescales(double T21_val, double T23_val) {
        T21_fixed = T21_val;
        T23_fixed = T23_val;
        fix_T21 = true;
        fix_T23 = true;
    }

    /**
     * @brief Allow T21 and T23 to be optimized during fitting
     */
    void optimizeTimescales() {
        fix_T21 = false;
        fix_T23 = false;
    }

    /**
     * @brief Check if any timescales are fixed
     */
    bool hasFixedTimescales() const {
        return fix_T21 || fix_T23;
    }
};

/**
 * @brief Reduced parameters for Repair Mediated model (S1->S2->S3)
 * 
 * FIVE parameters to be optimized (with E1=0 convention):
 * - e2 = E2/sigma: Repair state threshold (dimensionless)
 * - e3 = E3/sigma: Death state threshold (dimensionless)
 * - a = alpha/sigma: Damage per DSB (dimensionless)
 * - T21: Repair timescale in hours (S2->S1)
 * - T23: Death timescale in hours (S2->S3)
 * 
 * Physical interpretation:
 * - S1 (E=0): Healthy state
 * - S2 (E=e2*sigma): Arrested/repair state
 * - S3 (E=e3*sigma): Dead state
 * - After damage, cells can:
 *   1. Stay healthy (S1) if damage is low
 *   2. Enter repair (S2) if damage is moderate
 *   3. Die immediately (S3) if damage is high
 * - From S2, cells compete between repair (->S1) and death (->S3)
 */
struct ParamsS2 {
    double e2;      ///< E2/sigma: Repair state threshold (>0)
    double e3;      ///< E3/sigma: Death state threshold (>e2)
    double a;       ///< alpha/sigma: Damage per DSB (>0)
    double T21;     ///< Repair timescale in hours (>0) - OPTIMIZED
    double T23;     ///< Death timescale in hours (>0) - OPTIMIZED
    double k_error; ///< Misrepair rate per (DSB*hour) - stochastic misrepair channel

    ParamsS2() : e2(10.0), e3(25.0), a(1.0), T21(12.0), T23(48.0), k_error(0.0) {}
    ParamsS2(double e2_val, double e3_val, double a_val,
             double T21_val = 12.0, double T23_val = 48.0, double k_error_val = 0.0)
        : e2(e2_val), e3(e3_val), a(a_val), T21(T21_val), T23(T23_val), k_error(k_error_val) {}

    /**
     * @brief Check if parameters satisfy physical constraints
     */
    bool isValid() const {
        return (e2 > 0.0) && (e3 > e2) && (a > 0.0) && (T21 > 0.0) && (T23 > 0.0) && (k_error >= 0.0);
    }

    /**
     * @brief Print parameters
     */
    void print() const {
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "=== Repair Mediated Model Parameters (Reduced) ===\n";
        std::cout << "  e2      = E2/sigma = " << e2 << " (repair threshold)\n";
        std::cout << "  e3      = E3/sigma = " << e3 << " (death threshold)\n";
        std::cout << "  a       = alpha/sigma = " << a << " (damage per DSB)\n";
        std::cout << "  T21     = " << T21 << " hours (repair timescale) - OPTIMIZED\n";
        std::cout << "  T23     = " << T23 << " hours (death timescale) - OPTIMIZED\n";
        std::cout << "  k_error = " << k_error << " per (DSB*hour) (misrepair rate)\n";
        std::cout << "===============================================\n";
    }
};

/**
 * @brief Full physical parameters for Repair Mediated cell state model
 * 
 * These are the actual parameters needed by CellStateModel with S2 state:
 * - E1: Healthy state energy (reference, typically 0)
 * - E2: Repair/arrested state threshold
 * - E3: Death state threshold
 * - sigma: Gaussian width for state transitions
 * - alpha: Damage per DSB
 * - T21: Repair timescale (hours) - FROM OPTIMIZATION
 * - T23: Death-from-arrest timescale (hours) - FROM OPTIMIZATION
 */
struct CellStateModelParamsS2 {
    double E1;      ///< Healthy state energy (typically 0)
    double E2;      ///< Repair state threshold
    double E3;      ///< Death state threshold
    double sigma;   ///< Gaussian width for transitions
    double alpha;   ///< Damage per DSB
    double kappa;   ///< DSB per Gy
    double T21;     ///< Repair timescale (hours) - from optimization
    double T23;     ///< Death timescale (hours) - from optimization
    double k_error; ///< Misrepair rate per (DSB*hour) - stochastic misrepair channel

    CellStateModelParamsS2()
        : E1(0.0), E2(0.0), E3(0.0), sigma(1.0), alpha(0.0)
        , kappa(40.0), T21(12.0), T23(48.0), k_error(0.0) {}

    /**
     * @brief Create from reduced parameters by specifying sigma
     * @param reduced Fitted parameters (e2, e3, a, T21, T23, k_error) - all from optimization
     * @param sigma_value Chosen sigma value
     * @param kappa_value DSB per Gy
     * @return Full physical parameters
     */
    static CellStateModelParamsS2 fromReduced(const ParamsS2& reduced,
                                               double sigma_value,
                                               double kappa_value = 40.0) {
        CellStateModelParamsS2 params;
        params.E1 = 0.0;
        params.E2 = reduced.e2 * sigma_value;
        params.E3 = reduced.e3 * sigma_value;
        params.sigma = sigma_value;
        params.alpha = reduced.a * sigma_value;
        params.kappa = kappa_value;
        params.T21 = reduced.T21;      // From optimization result
        params.T23 = reduced.T23;      // From optimization result
        params.k_error = reduced.k_error;  // From optimization result (dimension-independent)
        return params;
    }

    /**
     * @brief Get critical DSB count for immediate death
     */
    double criticalDSBCount() const {
        return (alpha > 0) ? E3 / alpha : 0.0;
    }

    /**
     * @brief Get DSB count for entering repair state
     */
    double repairDSBCount() const {
        return (alpha > 0) ? E2 / alpha : 0.0;
    }

    /**
     * @brief Print parameters
     */
    void print() const {
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "=== Repair Mediated Cell State Model Parameters ===\n";
        std::cout << "  E1 (healthy state):     " << E1 << "\n";
        std::cout << "  E2 (repair threshold):  " << E2 << "\n";
        std::cout << "  E3 (death threshold):   " << E3 << "\n";
        std::cout << "  sigma (transition width): " << sigma << "\n";
        std::cout << "  alpha (damage per DSB): " << alpha << "\n";
        std::cout << "  kappa (DSB/Gy):         " << kappa << "\n";
        std::cout << "  T21 (repair time):      " << T21 << " hours\n";
        std::cout << "  T23 (death time):       " << T23 << " hours\n";
        std::cout << "  k_error (misrepair):    " << k_error << " per (DSB*hour)\n";
        std::cout << "\nDerived quantities:\n";
        if (alpha > 1e-9) {
            double repair_dsb = E2 / alpha;
            double death_dsb = E3 / alpha;
            std::cout << "  Repair DSB threshold (E2/alpha): " << std::setprecision(2) << repair_dsb << " DSBs\n";
            std::cout << "  Death DSB threshold (E3/alpha):  " << death_dsb << " DSBs\n";
            if (kappa > 1e-9) {
                std::cout << "  Repair dose (E2/alpha/kappa):    " << std::setprecision(3) << repair_dsb / kappa << " Gy\n";
                std::cout << "  Death dose (E3/alpha/kappa):     " << death_dsb / kappa << " Gy\n";
            }
        }
        std::cout << "================================================\n";
    }
};

} // namespace CellStateCalibration

#endif // DATA_TYPES_HH
