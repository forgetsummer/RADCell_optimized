/**
 * @file PIDEDataReader.hh
 * @brief Reader for Particle Irradiation Data Ensemble (PIDE) database
 * 
 * This module provides functionality to read and parse cell survival data from
 * the GSI PIDE database (version 3.4). The database contains in-vitro cell 
 * survival experiments with Linear-Quadratic (LQ) parameters.
 * 
 * Reference:
 * - Friedrich T, et al. J Radiat Res. 2013;54(3):494-514
 * - Friedrich T, et al. J Radiat Res. 2021;62(4):645-655
 * 
 * @author RADCellSimulation Project
 * @date 2026
 */

#ifndef PIDE_DATA_READER_HH
#define PIDE_DATA_READER_HH

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <fstream>
#include <sstream>
#include <cmath>

namespace PIDE {

/**
 * @struct LQParameters
 * @brief Linear-Quadratic model parameters
 * 
 * Survival fraction S(D) = exp(-(alpha*D + beta*D^2))
 */
struct LQParameters {
    double alpha = 0.0;  ///< Linear coefficient (Gy^-1)
    double beta = 0.0;   ///< Quadratic coefficient (Gy^-2)
    bool isValid = false;
    
    /// Calculate survival fraction at given dose
    double survivalFraction(double dose_Gy) const {
        if (!isValid) return -1.0;
        return std::exp(-(alpha * dose_Gy + beta * dose_Gy * dose_Gy));
    }
    
    /// Calculate D10 (dose for 10% survival)
    double D10() const {
        if (!isValid || (alpha == 0 && beta == 0)) return -1.0;
        // Solve: exp(-(alpha*D + beta*D^2)) = 0.1
        // alpha*D + beta*D^2 = ln(10) ≈ 2.303
        return solveForDose(0.1);
    }
    
    /// Calculate D37 (dose for 37% survival, mean lethal dose)
    double D37() const {
        if (!isValid || (alpha == 0 && beta == 0)) return -1.0;
        return solveForDose(0.37);
    }
    
    /// Calculate D50 (dose for 50% survival)
    double D50() const {
        if (!isValid || (alpha == 0 && beta == 0)) return -1.0;
        return solveForDose(0.5);
    }
    
private:
    double solveForDose(double survivalTarget) const {
        // Solve: alpha*D + beta*D^2 = -ln(survivalTarget)
        double target = -std::log(survivalTarget);
        if (beta == 0) return target / alpha;
        // Quadratic formula: D = (-alpha + sqrt(alpha^2 + 4*beta*target)) / (2*beta)
        double discriminant = alpha * alpha + 4.0 * beta * target;
        if (discriminant < 0) return -1.0;
        return (-alpha + std::sqrt(discriminant)) / (2.0 * beta);
    }
};

/**
 * @struct DoseSurvivalPoint
 * @brief A single dose-survival data point
 */
struct DoseSurvivalPoint {
    double dose_Gy = 0.0;
    double survivalFraction = 0.0;
};

/**
 * @struct RawSurvivalCurve
 * @brief Raw dose-response data from experiments
 */
struct RawSurvivalCurve {
    int expID = 0;
    int publicationID = 0;
    std::string publicationName;
    int experimentNum = 0;
    int photonExpNum = 0;  // Only for ion data - reference photon experiment
    std::vector<DoseSurvivalPoint> dataPoints;
    bool hasData = false;
    
    /// Get dose range
    std::pair<double, double> getDoseRange() const {
        if (dataPoints.empty()) return {0.0, 0.0};
        double minD = dataPoints[0].dose_Gy, maxD = dataPoints[0].dose_Gy;
        for (const auto& pt : dataPoints) {
            if (pt.dose_Gy < minD) minD = pt.dose_Gy;
            if (pt.dose_Gy > maxD) maxD = pt.dose_Gy;
        }
        return {minD, maxD};
    }
};

/**
 * @struct ExperimentEntry
 * @brief Complete experiment entry from PIDE database
 */
struct ExperimentEntry {
    // Identifiers
    int expID = 0;
    int publicationID = 0;
    std::string publicationName;
    int ionExpNum = 0;
    int photonExpNum = 0;
    
    // Cell information
    std::string cellLine;
    char cellClass = 'n';     ///< 't' = tumor, 'n' = normal
    char cellOrigin = 'h';    ///< 'h' = human, 'r' = rodent
    std::string cellCycle;    ///< Phase or 'a' for asynchronous
    double dnaContent = 6.0;  ///< Genomic length (10^9 bp)
    
    // Radiation information
    std::string photonRadiation;
    std::string ion;
    int charge = 0;
    char irradiationConditions = 'm';  ///< 'm' = monoenergetic, 's' = SOBP
    double LET_keV_um = 0.0;
    double energy_MeV_u = 0.0;
    
    // LQ parameters from publication
    LQParameters photonLQ_paper;
    LQParameters ionLQ_paper;
    
    // LQ parameters from PIDE fit
    LQParameters photonLQ_fit;
    LQParameters ionLQ_fit;
    
    /// Calculate RBE at given survival level using D10, D37, or D50
    double RBE_D10(bool useFitParams = true) const {
        const auto& photonLQ = useFitParams ? photonLQ_fit : photonLQ_paper;
        const auto& ionLQ = useFitParams ? ionLQ_fit : ionLQ_paper;
        
        if (!photonLQ.isValid || !ionLQ.isValid) return -1.0;
        
        double D_photon = photonLQ.D10();
        double D_ion = ionLQ.D10();
        
        if (D_ion <= 0) return -1.0;
        return D_photon / D_ion;
    }
    
    double RBE_D37(bool useFitParams = true) const {
        const auto& photonLQ = useFitParams ? photonLQ_fit : photonLQ_paper;
        const auto& ionLQ = useFitParams ? ionLQ_fit : ionLQ_paper;
        
        if (!photonLQ.isValid || !ionLQ.isValid) return -1.0;
        
        double D_photon = photonLQ.D37();
        double D_ion = ionLQ.D37();
        
        if (D_ion <= 0) return -1.0;
        return D_photon / D_ion;
    }
};

/**
 * @class PIDEDataReader
 * @brief Main class for reading and querying PIDE database
 */
class PIDEDataReader {
public:
    PIDEDataReader();
    ~PIDEDataReader() = default;
    
    /**
     * @brief Load PIDE data from directory
     * @param pideDirectory Path to PIDE3.4 folder
     * @return true if successful
     */
    bool loadData(const std::string& pideDirectory);
    
    /**
     * @brief Load main database from CSV file
     * @param csvFile Path to CSV file (converted from xlsx)
     * @return true if successful
     */
    bool loadMainDatabase(const std::string& csvFile);
    
    /**
     * @brief Load photon raw data
     * @param datFile Path to PhotonRawData.dat file
     * @return true if successful
     */
    bool loadPhotonRawData(const std::string& datFile);
    
    /**
     * @brief Load ion raw data
     * @param datFile Path to IonRawData.dat file
     * @return true if successful
     */
    bool loadIonRawData(const std::string& datFile);
    
    // Query methods
    
    /**
     * @brief Get all experiments for a specific cell line
     */
    std::vector<ExperimentEntry> getExperimentsByCellLine(const std::string& cellLine) const;
    
    /**
     * @brief Get all experiments for a specific ion
     */
    std::vector<ExperimentEntry> getExperimentsByIon(const std::string& ion) const;
    
    /**
     * @brief Get experiments by LET range
     */
    std::vector<ExperimentEntry> getExperimentsByLETRange(double minLET, double maxLET) const;
    
    /**
     * @brief Get experiments by cell type (tumor/normal)
     */
    std::vector<ExperimentEntry> getExperimentsByCellClass(char cellClass) const;
    
    /**
     * @brief Get experiments by cell origin (human/rodent)
     */
    std::vector<ExperimentEntry> getExperimentsByCellOrigin(char cellOrigin) const;
    
    /**
     * @brief Get photon-only experiments (for reference)
     */
    std::vector<ExperimentEntry> getPhotonExperiments() const;
    
    /**
     * @brief Get raw survival curve for photon experiment
     */
    RawSurvivalCurve getPhotonRawData(int publicationID, int photonExpNum) const;
    
    /**
     * @brief Get raw survival curve for ion experiment
     */
    RawSurvivalCurve getIonRawData(int expID) const;
    
    /**
     * @brief Get all unique cell lines in database
     */
    std::vector<std::string> getUniqueCellLines() const;
    
    /**
     * @brief Get all unique ions in database
     */
    std::vector<std::string> getUniqueIons() const;
    
    /**
     * @brief Get experiment by ID
     */
    const ExperimentEntry* getExperimentByID(int expID) const;
    
    /**
     * @brief Get total number of experiments
     */
    size_t getNumExperiments() const { return experiments_.size(); }
    
    /**
     * @brief Check if data is loaded
     */
    bool isLoaded() const { return dataLoaded_; }
    
    // Validation methods
    
    /**
     * @brief Compare simulation results with PIDE data
     * @param simDoses Vector of simulated doses
     * @param simSF Vector of simulated survival fractions
     * @param expID PIDE experiment ID for comparison
     * @return Mean squared error, or -1 if comparison not possible
     */
    double compareSurvivalCurve(const std::vector<double>& simDoses,
                                 const std::vector<double>& simSF,
                                 int expID) const;
    
    /**
     * @brief Fit LQ parameters to simulation data
     * @param doses Vector of doses
     * @param survivalFractions Vector of survival fractions
     * @return Fitted LQ parameters
     */
    static LQParameters fitLQParameters(const std::vector<double>& doses,
                                        const std::vector<double>& survivalFractions);
    
    /**
     * @brief Print summary of loaded data
     */
    void printSummary() const;
    
private:
    std::vector<ExperimentEntry> experiments_;
    std::map<std::pair<int, int>, RawSurvivalCurve> photonRawData_;  // (pubID, photonExpNum) -> data
    std::map<int, RawSurvivalCurve> ionRawData_;  // expID -> data
    bool dataLoaded_ = false;
    
    // Helper methods
    std::vector<std::string> splitLine(const std::string& line, char delimiter) const;
    double parseDouble(const std::string& str) const;
    int parseInt(const std::string& str) const;
    std::string trim(const std::string& str) const;
};

} // namespace PIDE

#endif // PIDE_DATA_READER_HH
