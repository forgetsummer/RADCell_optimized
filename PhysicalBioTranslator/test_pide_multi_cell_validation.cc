/**
 * @file test_pide_multi_cell_validation.cc
 * @brief Validation test using PIDE database 60Co photon survival data
 * 
 * This test:
 * 1. Loads 60Co photon survival data from PIDE database
 * 2. For each cell line with valid LQ parameters:
 *    - Calibrates the Repair Mediated model (S1→S2→S3) to LQ data
 *    - Runs the full cell state simulation with calibrated parameters
 *    - Compares simulated SF with LQ-generated data to validate calibration
 * 3. Outputs summary of validation results for all cell types
 * 
 * @author RADCellSimulation Project
 * @date 2026
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <omp.h>
#include <random>

#include "Cell.hh"
#include "CellStateModel.hh"
#include "CellLayoutInitializer.hh"

// Include calibration headers (with misrepair)
#include "RepairMediatedCalibrator_with_k_error.hh"
#include "DataTypes.hh"

// Include PIDE data reader
#include "PIDEDataReader.hh"

using namespace std;
using namespace CellStateCalibration;
using namespace PIDE;

//============================================================================
// LQ Model for generating survival data
//============================================================================
double LQ_survival(double dose, double alpha_LQ, double beta_LQ) {
    return exp(-alpha_LQ * dose - beta_LQ * dose * dose);
}

//============================================================================
// Structure to store validation results for each cell type
//============================================================================
struct CellValidationResult {
    string cellLine;
    string photonRadiation;
    double alpha_LQ;
    double beta_LQ;
    double alpha_beta_ratio;
    
    // Calibrated parameters
    double e2, e3, a, k_error;
    double T21, T23;
    bool calibrationConverged;
    
    // Validation results
    vector<double> doses;
    vector<double> sf_lq;
    vector<double> sf_calibrated;  // Analytical calibrated model prediction
    vector<double> sf_sim;
    vector<double> sf_sim_sd;
    double avg_log10_error;
    double max_log10_error;
};

//============================================================================
// Helper function to create directory
//============================================================================
void createDirectoryIfNotExists(const char* path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0755);
        cout << "Created directory: " << path << endl;
    }
}

//============================================================================
// Helper function to sanitize cell line name for use in filenames
// Replaces "/" with "_" and other problematic characters
//============================================================================
string sanitizeFilename(const string& name) {
    string result = name;
    for (size_t i = 0; i < result.length(); i++) {
        if (result[i] == '/' || result[i] == '\\' || result[i] == ':' || 
            result[i] == '*' || result[i] == '?' || result[i] == '"' ||
            result[i] == '<' || result[i] == '>' || result[i] == '|') {
            result[i] = '_';
        }
    }
    return result;
}

//============================================================================
// Helper structure for running statistics
//============================================================================
struct RunningStats {
    int n = 0;
    double sum = 0.0;
    double sumsq = 0.0;
    void add(double x) { n++; sum += x; sumsq += x * x; }
    double mean() const { return (n > 0) ? (sum / n) : 0.0; }
    double var() const {
        if (n < 2) return 0.0;
        const double m = mean();
        return std::max(0.0, (sumsq / n) - m * m);
    }
    double sd() const { return std::sqrt(var()); }
};

//============================================================================
// Run simulation for a single cell type with given LQ parameters
//============================================================================
CellValidationResult runValidationForCellType(
    const string& cellLine,
    const string& photonRadiation,
    double alpha_LQ,
    double beta_LQ,
    int cellNum,
    int numReplicates,
    bool checkpointEnabled,
    bool softSaturationEnabled,
    double repairRateScale,
    const string& outputDir,
    unsigned long long baseSeed,
    bool poissonDSBEnabled,
    bool useFittedTimescales,
    int checkpointMode,
    double lambdaCommit)
{
    CellValidationResult result;
    result.cellLine = cellLine;
    result.photonRadiation = photonRadiation;
    result.alpha_LQ = alpha_LQ;
    result.beta_LQ = beta_LQ;
    result.alpha_beta_ratio = (beta_LQ > 0) ? (alpha_LQ / beta_LQ) : 0.0;
    
    cout << "\n========================================" << endl;
    cout << "  Validating: " << cellLine << endl;
    cout << "  Radiation: " << photonRadiation << endl;
    cout << "  α = " << alpha_LQ << " Gy⁻¹" << endl;
    cout << "  β = " << beta_LQ << " Gy⁻²" << endl;
    cout << "  α/β = " << result.alpha_beta_ratio << " Gy" << endl;
    cout << "========================================" << endl;
    
    //------------------------------------------------------------------------
    // Generate fake data for calibration
    //------------------------------------------------------------------------
    vector<double> fakeDoses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
    vector<DataPoint> fakeData;
    
    for (double dose : fakeDoses) {
        double sf = LQ_survival(dose, alpha_LQ, beta_LQ);
        fakeData.push_back(DataPoint(dose, sf, 0.0));
    }
    
    //------------------------------------------------------------------------
    // Calibrate Repair Mediated Model WITH MISREPAIR
    //------------------------------------------------------------------------
    const double T_cellCycle = 10.0;  // hours
    
    FitConfigS2 cfg;
    cfg.kappa = 40.0;           // DSB/Gy
    cfg.dose_max_fit = 6.0;     // Fit all data
    cfg.T_assay_h = 27.8;       // Match simulation length
    cfg.setTimescalesToCellCycle(T_cellCycle);
    
    RepairMediatedMisrepairCalibrator calibrator(cfg);
    calibrator.setData(fakeData);
    calibrator.setDefaultUncertainties(0.1);
    
    // Initial guess
    ParamsS2 x0;
    x0.e2 = 2.0;
    x0.e3 = 10.0;
    x0.a = 0.05;
    x0.T21 = T_cellCycle;
    x0.T23 = T_cellCycle;
    x0.k_error = 1e-4;
    
    ParamsS2 step;
    step.e2 = 0.5;
    step.e3 = 2.0;
    step.a = 0.01;
    step.T21 = 1.0;
    step.T23 = 1.0;
    step.k_error = 5e-5;
    
    RepairMediatedMisrepairFitResult fitResult = calibrator.fit(x0, step);
    
    result.e2 = fitResult.params.e2;
    result.e3 = fitResult.params.e3;
    result.a = fitResult.params.a;
    result.k_error = fitResult.params.k_error;
    result.T21 = fitResult.params.T21;
    result.T23 = fitResult.params.T23;
    result.calibrationConverged = fitResult.converged;
    
    cout << "  Calibration " << (fitResult.converged ? "converged" : "FAILED") << endl;
    cout << "  e2=" << result.e2 << ", e3=" << result.e3 << ", a=" << result.a << endl;
    
    // Convert to physical parameters
    double sigma_chosen = 10.0;
    CellStateModelParamsS2 calibratedParams = CellStateModelParamsS2::fromReduced(
        fitResult.params, sigma_chosen, cfg.kappa);
    calibratedParams.k_error = fitResult.params.k_error;
    
    //------------------------------------------------------------------------
    // Run Cell State Simulation
    //------------------------------------------------------------------------
    double xDim = 5.0, yDim = 5.0, zDim = 0.0;
    double d = 0.05;
    
    double deltaT_cellPhaseUpdate = 60;
    double deltaT_cellStateUpdate = 60;
    double T = 10;
    int totalTimeStepNum = 10000;  // ~27.8 hours
    
    // Cell cycle parameters
    double tG1 = 1.5, sigmaG1 = 0.25;
    double tS = 6.0, sigmaS = 0.25;
    double tG2 = 1.5, sigmaG2 = 0.25;
    double tM = 1.0, sigmaM = 0.25;
    
    double alpha_model = calibratedParams.alpha;
    double E1 = calibratedParams.E1;
    double E2 = calibratedParams.E2;
    double E3 = calibratedParams.E3;
    double sigma_model = calibratedParams.sigma;
    
    double beta_model = 25.71;
    double fG1 = 1.0, fG2 = 1.015306122, fM = 1.015306122, fS = 0.816326;
    double f1 = 0.62, f2 = 0.38;
    double lambda1 = 3.31, lambda2 = 0.14;
    double lambda1_scaled = lambda1 * repairRateScale;
    double lambda2_scaled = lambda2 * repairRateScale;
    
    double tau1 = 1.0 / lambda1_scaled;
    double tau2 = 1.0 / lambda2_scaled;
    double effective_repair_time = f1 * tau1 + f2 * tau2;
    
    // If using fitted timescales, use calibrated T21/T23 instead of effective_repair_time
    double T21_to_use = useFittedTimescales ? fitResult.params.T21 : effective_repair_time;
    double T23_to_use = useFittedTimescales ? fitResult.params.T23 : effective_repair_time;
    
    double DSB_per_Gy = cfg.kappa;
    
    vector<double> simDoses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
    
    // Create Cell Layout
    CellLayoutInitializer layout;
    layout.SetRandomSeed(baseSeed);
    layout.SetCellHomeParamter(d, d, d);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);
    
    // Results storage
    result.doses = simDoses;
    result.sf_lq.resize(simDoses.size());
    result.sf_calibrated.resize(simDoses.size());
    result.sf_sim.resize(simDoses.size());
    result.sf_sim_sd.resize(simDoses.size());
    
    // Get calibrated model (analytical) predictions
    vector<double> sf_calibrated_vec = calibrator.predictSurvival(fitResult.params, simDoses);
    for (size_t i = 0; i < simDoses.size(); i++) {
        result.sf_calibrated[i] = sf_calibrated_vec[i];
    }
    
    const double sf_eps = std::max(1e-12, 0.5 / (static_cast<double>(cellNum) * static_cast<double>(numReplicates)));
    
    // Main parallel dose loop
    #pragma omp parallel for schedule(dynamic)
    for (size_t n = 0; n < simDoses.size(); n++)
    {
        double dose = simDoses[n];
        
        // Precompute initial DSB
        std::vector<int> initialDSB(cellNum, 0);
        double meanDSB = dose * DSB_per_Gy;
        
        if (poissonDSBEnabled && meanDSB > 0) {
            // Poisson DSB sampling: each cell gets a different DSB count
            // from Poisson distribution with mean = dose * DSB_per_Gy
            std::mt19937_64 dsbRng(baseSeed + n * 999999ULL);
            std::poisson_distribution<int> poissonDist(meanDSB);
            for (int cid = 0; cid < cellNum; cid++) {
                initialDSB[cid] = poissonDist(dsbRng);
            }
        } else {
            // Uniform DSB: all cells get the same DSB count
            for (int cid = 0; cid < cellNum; cid++) {
                initialDSB[cid] = static_cast<int>(meanDSB);
            }
        }
        
        RunningStats sfStats;
        
        for (int rep = 0; rep < numReplicates; rep++)
        {
            double clock_cellPhaseUpdate = 0;
            double clock_cellStateUpdate = 0;
            int cycle_cellPhaseUpdate = 0;
            int cycle_cellStateUpdate = 0;
            
            Cell testCell;
            testCell.CellConstruct("Epithelia", "Cytoplasm Nucleus", "Sphere", "Blue Green");
            
            CellCycleParameter cyclePara;
            cyclePara.mTG1 = tG1; cyclePara.sigmaTG1 = sigmaG1;
            cyclePara.mTS = tS; cyclePara.sigmaTS = sigmaS;
            cyclePara.mTG2 = tG2; cyclePara.sigmaTG2 = sigmaG2;
            cyclePara.mTM = tM; cyclePara.sigmaTM = sigmaM;
            cyclePara.fG1 = fG1; cyclePara.fS = fS; cyclePara.fG2 = fG2; cyclePara.fM = fM;
            
            CellStateParameter statePara;
            statePara.E1 = E1;
            statePara.E2 = E2;
            statePara.E3 = E3;
            statePara.sigma = sigma_model;
            statePara.alpha = alpha_model;
            statePara.beta = beta_model;
            statePara.To12 = T_cellCycle;
            statePara.To13 = T_cellCycle;
            statePara.To21 = T21_to_use;
            statePara.To23 = T23_to_use;
            
            CellDNADamageRepairParameter dnaRepairPara;
            dnaRepairPara.f1 = f1; dnaRepairPara.f2 = f2;
            dnaRepairPara.lambda1 = lambda1_scaled; dnaRepairPara.lambda2 = lambda2_scaled;
            
            CellBystanderSignalParameter bystanderPara;
            bystanderPara.Rt = 5000; bystanderPara.mu = 200;
            
            testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);
            
            CellStateModel myCellState;
            myCellState.SetRandomSeed(baseSeed + static_cast<unsigned long long>(n) * 100000ULL +
                                      static_cast<unsigned long long>(rep) * 1000ULL);
            myCellState.TissueGeometryInitialization(xDim, yDim, zDim, d);
            myCellState.CellStateModelParameterSetup(testCell);
            myCellState.SetMisrepairRate("Epithelia", calibratedParams.k_error);
            myCellState.SetUpContactInhibition(true);
            myCellState.SetCheckpointEnabled(checkpointEnabled);
            myCellState.SetSoftSaturationEnabled(softSaturationEnabled);
            myCellState.SetCheckpointMode(checkpointMode);
            myCellState.SetCommitmentHazardRate(lambdaCommit);
            
            // Initialize cells
            for (int i = 0; i < cellNum; i++) {
                double cX = layout.GetCellPositionX(i) + xDim / 2.0;
                double cY = layout.GetCellPositionY(i) + yDim / 2.0;
                double cZ = layout.GetCellPositionZ(i) + zDim / 2.0;
                myCellState.CellPositionInitialization(i, cX, cY, cZ);
                myCellState.CellTypeInitialiation(i, testCell);
                myCellState.CellPhaseInitializationRandom(i);
                myCellState.CellStateInitialization(i, "S1");
            }
            
            // Run simulation
            for (int i = 0; i < totalTimeStepNum; i++)
            {
                clock_cellStateUpdate += T;
                clock_cellPhaseUpdate += T;
                
                map<int, string> phaseMap = myCellState.GetCellPhase();
                
                if (clock_cellPhaseUpdate >= deltaT_cellPhaseUpdate) {
                    cycle_cellPhaseUpdate++;
                    for (auto& cell : phaseMap) {
                        myCellState.CellPhaseUpdate(cell.first, true, deltaT_cellPhaseUpdate, 1);
                    }
                    clock_cellPhaseUpdate = 0;
                }
                
                if (clock_cellStateUpdate >= deltaT_cellStateUpdate) {
                    cycle_cellStateUpdate++;
                    map<int, string> stateMap = myCellState.GetCellState();
                    
                    for (auto& cell : stateMap) {
                        int currentDSB = 0;
                        if (cycle_cellStateUpdate == 1) {
                            const int cid = cell.first;
                            if (cid >= 0 && cid < cellNum) currentDSB = initialDSB[cid];
                        }
                        myCellState.CellStateUpdate(cell.first, currentDSB, 0, deltaT_cellStateUpdate, 1);
                    }
                    clock_cellStateUpdate = 0;
                }
            }
            
            // Colony formation assay
            map<int, string> finalStateMap = myCellState.GetCellState();
            map<int, int> cellAncestryIDMap = myCellState.GetCellAncestryID();
            map<int, int> cellColonySizeMap;
            for (int i = 0; i < cellNum; i++) cellColonySizeMap[i] = 0;
            
            for (const auto& cell : finalStateMap)
            {
                const int cellID = cell.first;
                const string& st = cell.second;
                if (st != "S1") continue;
                
                auto itA = cellAncestryIDMap.find(cellID);
                if (itA == cellAncestryIDMap.end()) continue;
                
                const int ancestryID = itA->second;
                cellColonySizeMap[ancestryID]++;
            }
            
            int colonySizeThreshold = 4;
            int colonyNumber = 0;
            for (const auto& colony : cellColonySizeMap) {
                if (colony.second >= colonySizeThreshold) colonyNumber++;
            }
            
            const double survivalFraction = static_cast<double>(colonyNumber) / static_cast<double>(cellNum);
            sfStats.add(survivalFraction);
        }
        
        // Store results
        result.sf_lq[n] = LQ_survival(dose, alpha_LQ, beta_LQ);
        result.sf_sim[n] = sfStats.mean();
        result.sf_sim_sd[n] = sfStats.sd();
    }
    
    // Calculate errors
    double total_log10_error = 0.0;
    result.max_log10_error = 0.0;
    
    for (size_t n = 0; n < simDoses.size(); n++) {
        double log10_lq = std::log10(result.sf_lq[n] + sf_eps);
        double log10_sim = std::log10(result.sf_sim[n] + sf_eps);
        double log10_error = std::fabs(log10_sim - log10_lq);
        
        total_log10_error += log10_error;
        if (log10_error > result.max_log10_error) {
            result.max_log10_error = log10_error;
        }
    }
    
    result.avg_log10_error = total_log10_error / simDoses.size();
    
    cout << "  Avg |Δlog10(SF)|: " << fixed << setprecision(3) << result.avg_log10_error << endl;
    cout << "  Max |Δlog10(SF)|: " << result.max_log10_error << endl;
    
    // Save cell-specific results (sanitize filename to handle "/" in cell names like C3H10T1/2)
    string safeCellLine = sanitizeFilename(cellLine);
    string cellOutputFile = outputDir + "/pide_validation_" + safeCellLine + ".csv";
    ofstream cellFile(cellOutputFile);
    cellFile << "Dose,SF_LQ,SF_Calibrated,SF_Simulated,SF_Simulated_SD,AbsLog10Error" << endl;
    
    for (size_t n = 0; n < simDoses.size(); n++) {
        double log10_error = std::fabs(std::log10(result.sf_sim[n] + sf_eps) - std::log10(result.sf_lq[n] + sf_eps));
        cellFile << simDoses[n] << "," << result.sf_lq[n] << "," << result.sf_calibrated[n] << "," 
                 << result.sf_sim[n] << "," << result.sf_sim_sd[n] << "," << log10_error << endl;
    }
    cellFile.close();
    
    return result;
}

//============================================================================
// Main function
//============================================================================
int main(int argc, char** argv)
{
    // Parse command-line arguments
    bool checkpointEnabled = true;
    bool softSaturationEnabled = true;
    double repairRateScale = 0.1;
    int numReplicates = 3;  // Fewer replicates for faster multi-cell testing
    int cellNum = 500;      // Fewer cells for faster testing
    int maxCellTypes = 10;  // Maximum number of cell types to test
    string pideDir = "../PIDE3.4";
    string outputDir = "./pide_multi_cell_validation_output";
    bool poissonDSBEnabled = false;      // Ablation: Poisson DSB sampling per cell
    bool useFittedTimescales = false;    // Ablation: Use calibrated T21/T23
    int checkpointMode = 0;              // 0 = gating-catastrophe, 1 = attempt-based
    double lambdaCommit = 0.0;           // Background commitment hazard (mode 1 only)
    for (int i = 1; i < argc; i++) {
        string arg = string(argv[i]);
        
        if (arg == "--pide-dir" && i + 1 < argc) {
            pideDir = argv[++i];
        } else if (arg == "--output-dir" && i + 1 < argc) {
            outputDir = argv[++i];
        } else if (arg == "--replicates" && i + 1 < argc) {
            numReplicates = atoi(argv[++i]);
        } else if (arg == "--cells" && i + 1 < argc) {
            cellNum = atoi(argv[++i]);
        } else if (arg == "--max-cell-types" && i + 1 < argc) {
            maxCellTypes = atoi(argv[++i]);
        } else if (arg == "--no-checkpoint") {
            checkpointEnabled = false;
        } else if (arg == "--repair-scale" && i + 1 < argc) {
            repairRateScale = atof(argv[++i]);
        } else if (arg == "--poisson-dsb") {
            poissonDSBEnabled = true;
        } else if (arg == "--use-fitted-timescales") {
            useFittedTimescales = true;
        } else if (arg == "--checkpoint-mode" && i + 1 < argc) {
            checkpointMode = atoi(argv[++i]);
        } else if (arg == "--lambda-commit" && i + 1 < argc) {
            lambdaCommit = atof(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            cout << "Usage: " << argv[0] << " [OPTIONS]" << endl;
            cout << "Options:" << endl;
            cout << "  --pide-dir PATH: Path to PIDE3.4 directory (default: ../PIDE3.4)" << endl;
            cout << "  --output-dir PATH: Output directory (default: ./pide_multi_cell_validation_output)" << endl;
            cout << "  --replicates N: Replicates per dose (default: 3)" << endl;
            cout << "  --cells N: Number of cells (default: 500)" << endl;
            cout << "  --max-cell-types N: Maximum cell types to test (default: 10)" << endl;
            cout << "  --no-checkpoint: Disable checkpoint" << endl;
            cout << "  --repair-scale VALUE: Repair rate scale factor (default: 0.1)" << endl;
            cout << "  --poisson-dsb: Enable Poisson DSB sampling per cell" << endl;
            cout << "  --use-fitted-timescales: Use calibrated T21/T23 instead of fixed values" << endl;
            cout << "  --checkpoint-mode N: Checkpoint mode (0=gating-catastrophe, 1=attempt-based)" << endl;
            cout << "  --lambda-commit VALUE: Background commitment hazard rate (mode 1 only, default: 0.0)" << endl;
            return 0;
        }
    }
    
    cout << "============================================" << endl;
    cout << "  PIDE Multi-Cell Validation Test" << endl;
    cout << "  60Co Photon Survival Data" << endl;
    cout << "============================================" << endl;
    cout << "Configuration:" << endl;
    cout << "  PIDE directory: " << pideDir << endl;
    cout << "  Output directory: " << outputDir << endl;
    cout << "  Replicates per dose: " << numReplicates << endl;
    cout << "  Cells per simulation: " << cellNum << endl;
    cout << "  Max cell types: " << maxCellTypes << endl;
    cout << "  Checkpoint: " << (checkpointEnabled ? "enabled" : "disabled") << endl;
    cout << "  Repair rate scale: " << repairRateScale << endl;
    cout << "  Poisson DSB: " << (poissonDSBEnabled ? "enabled" : "disabled") << endl;
    cout << "  Fitted timescales: " << (useFittedTimescales ? "enabled" : "disabled") << endl;
    cout << "  Checkpoint mode: " << checkpointMode << " (" 
         << (checkpointMode == 0 ? "gating-catastrophe" : "attempt-based") << ")" << endl;
    cout << "  Lambda commit: " << lambdaCommit << endl;
    
    // Create output directory
    createDirectoryIfNotExists(outputDir.c_str());
    
    //========================================================================
    // Load PIDE data
    //========================================================================
    cout << "\n--- Loading PIDE Database ---" << endl;
    
    PIDEDataReader reader;
    if (!reader.loadData(pideDir)) {
        cerr << "ERROR: Failed to load PIDE data from: " << pideDir << endl;
        cerr << "Make sure to run convert_pide_xlsx.py first to generate the CSV file." << endl;
        return 1;
    }
    
    reader.printSummary();
    
    //========================================================================
    // Find 60Co experiments with valid LQ parameters
    //========================================================================
    cout << "\n--- Finding 60Co Experiments ---" << endl;
    
    // Get all experiments and filter for 60Co photon data
    auto allExps = reader.getExperimentsByCellLine("");  // Get all
    
    // We need to iterate through all experiments since getPhotonExperiments may not exist
    // Instead, filter manually for photon radiation containing "60Co" or "Co60" or "Co-60"
    vector<ExperimentEntry> co60Exps;
    map<string, ExperimentEntry> uniqueCellLines;  // cellLine -> first valid experiment
    
    for (size_t i = 0; i < reader.getNumExperiments(); i++) {
        const ExperimentEntry* exp = reader.getExperimentByID(i + 1);  // ExpIDs are 1-based
        if (!exp) continue;
        
        // Check if photon radiation is 60Co
        string photonRad = exp->photonRadiation;
        bool is60Co = (photonRad.find("60Co") != string::npos ||
                       photonRad.find("Co60") != string::npos ||
                       photonRad.find("Co-60") != string::npos);
        
        if (!is60Co) continue;
        
        // Check for valid photon LQ parameters
        bool hasValidLQ = exp->photonLQ_paper.isValid || exp->photonLQ_fit.isValid;
        if (!hasValidLQ) continue;
        
        // Use unique cell lines only (first occurrence)
        if (uniqueCellLines.find(exp->cellLine) == uniqueCellLines.end()) {
            uniqueCellLines[exp->cellLine] = *exp;
            co60Exps.push_back(*exp);
        }
    }
    
    cout << "Found " << co60Exps.size() << " unique cell lines with 60Co photon data" << endl;
    
    if (co60Exps.empty()) {
        cerr << "ERROR: No 60Co experiments found with valid LQ parameters!" << endl;
        return 1;
    }
    
    // Limit to maxCellTypes
    if (co60Exps.size() > static_cast<size_t>(maxCellTypes)) {
        co60Exps.resize(maxCellTypes);
        cout << "Limited to " << maxCellTypes << " cell types for testing" << endl;
    }
    
    // Print found experiments
    cout << "\nCell lines to validate:" << endl;
    for (const auto& exp : co60Exps) {
        const LQParameters& lq = exp.photonLQ_paper.isValid ? exp.photonLQ_paper : exp.photonLQ_fit;
        cout << "  - " << exp.cellLine << ": α=" << lq.alpha << ", β=" << lq.beta << endl;
    }
    
    //========================================================================
    // Run validation for each cell type
    //========================================================================
    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    cout << "\n--- Running Validation (OpenMP threads: " << numThreads << ") ---" << endl;
    
    double wall_start = omp_get_wtime();
    
    vector<CellValidationResult> allResults;
    unsigned long long baseSeed = 12345ULL;
    
    for (size_t i = 0; i < co60Exps.size(); i++) {
        const auto& exp = co60Exps[i];
        const LQParameters& lq = exp.photonLQ_paper.isValid ? exp.photonLQ_paper : exp.photonLQ_fit;
        
        CellValidationResult result = runValidationForCellType(
            exp.cellLine,
            exp.photonRadiation,
            lq.alpha,
            lq.beta,
            cellNum,
            numReplicates,
            checkpointEnabled,
            softSaturationEnabled,
            repairRateScale,
            outputDir,
            baseSeed + i * 1000000ULL,
            poissonDSBEnabled,
            useFittedTimescales,
            checkpointMode,
            lambdaCommit
        );
        
        allResults.push_back(result);
    }
    
    double wall_end = omp_get_wtime();
    
    //========================================================================
    // Save summary results
    //========================================================================
    cout << "\n--- Saving Summary Results ---" << endl;
    
    // Summary CSV
    string summaryFile = outputDir + "/pide_validation_summary.csv";
    ofstream summaryCSV(summaryFile);
    summaryCSV << "CellLine,PhotonRadiation,Alpha_LQ,Beta_LQ,AlphaBetaRatio,"
               << "e2,e3,a,k_error,T21,T23,CalibConverged,"
               << "AvgLog10Error,MaxLog10Error" << endl;
    
    for (const auto& r : allResults) {
        summaryCSV << r.cellLine << "," << r.photonRadiation << ","
                   << r.alpha_LQ << "," << r.beta_LQ << "," << r.alpha_beta_ratio << ","
                   << r.e2 << "," << r.e3 << "," << r.a << "," << r.k_error << ","
                   << r.T21 << "," << r.T23 << "," << (r.calibrationConverged ? "yes" : "no") << ","
                   << r.avg_log10_error << "," << r.max_log10_error << endl;
    }
    summaryCSV.close();
    
    // Parameters text file
    string paramsFile = outputDir + "/pide_validation_params.txt";
    ofstream paramsTxt(paramsFile);
    paramsTxt << "=== PIDE Multi-Cell Validation Parameters ===" << endl;
    paramsTxt << "\nConfiguration:" << endl;
    paramsTxt << "  PIDE directory: " << pideDir << endl;
    paramsTxt << "  Replicates per dose: " << numReplicates << endl;
    paramsTxt << "  Cells per simulation: " << cellNum << endl;
    paramsTxt << "  Checkpoint enabled: " << (checkpointEnabled ? "yes" : "no") << endl;
    paramsTxt << "  Checkpoint mode: " << checkpointMode << " (" 
              << (checkpointMode == 0 ? "gating-catastrophe" : "attempt-based") << ")" << endl;
    paramsTxt << "  Lambda commit: " << lambdaCommit << endl;
    paramsTxt << "  Repair rate scale: " << repairRateScale << endl;
    paramsTxt << "\nTotal wall time: " << fixed << setprecision(2) << (wall_end - wall_start) << " seconds" << endl;
    
    paramsTxt << "\n=== Results by Cell Line ===" << endl;
    for (const auto& r : allResults) {
        paramsTxt << "\n--- " << r.cellLine << " ---" << endl;
        paramsTxt << "  Radiation: " << r.photonRadiation << endl;
        paramsTxt << "  LQ: α=" << r.alpha_LQ << " Gy⁻¹, β=" << r.beta_LQ << " Gy⁻²" << endl;
        paramsTxt << "  α/β = " << r.alpha_beta_ratio << " Gy" << endl;
        paramsTxt << "  Calibrated: e2=" << r.e2 << ", e3=" << r.e3 << ", a=" << r.a << endl;
        paramsTxt << "  k_error=" << r.k_error << ", T21=" << r.T21 << ", T23=" << r.T23 << endl;
        paramsTxt << "  Converged: " << (r.calibrationConverged ? "yes" : "no") << endl;
        paramsTxt << "  Avg |Δlog10(SF)|: " << r.avg_log10_error << endl;
        paramsTxt << "  Max |Δlog10(SF)|: " << r.max_log10_error << endl;
    }
    paramsTxt.close();
    
    //========================================================================
    // Print final summary
    //========================================================================
    cout << "\n============================================" << endl;
    cout << "  VALIDATION COMPLETE" << endl;
    cout << "============================================" << endl;
    cout << "Cell types validated: " << allResults.size() << endl;
    cout << "Total wall time: " << fixed << setprecision(2) << (wall_end - wall_start) << " seconds" << endl;
    
    cout << "\nSummary of errors:" << endl;
    cout << setw(15) << "Cell Line" << setw(12) << "α" << setw(12) << "β" 
         << setw(15) << "Avg|Δlog10|" << setw(15) << "Max|Δlog10|" << endl;
    cout << string(69, '-') << endl;
    
    double totalAvgError = 0.0;
    for (const auto& r : allResults) {
        cout << setw(15) << r.cellLine 
             << setw(12) << fixed << setprecision(4) << r.alpha_LQ 
             << setw(12) << r.beta_LQ
             << setw(15) << setprecision(3) << r.avg_log10_error
             << setw(15) << r.max_log10_error << endl;
        totalAvgError += r.avg_log10_error;
    }
    
    cout << string(69, '-') << endl;
    double overallAvg = totalAvgError / allResults.size();
    cout << "Overall average |Δlog10(SF)|: " << fixed << setprecision(3) << overallAvg 
         << " (≈×" << setprecision(2) << std::pow(10.0, overallAvg) << ")" << endl;
    
    cout << "\nOutput files:" << endl;
    cout << "  - " << summaryFile << endl;
    cout << "  - " << paramsFile << endl;
    cout << "  - " << outputDir << "/pide_validation_*.csv (per cell line)" << endl;
    
    return 0;
}

