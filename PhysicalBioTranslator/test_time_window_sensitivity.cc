/**
 * @file test_time_window_sensitivity.cc
 * @brief Time window sensitivity test for validating repair-time mechanism
 * 
 * This test validates the hypothesis that high-dose SF divergence is caused by
 * S2→S1 repair recovery over time. It runs simulations at different assay time
 * windows (12h, 18h, 27.8h, 40h) with checkpoint on/off.
 * 
 * Prediction: If mechanism is correct, high-dose SF bias should increase with
 * assay time (more chance for S2->S1 recovery), and this effect should be
 * amplified when checkpoint is enabled.
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

using namespace std;
using namespace CellStateCalibration;

//============================================================================
// LQ Model for generating survival data
//============================================================================
double LQ_survival(double dose, double alpha_LQ, double beta_LQ) {
    return exp(-alpha_LQ * dose - beta_LQ * dose * dose);
}

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
// Structure to store results for a single assay time run
//============================================================================
struct TimeSensitivityResult {
    double assayTimeHours;
    bool checkpointEnabled;
    string cellLine;
    double alpha_LQ;
    double beta_LQ;
    
    // Results per dose
    vector<double> doses;
    vector<double> sf_lq;
    vector<double> sf_calibrated;
    vector<double> sf_sim;
    vector<double> sf_sim_sd;
    
    // Error metrics
    double avg_log10_error;
    double max_log10_error;
    double high_dose_avg_error;  // Error for doses >= 4 Gy
    double low_dose_avg_error;   // Error for doses < 4 Gy
};

//============================================================================
// Run simulation for given assay time
//============================================================================
TimeSensitivityResult runTimeSensitivityTest(
    double assayTimeHours,
    double alpha_LQ,
    double beta_LQ,
    const string& cellLine,
    int cellNum,
    int numReplicates,
    bool checkpointEnabled,
    double repairRateScale,
    unsigned long long baseSeed)
{
    TimeSensitivityResult result;
    result.assayTimeHours = assayTimeHours;
    result.checkpointEnabled = checkpointEnabled;
    result.cellLine = cellLine;
    result.alpha_LQ = alpha_LQ;
    result.beta_LQ = beta_LQ;
    
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
    // Calibrate Repair Mediated Model
    //------------------------------------------------------------------------
    const double T_cellCycle = 10.0;  // hours
    
    FitConfigS2 cfg;
    cfg.kappa = 40.0;           // DSB/Gy
    cfg.dose_max_fit = 6.0;     // Fit all data
    cfg.T_assay_h = assayTimeHours;  // Use the specified assay time!
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
    double T = 10;  // seconds per time step
    
    // Calculate time steps from assay time
    int totalTimeStepNum = static_cast<int>(assayTimeHours * 3600.0 / T);
    
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
        
        // Precompute initial DSB (uniform for all cells)
        std::vector<int> initialDSB(cellNum, 0);
        double meanDSB = dose * DSB_per_Gy;
        for (int cid = 0; cid < cellNum; cid++) {
            initialDSB[cid] = static_cast<int>(meanDSB);
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
            statePara.To21 = effective_repair_time;
            statePara.To23 = effective_repair_time;
            
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
            myCellState.SetSoftSaturationEnabled(true);
            
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
            
            // Run simulation for the specified assay time
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
    double high_dose_error = 0.0;
    double low_dose_error = 0.0;
    int high_dose_count = 0;
    int low_dose_count = 0;
    result.max_log10_error = 0.0;
    
    for (size_t n = 0; n < simDoses.size(); n++) {
        double log10_lq = std::log10(result.sf_lq[n] + sf_eps);
        double log10_sim = std::log10(result.sf_sim[n] + sf_eps);
        double log10_error = std::fabs(log10_sim - log10_lq);
        
        total_log10_error += log10_error;
        if (log10_error > result.max_log10_error) {
            result.max_log10_error = log10_error;
        }
        
        // Separate high-dose (>= 4 Gy) and low-dose (< 4 Gy) errors
        if (simDoses[n] >= 4.0) {
            high_dose_error += log10_error;
            high_dose_count++;
        } else {
            low_dose_error += log10_error;
            low_dose_count++;
        }
    }
    
    result.avg_log10_error = total_log10_error / simDoses.size();
    result.high_dose_avg_error = (high_dose_count > 0) ? high_dose_error / high_dose_count : 0.0;
    result.low_dose_avg_error = (low_dose_count > 0) ? low_dose_error / low_dose_count : 0.0;
    
    return result;
}

//============================================================================
// Main function
//============================================================================
int main(int argc, char** argv)
{
    // Default parameters (HS-23 cell line)
    double alpha_LQ = 0.34;   // HS-23 alpha
    double beta_LQ = 0.028;   // HS-23 beta
    string cellLine = "HS-23";
    
    // Parse command-line arguments
    bool checkpointEnabled = true;
    double repairRateScale = 0.1;
    int numReplicates = 5;
    int cellNum = 1000;
    double assayTimeHours = 27.8;  // Default
    string outputDir = "./time_sensitivity_output";
    
    for (int i = 1; i < argc; i++) {
        string arg = string(argv[i]);
        
        if (arg == "--assay-time" && i + 1 < argc) {
            assayTimeHours = atof(argv[++i]);
        } else if (arg == "--output-dir" && i + 1 < argc) {
            outputDir = argv[++i];
        } else if (arg == "--replicates" && i + 1 < argc) {
            numReplicates = atoi(argv[++i]);
        } else if (arg == "--cells" && i + 1 < argc) {
            cellNum = atoi(argv[++i]);
        } else if (arg == "--no-checkpoint") {
            checkpointEnabled = false;
        } else if (arg == "--repair-scale" && i + 1 < argc) {
            repairRateScale = atof(argv[++i]);
        } else if (arg == "--alpha" && i + 1 < argc) {
            alpha_LQ = atof(argv[++i]);
        } else if (arg == "--beta" && i + 1 < argc) {
            beta_LQ = atof(argv[++i]);
        } else if (arg == "--cell-line" && i + 1 < argc) {
            cellLine = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            cout << "Usage: " << argv[0] << " [OPTIONS]" << endl;
            cout << "Options:" << endl;
            cout << "  --assay-time HOURS: Assay time window in hours (default: 27.8)" << endl;
            cout << "  --output-dir PATH: Output directory (default: ./time_sensitivity_output)" << endl;
            cout << "  --replicates N: Replicates per dose (default: 5)" << endl;
            cout << "  --cells N: Number of cells (default: 1000)" << endl;
            cout << "  --no-checkpoint: Disable checkpoint" << endl;
            cout << "  --repair-scale VALUE: Repair rate scale factor (default: 0.1)" << endl;
            cout << "  --alpha VALUE: LQ alpha parameter (default: 0.34 for HS-23)" << endl;
            cout << "  --beta VALUE: LQ beta parameter (default: 0.028 for HS-23)" << endl;
            cout << "  --cell-line NAME: Cell line name (default: HS-23)" << endl;
            return 0;
        }
    }
    
    // Calculate time steps for display
    double T = 10;  // seconds per step
    int totalTimeStepNum = static_cast<int>(assayTimeHours * 3600.0 / T);
    
    cout << "============================================" << endl;
    cout << "  Time Window Sensitivity Test" << endl;
    cout << "  Validating Repair-Time Mechanism" << endl;
    cout << "============================================" << endl;
    cout << "Configuration:" << endl;
    cout << "  Cell line: " << cellLine << endl;
    cout << "  LQ parameters: α=" << alpha_LQ << " Gy⁻¹, β=" << beta_LQ << " Gy⁻²" << endl;
    cout << "  Assay time: " << assayTimeHours << " hours (" << totalTimeStepNum << " time steps)" << endl;
    cout << "  Checkpoint: " << (checkpointEnabled ? "enabled" : "disabled") << endl;
    cout << "  Cells: " << cellNum << endl;
    cout << "  Replicates: " << numReplicates << endl;
    cout << "  Repair rate scale: " << repairRateScale << endl;
    cout << "  Output directory: " << outputDir << endl;
    
    // Create output directory
    createDirectoryIfNotExists(outputDir.c_str());
    
    // Set up OpenMP
    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    cout << "\nRunning with " << numThreads << " OpenMP threads" << endl;
    
    double wall_start = omp_get_wtime();
    
    // Run the test
    unsigned long long baseSeed = 12345ULL;
    
    TimeSensitivityResult result = runTimeSensitivityTest(
        assayTimeHours,
        alpha_LQ,
        beta_LQ,
        cellLine,
        cellNum,
        numReplicates,
        checkpointEnabled,
        repairRateScale,
        baseSeed
    );
    
    double wall_end = omp_get_wtime();
    
    //------------------------------------------------------------------------
    // Save results
    //------------------------------------------------------------------------
    string checkpointStr = checkpointEnabled ? "on" : "off";
    string baseFilename = "time_" + to_string(static_cast<int>(assayTimeHours)) + "h_checkpoint_" + checkpointStr;
    
    // Save detailed CSV
    string csvFile = outputDir + "/" + baseFilename + ".csv";
    ofstream csv(csvFile);
    csv << "Dose,SF_LQ,SF_Calibrated,SF_Simulated,SF_Simulated_SD,AbsLog10Error" << endl;
    
    const double sf_eps = 1e-12;
    for (size_t n = 0; n < result.doses.size(); n++) {
        double log10_error = std::fabs(std::log10(result.sf_sim[n] + sf_eps) - std::log10(result.sf_lq[n] + sf_eps));
        csv << result.doses[n] << "," << result.sf_lq[n] << "," << result.sf_calibrated[n] << ","
            << result.sf_sim[n] << "," << result.sf_sim_sd[n] << "," << log10_error << endl;
    }
    csv.close();
    
    // Save summary to a single-line CSV for easy aggregation
    string summaryFile = outputDir + "/" + baseFilename + "_summary.csv";
    ofstream summary(summaryFile);
    summary << "AssayTime_h,Checkpoint,CellLine,Alpha,Beta,AvgLog10Error,MaxLog10Error,HighDoseError,LowDoseError" << endl;
    summary << assayTimeHours << "," << checkpointStr << "," << cellLine << ","
            << alpha_LQ << "," << beta_LQ << ","
            << result.avg_log10_error << "," << result.max_log10_error << ","
            << result.high_dose_avg_error << "," << result.low_dose_avg_error << endl;
    summary.close();
    
    //------------------------------------------------------------------------
    // Print results
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  RESULTS" << endl;
    cout << "============================================" << endl;
    cout << "Assay time: " << assayTimeHours << " hours" << endl;
    cout << "Checkpoint: " << checkpointStr << endl;
    cout << "\nError metrics:" << endl;
    cout << "  Overall Avg |Δlog10(SF)|: " << fixed << setprecision(4) << result.avg_log10_error << endl;
    cout << "  Max |Δlog10(SF)|:         " << result.max_log10_error << endl;
    cout << "  Low-dose (<4 Gy) error:   " << result.low_dose_avg_error << endl;
    cout << "  High-dose (>=4 Gy) error: " << result.high_dose_avg_error << endl;
    
    cout << "\nDose-by-dose results:" << endl;
    cout << setw(8) << "Dose" << setw(12) << "SF_LQ" << setw(12) << "SF_Sim" << setw(12) << "Error" << endl;
    cout << string(44, '-') << endl;
    
    for (size_t n = 0; n < result.doses.size(); n++) {
        double log10_error = std::fabs(std::log10(result.sf_sim[n] + sf_eps) - std::log10(result.sf_lq[n] + sf_eps));
        cout << setw(8) << fixed << setprecision(1) << result.doses[n] 
             << setw(12) << setprecision(4) << result.sf_lq[n]
             << setw(12) << result.sf_sim[n]
             << setw(12) << setprecision(4) << log10_error << endl;
    }
    
    cout << "\nWall time: " << fixed << setprecision(2) << (wall_end - wall_start) << " seconds" << endl;
    cout << "\nOutput files:" << endl;
    cout << "  - " << csvFile << endl;
    cout << "  - " << summaryFile << endl;
    
    return 0;
}
