/**
 * @file test_calibration_validation.cc
 * @brief Validation test for calibration methods using LQ model generated data
 * 
 * This test:
 * 1. Generates "fake" survival fraction data using LQ model (α=0.75, β=0.11)
 * 2. Calibrates the Repair Mediated model (S1→S2→S3) to this data
 * 3. Runs the full cell state simulation with calibrated parameters
 * 4. Compares simulated SF with LQ-generated data to validate calibration
 * 
 * @author Generated for calibration validation
 * @date 2025
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

#include "Cell.hh"
#include "CellStateModel.hh"
#include "CellLayoutInitializer.hh"

// Include calibration headers (with misrepair)
#include "RepairMediatedCalibrator_with_k_error.hh"
#include "DataTypes.hh"

using namespace std;
using namespace CellStateCalibration;

//============================================================================
// LQ Model for generating fake survival data
//============================================================================
double LQ_survival(double dose, double alpha_LQ, double beta_LQ) {
    return exp(-alpha_LQ * dose - beta_LQ * dose * dose);
}

//============================================================================
// Class to read radiation transport results (copied from test_phase_state_dose_parallel.cc)
// This reads dose and DSB data from Monte Carlo simulation output files
//============================================================================
class RadiationTransportData {
private:
    std::map<int, double> cellDoseMap;      // cellID -> dose
    std::map<int, double> cellDoseStdMap;   // cellID -> dose std
    std::map<int, double> cellDoseFractionMap; // cellID -> dose fraction relative to mean
    double singleCellDose;      // reference single cell dose
    double singleCellDoseStd;
    double singleCellDSB;       // reference single cell DSB
    double singleCellDSBStd;
    double meanDose;
    bool dataLoaded;

public:
    RadiationTransportData() : singleCellDose(0), singleCellDoseStd(0), 
                                singleCellDSB(0), singleCellDSBStd(0), 
                                meanDose(0), dataLoaded(false) {}
    
    // Read dose data from MC simulation output
    void ReadDoseTallyOutput(const string& doseFileName) {
        ifstream file(doseFileName.c_str());
        if (!file.is_open()) {
            cerr << "Warning: Could not open dose file: " << doseFileName << endl;
            return;
        }
        
        string line;
        int lineNum = 0;
        
        while (getline(file, line)) {
            if (lineNum >= 1) { // Skip header
                istringstream ss(line);
                string field;
                int elementNum = 0;
                int cellID = 0;
                double totalDose = 0, doseStd = 0;
                
                while (getline(ss, field, ',')) {
                    if (elementNum == 0) cellID = stoi(field);
                    else if (elementNum == 1) totalDose = stod(field);
                    else if (elementNum == 2) doseStd = stod(field);
                    elementNum++;
                }
                
                cellDoseMap[cellID] = totalDose;
                cellDoseStdMap[cellID] = doseStd;
            }
            lineNum++;
        }
        file.close();
        
        // Calculate mean dose
        if (!cellDoseMap.empty()) {
            double sumDose = 0;
            for (auto& cell : cellDoseMap) {
                sumDose += cell.second;
            }
            meanDose = sumDose / cellDoseMap.size();
            
            // Calculate dose fraction for each cell
            for (auto& cell : cellDoseMap) {
                cellDoseFractionMap[cell.first] = cell.second / meanDose;
            }
        }
        
        cout << "  Loaded dose data for " << cellDoseMap.size() << " cells, mean dose = " << meanDose << endl;
    }
    
    // Read single cell reference data (for DSB calculation)
    void ReadSingleCellDoseAsReference(const string& doseFileName) {
        ifstream file(doseFileName.c_str());
        if (!file.is_open()) {
            cerr << "Warning: Could not open single cell dose file: " << doseFileName << endl;
            return;
        }
        
        string line;
        int lineNum = 0;
        
        while (getline(file, line)) {
            if (lineNum >= 1) { // Skip header
                istringstream ss(line);
                string field;
                int elementNum = 0;
                
                while (getline(ss, field, ',')) {
                    if (elementNum == 1) singleCellDose = stod(field);
                    else if (elementNum == 2) singleCellDoseStd = stod(field);
                    elementNum++;
                }
                break; // Only read first data line
            }
            lineNum++;
        }
        file.close();
        cout << "  Single cell reference dose = " << singleCellDose << endl;
    }
    
    // Read single cell DNA damage reference
    void ReadSingleCellDNADamageAsReference(const string& dnaFileName) {
        ifstream file(dnaFileName.c_str());
        if (!file.is_open()) {
            cerr << "Warning: Could not open single cell DNA damage file: " << dnaFileName << endl;
            return;
        }
        
        string line;
        int lineNum = 0;
        
        while (getline(file, line)) {
            if (lineNum >= 1) { // Skip header
                istringstream ss(line);
                string field;
                int elementNum = 0;
                
                while (getline(ss, field, ',')) {
                    if (elementNum == 3) singleCellDSB = stod(field);  // DSB column
                    else if (elementNum == 4) singleCellDSBStd = stod(field);
                    elementNum++;
                }
                break;
            }
            lineNum++;
        }
        file.close();
        cout << "  Single cell reference DSB = " << singleCellDSB << endl;
        dataLoaded = true;
    }
    
    // Get absolute DSB for a specific cell at a given prescribed dose
    // Formula: absDSB = cellDoseFraction * prescribedDose / singleCellDose * singleCellDSB
    double GetAbsDSBOfCell(int cellID, double prescribedDose) const {
        if (!dataLoaded || cellDoseFractionMap.find(cellID) == cellDoseFractionMap.end()) {
            // If no data, return uniform DSB
            return prescribedDose * 40.0;  // Default: 40 DSB per Gy
        }
        
        double doseFraction = cellDoseFractionMap.at(cellID);
        double absDSB = doseFraction * prescribedDose / singleCellDose * singleCellDSB;
        return absDSB;
    }
    
    // Get dose fraction for a cell (for debugging)
    double GetDoseFractionOfCell(int cellID) const {
        if (cellDoseFractionMap.find(cellID) == cellDoseFractionMap.end()) {
            return 1.0;
        }
        return cellDoseFractionMap.at(cellID);
    }
    
    bool IsDataLoaded() const { return dataLoaded; }
    int GetCellCount() const { return cellDoseMap.size(); }
};

//============================================================================
// Structure to store results
//============================================================================
struct DoseResult {
    double dose;
    int G1_num, S_num, G2_num, M_num, G0_num;
    int S1_num, S2_num, S3_num;
    int totalCells;
};

struct CellSFResult {
    double dose;
    int colonyNumber;
    double sf;
};

struct CheckpointDiagnostics {
    double dose;
    int checkpointHoldCount;
    double checkpointHoldPercent;
    int S1_FromCheckpoint;
    int S2_FromCheckpoint;
    int S3_FromCheckpoint;
    double avgHoldAge_h;
    double maxHoldAge_h;
    double sf_simulated;
    double sf_lq;
    double overpredictionFactor;
};

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

// Helper function to create directory
void createDirectoryIfNotExists(const char* path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0755);
        cout << "Created directory: " << path << endl;
    }
}

int main(int argc, char** argv)
{
    // Parse command-line arguments for checkpoint and soft saturation enable/disable
    bool checkpointEnabled = true;  // Default: enabled
    bool softSaturationEnabled = true;  // Default: enabled (soft saturation)
    double repairRateScale = 0.1;  // Default: 0.1 (10x slower repair)
    int numReplicates = 5;  // Default: match historical plots (5 replicates/dose)
    string outputSuffix = "";
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        string arg = string(argv[i]);
        
        // Checkpoint options
        if (arg == "--no-checkpoint" || arg == "-n") {
            checkpointEnabled = false;
            outputSuffix = "_no_checkpoint";
        } else if (arg == "--checkpoint" || arg == "-c") {
            checkpointEnabled = true;
            outputSuffix = "_with_checkpoint";
        }
        // Soft saturation options
        else if (arg == "--no-soft-saturation" || arg == "-ns") {
            softSaturationEnabled = false;
            outputSuffix = "_hard_saturation";
        } else if (arg == "--soft-saturation" || arg == "-s") {
            softSaturationEnabled = true;
            outputSuffix = "_soft_saturation";
        }
        // Replicates per dose
        else if (arg == "--replicates" && i + 1 < argc) {
            numReplicates = atoi(argv[++i]);
            if (numReplicates <= 0) {
                cerr << "Error: replicates must be > 0" << endl;
                return 1;
            }
        }
        // Repair rate scale option
        else if (arg == "--repair-scale" && i + 1 < argc) {
            repairRateScale = atof(argv[++i]);
            if (repairRateScale <= 0) {
                cerr << "Error: repair scale must be > 0" << endl;
                return 1;
            }
        }
        // Help
        else if (arg == "--help" || arg == "-h") {
            cout << "Usage: " << argv[0] << " [OPTIONS]" << endl;
            cout << "Options:" << endl;
            cout << "  --checkpoint, -c: Enable probabilistic checkpoint (default)" << endl;
            cout << "  --no-checkpoint, -n: Disable checkpoint" << endl;
            cout << "  --soft-saturation, -s: Enable soft saturation (default)" << endl;
            cout << "  --no-soft-saturation, -ns: Disable soft saturation (use hard saturation)" << endl;
            cout << "  --replicates N: Replicates per dose (default: 5)" << endl;
            cout << "  --repair-scale VALUE: Set repair rate scale factor (default: 0.1)" << endl;
            return 0;
        } else {
            cout << "Unknown option: " << arg << endl;
            cout << "Use --help for usage information" << endl;
            return 1;
        }
    }
    
    // Display configuration
    cout << "============================================" << endl;
    cout << "  Calibration Validation Test" << endl;
    cout << "  LQ Model → Calibration → Simulation" << endl;
    cout << "  CHECKPOINT: " << (checkpointEnabled ? "ENABLED (Probabilistic)" : "DISABLED") << endl;
    cout << "  SOFT SATURATION: " << (softSaturationEnabled ? "ENABLED" : "DISABLED") << endl;
    cout << "  REPAIR RATE SCALE: " << repairRateScale << "x" << endl;
    cout << "============================================" << endl;

    //========================================================================
    // STEP 1: Generate fake survival data using LQ model
    //========================================================================
    cout << "\n=== STEP 1: Generate Fake LQ Survival Data ===" << endl;
    
    // LQ model parameters (given)
    double alpha_LQ = 0.75;  // Gy^-1
    double beta_LQ = 0.11;   // Gy^-2
    double alpha_beta_ratio = alpha_LQ / beta_LQ;
    
    cout << "LQ Model Parameters:" << endl;
    cout << "  alpha = " << alpha_LQ << " Gy^-1" << endl;
    cout << "  beta  = " << beta_LQ << " Gy^-2" << endl;
    cout << "  alpha/beta = " << alpha_beta_ratio << " Gy" << endl;
    
    // Generate fake data at specific dose points
    vector<double> fakeDoses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
    vector<DataPoint> fakeData;
    
    cout << "\nGenerated LQ Survival Data:" << endl;
    cout << "  Dose (Gy)    SF_LQ" << endl;
    cout << "  ----------------------" << endl;
    
    for (double dose : fakeDoses) {
        double sf = LQ_survival(dose, alpha_LQ, beta_LQ);
        fakeData.push_back(DataPoint(dose, sf, 0.0));
        cout << "  " << setw(8) << fixed << setprecision(2) << dose 
             << "    " << setprecision(6) << sf << endl;
    }
    
    //========================================================================
    // STEP 2: Calibrate Repair Mediated Model WITH MISREPAIR to fake data
    //========================================================================
    cout << "\n=== STEP 2: Calibrate Repair Mediated Model WITH MISREPAIR ===" << endl;
    
    // Cell cycle time (needed for fixed T21/T23)
    const double T_cellCycle = 10.0;  // hours (tG1 + tS + tG2 + tM = 1.5 + 6 + 1.5 + 1)
    
    // Configuration for calibration
    FitConfigS2 cfg;
    cfg.kappa = 40.0;           // DSB/Gy
    cfg.dose_max_fit = 6.0;     // Fit all data
    cfg.T_assay_h = 27.8; // Match simulation length (~27.8 hours)
    
    // FIX T21 and T23 at T_cellCycle (using new API)
    cfg.setTimescalesToCellCycle(T_cellCycle);
    
    cout << "Calibration Configuration:" << endl;
    cout << "  kappa (DSB/Gy): " << cfg.kappa << endl;
    cout << "  dose_max_fit: " << cfg.dose_max_fit << " Gy" << endl;
    cout << "  T_assay: " << cfg.T_assay_h << " hours" << endl;
    cout << "  T21_fixed: " << cfg.T21_fixed << " hours (= T_cellCycle)" << endl;
    cout << "  T23_fixed: " << cfg.T23_fixed << " hours (= T_cellCycle)" << endl;
    cout << "  Model: Repair Mediated WITH stochastic misrepair (k_error)" << endl;
    
    // Create misrepair calibrator and fit
    RepairMediatedMisrepairCalibrator calibrator(cfg);
    calibrator.setData(fakeData);
    calibrator.setDefaultUncertainties(0.1);  // 10% uncertainty
    
    // Initial guess for 4 parameters: e2, e3, a, k_error (T21, T23 are FIXED)
    ParamsS2 x0;
    x0.e2 = 2.0;
    x0.e3 = 10.0;
    x0.a = 0.05;
    x0.T21 = T_cellCycle;  // Will use fixed value
    x0.T23 = T_cellCycle;  // Will use fixed value
    x0.k_error = 1e-4;      // per (DSB*hour), within plausible range
    
    ParamsS2 step;
    step.e2 = 0.5;
    step.e3 = 2.0;
    step.a = 0.01;
    step.T21 = 1.0;        // Doesn't matter (fixed)
    step.T23 = 1.0;        // Doesn't matter (fixed)
    step.k_error = 5e-5;
    
    cout << "\nFitting 4 parameters: e2, e3, a, k_error (T21, T23 FIXED at T_cellCycle)" << endl;
    cout << "Initial guess: e2=" << x0.e2 << ", e3=" << x0.e3 << ", a=" << x0.a 
         << ", k_error=" << x0.k_error << endl;
    RepairMediatedMisrepairFitResult result = calibrator.fit(x0, step);
    
    cout << "\n=== Calibration Results ===" << endl;
    cout << "4 Parameters OPTIMIZED from survival data:" << endl;
    cout << "  e2  = " << fixed << setprecision(4) << result.params.e2 << " (E2/sigma)" << endl;
    cout << "  e3  = " << result.params.e3 << " (E3/sigma)" << endl;
    cout << "  a   = " << result.params.a << " (alpha/sigma)" << endl;
    cout << "  k_error = " << result.params.k_error << " per (DSB*hour) - MISREPAIR RATE" << endl;
    cout << "\n2 Parameters FIXED (not optimized):" << endl;
    cout << "  T21 = " << result.params.T21 << " hours (repair timescale) - FIXED at T_cellCycle" << endl;
    cout << "  T23 = " << result.params.T23 << " hours (death timescale) - FIXED at T_cellCycle" << endl;
    cout << "\nOptimization:" << endl;
    cout << "  NLL = " << result.negLogLikelihood << endl;
    cout << "  Converged: " << (result.converged ? "Yes" : "No") << endl;
    
    // Convert to physical parameters (choose sigma = 10.0)
    double sigma_chosen = 10.0;
    CellStateModelParamsS2 calibratedParams = CellStateModelParamsS2::fromReduced(
        result.params, sigma_chosen, cfg.kappa);
    // Set k_error (not included in fromReduced)
    calibratedParams.k_error = result.params.k_error;
    
    cout << "\nPhysical Parameters (sigma = " << sigma_chosen << "):" << endl;
    calibratedParams.print();
    
    //========================================================================
    // STEP 3: Run Cell State Simulation with Calibrated Parameters
    //========================================================================
    cout << "\n=== STEP 3: Run Cell State Simulation ===" << endl;
    
    //------------------------------------------------------------------------
    // Simulation parameters
    //------------------------------------------------------------------------
    double xDim = 5.0;      // mm
    double yDim = 5.0;      // mm
    double zDim = 0.0;      // mm (2D)
    double d = 0.05;        // grid size mm
    
    // Time parameters
    double deltaT_cellPhaseUpdate = 60;  // seconds
    double deltaT_cellStateUpdate = 60;  // seconds
    double T = 10;                       // time step seconds
    int totalTimeStepNum = 10000;        // ~27.8 hours
    int cellNum = 1000;
    
    // Cell cycle parameters
    double tG1 = 1.5, sigmaG1 = 0.25;
    double tS = 6.0, sigmaS = 0.25;
    double tG2 = 1.5, sigmaG2 = 0.25;
    double tM = 1.0, sigmaM = 0.25;
    // T_cellCycle already defined earlier (= 10.0 hours)
    
    // Use CALIBRATED parameters for cell state model
    double alpha_model = calibratedParams.alpha;
    double E1 = calibratedParams.E1;
    double E2 = calibratedParams.E2;
    double E3 = calibratedParams.E3;
    double sigma_model = calibratedParams.sigma;
    
    // Other parameters (not from calibration)
    double beta_model = 25.71;  // bystander parameter
    double fG1 = 1.0, fG2 = 1.015306122, fM = 1.015306122, fS = 0.816326;
    double f1 = 0.62, f2 = 0.38;
    double lambda1 = 3.31, lambda2 = 0.14;
    
    // Repair rate scaling factor: < 1.0 slows down repair, > 1.0 speeds it up
    // Testing slower repair to see if it reduces low-dose SF toward LQ
    // repairRateScale is now set from command-line argument (default: 0.1)
    double lambda1_scaled = lambda1 * repairRateScale;
    double lambda2_scaled = lambda2 * repairRateScale;
    
    // Calculate effective repair time (mean repair time weighted by fractions)
    // τ_eff = f1 * (1/λ1) + f2 * (1/λ2)
    double tau1 = 1.0 / lambda1_scaled;  // Mean repair time for fast component
    double tau2 = 1.0 / lambda2_scaled;  // Mean repair time for slow component
    double effective_repair_time = f1 * tau1 + f2 * tau2;  // Weighted mean repair time
    cout << "\nEffective Repair Time Calculation:" << endl;
    cout << "  Fast component (f1=" << f1 << "): τ1 = 1/λ1 = " << tau1 << " hours" << endl;
    cout << "  Slow component (f2=" << f2 << "): τ2 = 1/λ2 = " << tau2 << " hours" << endl;
    cout << "  Effective repair time: τ_eff = " << effective_repair_time << " hours" << endl;
    
    // DSB per Gy (from calibration config)
    double DSB_per_Gy = cfg.kappa;
    
    // Dose levels for simulation
    vector<double> simDoses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};

	    // Replicates per dose (CLI-configurable)
	    const unsigned long long baseSeed = 1337ULL;
    
    int numThreads = 10;
    omp_set_num_threads(numThreads);
    
    cout << "\nSimulation Parameters:" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " mm" << endl;
    cout << "  Cell number: " << cellNum << endl;
    cout << "  Time steps: " << totalTimeStepNum << endl;
    cout << "  DSB per Gy: " << DSB_per_Gy << endl;
    cout << "  OpenMP threads: " << numThreads << endl;
	    cout << "  Replicates per dose: " << numReplicates << endl;
    cout << "  Repair rate scale: " << repairRateScale << "x (slower repair test)" << endl;
    cout << "    Original: λ1=" << lambda1 << " 1/h, λ2=" << lambda2 << " 1/h" << endl;
    cout << "    Scaled:  λ1=" << lambda1_scaled << " 1/h, λ2=" << lambda2_scaled << " 1/h" << endl;
    
    cout << "\nCalibrated Cell State Parameters:" << endl;
    cout << "  E1 = " << E1 << endl;
    cout << "  E2 = " << E2 << endl;
    cout << "  E3 = " << E3 << endl;
    cout << "  sigma = " << sigma_model << endl;
    cout << "  alpha = " << alpha_model << endl;
    cout << "  Misrepair enabled: simulation uses calibrated k_error" << endl;
    
    //------------------------------------------------------------------------
    // Create Cell Layout
    //------------------------------------------------------------------------
    CellLayoutInitializer layout;
    layout.SetRandomSeed(baseSeed);
    layout.SetCellHomeParamter(d, d, d);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);
    cout << "\nTotal cells: " << layout.GetCellNumber() << endl;
    
    //------------------------------------------------------------------------
    // Load radiation transport data (like test_phase_state_dose_parallel.cc)
    // This provides cell-specific DSB values from Monte Carlo simulation
    //------------------------------------------------------------------------
    cout << "\n--- Loading Radiation Transport Data ---" << endl;
    RadiationTransportData radiationData;
    
    // File paths (same as test_phase_state_dose_parallel.cc)
    string doseFileName = "/home/user/Geant4Projects/CellStateTransitionPaper/RAD-build/Mono_Electron_Plane_Eng_1PNum_1000000_dose.csv";
    string doseFileName_singleCell = "/home/user/Geant4Projects/CellStateTransitionPaper/RAD-build/Mono_Electron_Eng_1PNum_100000_dose.csv";
    string DNADamageFileName_singleCell = "/home/user/Geant4Projects/CellStateTransitionPaper/RAD-build/Mono_Electron_Eng_1PNum_100000_DNADamage.csv";
    
    radiationData.ReadDoseTallyOutput(doseFileName);
    radiationData.ReadSingleCellDoseAsReference(doseFileName_singleCell);
    radiationData.ReadSingleCellDNADamageAsReference(DNADamageFileName_singleCell);
    
    if (radiationData.IsDataLoaded()) {
        cout << "  Using cell-specific DSB from MC simulation (like old version)" << endl;
    } else {
        cout << "  WARNING: MC data not loaded, using uniform DSB = dose × " << DSB_per_Gy << endl;
    }
    
    //------------------------------------------------------------------------
    // Create output directory
    //------------------------------------------------------------------------
    // Create output directory with suffix
    string outputDir = "./calibration_validation_output" + outputSuffix;
    if (repairRateScale != 0.1) {  // Add repair scale suffix if not default
        string scaleStr = to_string(repairRateScale);
        // Remove trailing zeros
        scaleStr.erase(scaleStr.find_last_not_of('0') + 1, string::npos);
        scaleStr.erase(scaleStr.find_last_not_of('.') + 1, string::npos);
        outputDir += "_repair" + scaleStr;
    }
    createDirectoryIfNotExists(outputDir.c_str());
    
    //------------------------------------------------------------------------
    // Prepare results storage
    //------------------------------------------------------------------------
    vector<DoseResult> results(simDoses.size());
    vector<CellSFResult> cellSF_results(simDoses.size());
    vector<double> cellSF_sd(simDoses.size(), 0.0);
    vector<CheckpointDiagnostics> checkpointDiag(simDoses.size());
    
    //------------------------------------------------------------------------
    // Main parallel dose loop
    //------------------------------------------------------------------------
    cout << "\n--- Starting Simulation ---" << endl;
    double wall_start = omp_get_wtime();
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t n = 0; n < simDoses.size(); n++)
    {
        double dose = simDoses[n];
        
        #pragma omp critical
        {
            cout << "Thread " << omp_get_thread_num() << " processing dose " << dose << " Gy" << endl;
        }
        
        // Precompute initial DSB for the initial cell population (same for all replicates)
        std::vector<int> initialDSB(cellNum, 0);
        std::vector<double> initialDSB_double(cellNum, 0.0);
        for (int cid = 0; cid < cellNum; cid++) {
            double dsb_val = radiationData.GetAbsDSBOfCell(cid, dose);
            initialDSB_double[cid] = dsb_val;
            initialDSB[cid] = static_cast<int>(dsb_val);
        }
        
        // Diagnostic for 4 Gy: analyze initial DSB distribution
        bool is4Gy = (std::abs(dose - 4.0) < 0.01);
        if (is4Gy) {
            std::vector<double> dsb_sorted = initialDSB_double;
            std::sort(dsb_sorted.begin(), dsb_sorted.end());
            double dsb_median = dsb_sorted[dsb_sorted.size() / 2];
            double dsb_95th = dsb_sorted[static_cast<size_t>(dsb_sorted.size() * 0.95)];
            double dsb_max = dsb_sorted.back();
            double dsb_mean = 0.0;
            for (double d : initialDSB_double) dsb_mean += d;
            dsb_mean /= cellNum;
            
            // Calculate death threshold: E3/alpha (in DSB units)
            // E3 is in energy units, alpha converts DSB to energy: E = alpha * N_DSB
            // So N_DSB_death = E3 / alpha
            double death_threshold_DSB = (alpha_model > 0.0) ? (E3 / alpha_model) : 0.0;
            int cells_above_threshold = 0;
            for (double d : initialDSB_double) {
                if (d >= death_threshold_DSB) cells_above_threshold++;
            }
            
            cout << "\n=== 4 Gy DIAGNOSTIC: Initial DSB Distribution ===" << endl;
            cout << "  Mean DSB: " << dsb_mean << endl;
            cout << "  Median DSB: " << dsb_median << endl;
            cout << "  95th percentile DSB: " << dsb_95th << endl;
            cout << "  Max DSB: " << dsb_max << endl;
            cout << "  Death threshold (E3/alpha): " << death_threshold_DSB << " DSB" << endl;
            cout << "  Cells above death threshold: " << cells_above_threshold 
                 << " / " << cellNum << " (" << (100.0 * cells_above_threshold / cellNum) << "%)" << endl;
        }

        RunningStats sfStats;
        double colonyNumberSum = 0.0;
        
        // Accumulate checkpoint diagnostics across replicates
        RunningStats checkpointHoldCountStats;
        RunningStats checkpointHoldPercentStats;
        RunningStats S1_FromCheckpointStats;
        RunningStats S2_FromCheckpointStats;
        RunningStats S3_FromCheckpointStats;
        RunningStats avgHoldAgeStats;
        RunningStats maxHoldAgeStats;

        // Track and average fractions over replicates at each time point
        vector<double> trackTimes;
        vector<double> sumF1, sumF2, sumF3;
        vector<double> sumC1, sumC2, sumC3;
        vector<double> sumCheckpointHold;  // Track checkpoint hold count over time

        for (int rep = 0; rep < numReplicates; rep++)
        {
            // Reset clocks
            double clock_cellPhaseUpdate = 0;
            double clock_cellStateUpdate = 0;
            int cycle_cellPhaseUpdate = 0;
            int cycle_cellStateUpdate = 0;

            // Create thread-local Cell object with CALIBRATED parameters
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
            statePara.To21 = effective_repair_time;  // Use effective repair time for S2→S1 transition
            // Option 1: Increase To23 proportionally with To21 to slow death progression
            // This keeps repair and death timescales balanced
            statePara.To23 = effective_repair_time;  // Match To21 to slow S2→S3 death

            CellDNADamageRepairParameter dnaRepairPara;
            dnaRepairPara.f1 = f1; dnaRepairPara.f2 = f2;
            dnaRepairPara.lambda1 = lambda1_scaled; dnaRepairPara.lambda2 = lambda2_scaled;

            CellBystanderSignalParameter bystanderPara;
            bystanderPara.Rt = 5000; bystanderPara.mu = 200;

            testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);

	            // Initialize CellStateModel (seeded for this replicate)
	            CellStateModel myCellState;
	            // Deterministic seeding: do NOT include omp thread id.
	            // Using omp_get_thread_num() makes results depend on OpenMP scheduling (schedule(dynamic)),
	            // so repeated runs can produce different SF curves even with the same inputs.
	            myCellState.SetRandomSeed(baseSeed + static_cast<unsigned long long>(n) * 100000ULL +
	                                      static_cast<unsigned long long>(rep) * 1000ULL);
	            myCellState.TissueGeometryInitialization(xDim, yDim, zDim, d);
		            myCellState.CellStateModelParameterSetup(testCell);
		            myCellState.SetMisrepairRate("Epithelia", calibratedParams.k_error);
		            myCellState.SetUpContactInhibition(true);
	            myCellState.SetCheckpointEnabled(checkpointEnabled);  // Set checkpoint enable/disable
	            myCellState.SetSoftSaturationEnabled(softSaturationEnabled);  // Set soft saturation enable/disable

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

            // Track state distribution at specific time points
            vector<double> repTimes;
            vector<int> repS1, repS2, repS3;
            vector<int> repCheckpointHold;  // Track checkpoint hold count at each time point
            double trackInterval = 3600.0;
            double nextTrackTime = 0.0;
            double currentTime = 0.0;

            for (int i = 0; i < totalTimeStepNum; i++)
            {
                clock_cellStateUpdate += T;
                clock_cellPhaseUpdate += T;
                currentTime += T;

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

                if (currentTime >= nextTrackTime) {
                    map<int, string> stateMap = myCellState.GetCellState();
                    int s1 = 0, s2 = 0, s3 = 0;
                    for (auto& cell : stateMap) {
                        if (cell.second == "S1") s1++;
                        else if (cell.second == "S2") s2++;
                        else if (cell.second == "S3") s3++;
                    }
                    // Get checkpoint hold count at this time point
                    int checkpointHoldCount = myCellState.GetCheckpointHoldCount();
                    
                    repTimes.push_back(currentTime / 3600.0);
                    repS1.push_back(s1);
                    repS2.push_back(s2);
                    repS3.push_back(s3);
                    repCheckpointHold.push_back(checkpointHoldCount);
                    nextTrackTime += trackInterval;
                }
            }

            if (rep == 0) {
                trackTimes = repTimes;
                sumF1.assign(trackTimes.size(), 0.0);
                sumF2.assign(trackTimes.size(), 0.0);
                sumF3.assign(trackTimes.size(), 0.0);
                sumC1.assign(trackTimes.size(), 0.0);
                sumC2.assign(trackTimes.size(), 0.0);
                sumC3.assign(trackTimes.size(), 0.0);
                sumCheckpointHold.assign(trackTimes.size(), 0.0);
            }

            const size_t tN = std::min(trackTimes.size(), repTimes.size());
            for (size_t t = 0; t < tN; t++) {
                const int total = repS1[t] + repS2[t] + repS3[t];
                const double f1 = (total > 0) ? static_cast<double>(repS1[t]) / total : 0.0;
                const double f2 = (total > 0) ? static_cast<double>(repS2[t]) / total : 0.0;
                const double f3 = (total > 0) ? static_cast<double>(repS3[t]) / total : 0.0;
                sumF1[t] += f1; sumF2[t] += f2; sumF3[t] += f3;
                sumC1[t] += repS1[t]; sumC2[t] += repS2[t]; sumC3[t] += repS3[t];
                if (t < repCheckpointHold.size()) {
                    sumCheckpointHold[t] += repCheckpointHold[t];
                }
            }

            // Colony formation assay (snapshot proxy) for this replicate
            map<int, string> finalStateMap = myCellState.GetCellState();
            map<int, int> cellAncestryIDMap = myCellState.GetCellAncestryID();
            map<int, int> cellColonySizeMap;
            for (int i = 0; i < cellNum; i++) cellColonySizeMap[i] = 0;

            // Count final states
            int finalS1 = 0, finalS2 = 0, finalS3 = 0;
            for (const auto& cell : finalStateMap) {
                if (cell.second == "S1") finalS1++;
                else if (cell.second == "S2") finalS2++;
                else if (cell.second == "S3") finalS3++;
            }
            int finalTotal = finalS1 + finalS2 + finalS3;

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
            colonyNumberSum += colonyNumber;
            
            // Checkpoint hold diagnostics for all doses
            map<int, double> checkpointHoldAgeMap = myCellState.GetCheckpointHoldAge();
            int checkpointHoldCount = myCellState.GetCheckpointHoldCount();
            double checkpointHoldPercent = (cellNum > 0) ? (100.0 * checkpointHoldCount / cellNum) : 0.0;
            
            // Count cells in checkpoint hold by state
            int checkpointHold_S1 = 0, checkpointHold_S2 = 0, checkpointHold_S3 = 0;
            double avgHoldAge = 0.0, maxHoldAge = 0.0;
            if (checkpointHoldCount > 0) {
                for (const auto& hold : checkpointHoldAgeMap) {
                    int cid = hold.first;
                    double age = hold.second;
                    avgHoldAge += age;
                    maxHoldAge = std::max(maxHoldAge, age);
                    
                    auto itState = finalStateMap.find(cid);
                    if (itState != finalStateMap.end()) {
                        if (itState->second == "S1") checkpointHold_S1++;
                        else if (itState->second == "S2") checkpointHold_S2++;
                        else if (itState->second == "S3") checkpointHold_S3++;
                    }
                }
                avgHoldAge /= checkpointHoldCount;
            }
            
            // Accumulate checkpoint diagnostics
            checkpointHoldCountStats.add(static_cast<double>(checkpointHoldCount));
            checkpointHoldPercentStats.add(checkpointHoldPercent);
            S1_FromCheckpointStats.add(static_cast<double>(checkpointHold_S1));
            S2_FromCheckpointStats.add(static_cast<double>(checkpointHold_S2));
            S3_FromCheckpointStats.add(static_cast<double>(checkpointHold_S3));
            avgHoldAgeStats.add(avgHoldAge);
            maxHoldAgeStats.add(maxHoldAge);
            
            // Diagnostic for 4 Gy: final state breakdown
            if (is4Gy && rep == 0) {
                cout << "\n=== 4 Gy DIAGNOSTIC: Final State Distribution ===" << endl;
                cout << "  Final S1 count: " << finalS1 << " (" << (100.0 * finalS1 / cellNum) << "%)" << endl;
                cout << "  Final S2 count: " << finalS2 << " (" << (100.0 * finalS2 / cellNum) << "%)" << endl;
                cout << "  Final S3 count: " << finalS3 << " (" << (100.0 * finalS3 / cellNum) << "%)" << endl;
                cout << "  Final total cells: " << finalTotal << " / " << cellNum << endl;
                cout << "  Colony number (S1 only, size >= " << colonySizeThreshold << "): " << colonyNumber << endl;
                cout << "  Survival fraction: " << survivalFraction << endl;
                cout << "  NOTE: S2 cells are NON-CLONOGENIC (not counted in SF)" << endl;
                if (finalS2 > 0) {
                    cout << "  WARNING: " << finalS2 << " cells remain in S2 (repair state) at 27.8h" << endl;
                    cout << "    These cells are alive but non-clonogenic, contributing to SF mismatch" << endl;
                }
                cout << "\n=== 4 Gy DIAGNOSTIC: Checkpoint Hold ===" << endl;
                cout << "  Cells that entered checkpoint hold: " << checkpointHoldCount << endl;
                cout << "  Final checkpoint hold distribution: S1=" << checkpointHold_S1 
                     << ", S2=" << checkpointHold_S2 << ", S3=" << checkpointHold_S3 << endl;
                cout << "  Average hold age: " << fixed << setprecision(2) << avgHoldAge << " hours" << endl;
                cout << "  Max hold age: " << maxHoldAge << " hours" << endl;
            }
            
            // Checkpoint diagnostics for doses >= 2.5 Gy
            if (dose >= 2.5 && rep == 0) {
                cout << "\n=== CHECKPOINT DIAGNOSTIC: Dose " << dose << " Gy ===" << endl;
                cout << "  Cells that entered checkpoint hold: " << checkpointHoldCount 
                     << " (" << (100.0 * checkpointHoldCount / cellNum) << "% of initial cells)" << endl;
                cout << "  Final checkpoint hold: S1=" << checkpointHold_S1 
                     << ", S2=" << checkpointHold_S2 << ", S3=" << checkpointHold_S3 << endl;
                cout << "  Average hold age: " << fixed << setprecision(2) << avgHoldAge << " hours" << endl;
                cout << "  Max hold age: " << maxHoldAge << " hours" << endl;
                cout << "  SF: " << survivalFraction << " (LQ target: " 
                     << LQ_survival(dose, alpha_LQ, beta_LQ) << ")" << endl;
            }
        }

        // Save averaged state distribution over time for this dose
        {
            ostringstream fname;
            fname << outputDir << "/state_evolution_dose_"
                  << fixed << setprecision(1) << dose << "Gy.csv";
            ofstream fout(fname.str());
            fout << "Time_hours,S1_count,S2_count,S3_count,S1_fraction,S2_fraction,S3_fraction,CheckpointHold_count,CheckpointHold_fraction" << endl;
            for (size_t t = 0; t < trackTimes.size(); t++) {
                const double c1 = sumC1[t] / numReplicates;
                const double c2 = sumC2[t] / numReplicates;
                const double c3 = sumC3[t] / numReplicates;
                const double f1 = sumF1[t] / numReplicates;
                const double f2 = sumF2[t] / numReplicates;
                const double f3 = sumF3[t] / numReplicates;
                const double checkpointCount = sumCheckpointHold[t] / numReplicates;
                const double totalCells = c1 + c2 + c3;
                const double checkpointFraction = (totalCells > 0) ? checkpointCount / totalCells : 0.0;
                fout << trackTimes[t] << "," << c1 << "," << c2 << "," << c3
                     << "," << f1 << "," << f2 << "," << f3
                     << "," << checkpointCount << "," << checkpointFraction << endl;
            }
            fout.close();
        }


        // Store averaged SF results
        cellSF_results[n].dose = dose;
        cellSF_results[n].colonyNumber = static_cast<int>(std::round(colonyNumberSum / numReplicates));
        cellSF_results[n].sf = sfStats.mean();
        cellSF_sd[n] = sfStats.sd();
        
        // Store checkpoint diagnostics
        double sf_lq = LQ_survival(dose, alpha_LQ, beta_LQ);
        double sf_sim = sfStats.mean();
        double overpredictionFactor = (sf_lq > 1e-10) ? (sf_sim / sf_lq) : 0.0;
        
        checkpointDiag[n].dose = dose;
        checkpointDiag[n].checkpointHoldCount = static_cast<int>(std::round(checkpointHoldCountStats.mean()));
        checkpointDiag[n].checkpointHoldPercent = checkpointHoldPercentStats.mean();
        checkpointDiag[n].S1_FromCheckpoint = static_cast<int>(std::round(S1_FromCheckpointStats.mean()));
        checkpointDiag[n].S2_FromCheckpoint = static_cast<int>(std::round(S2_FromCheckpointStats.mean()));
        checkpointDiag[n].S3_FromCheckpoint = static_cast<int>(std::round(S3_FromCheckpointStats.mean()));
        checkpointDiag[n].avgHoldAge_h = avgHoldAgeStats.mean();
        checkpointDiag[n].maxHoldAge_h = maxHoldAgeStats.mean();
        checkpointDiag[n].sf_simulated = sf_sim;
        checkpointDiag[n].sf_lq = sf_lq;
        checkpointDiag[n].overpredictionFactor = overpredictionFactor;

        #pragma omp critical
        {
            cout << "  Dose " << dose << " Gy: SF_sim(avg)=" << fixed << setprecision(4) << cellSF_results[n].sf
                 << " ± " << setprecision(4) << cellSF_sd[n] << " (sd over " << numReplicates << " reps)" << endl;
        }
    }
    
    double wall_end = omp_get_wtime();
    
    //========================================================================
    // STEP 4: Compare Results
    //========================================================================
    cout << "\n=== STEP 4: Compare LQ Data vs Simulation ===" << endl;
    cout << "\n" << setw(10) << "Dose (Gy)"
         << setw(12) << "SF_LQ"
         << setw(12) << "SF_Calib"
         << setw(12) << "SF_Sim"
         << setw(14) << "|Δlog10(SF)|" << endl;
    cout << string(60, '-') << endl;
    
    double total_log10_error = 0.0;
    int count = 0;
    
    // Get calibrated model predictions
    vector<double> sf_calibrated = calibrator.predictSurvival(result.params, simDoses);
    
    // Log-space error uses a small epsilon to handle SF=0 cases.
    // Use half a "counting" unit as epsilon based on assay resolution.
    // (With 1000 cells and numReplicates runs, smallest nonzero SF is ~1/(cellNum*numReplicates).)
    const double sf_eps = std::max(1e-12, 0.5 / (static_cast<double>(cellNum) * static_cast<double>(numReplicates)));

    for (size_t n = 0; n < simDoses.size(); n++) {
        double dose = simDoses[n];
        double sf_lq = LQ_survival(dose, alpha_LQ, beta_LQ);
        double sf_calib = sf_calibrated[n];
        double sf_sim = cellSF_results[n].sf;
        
        const double log10_lq  = std::log10(sf_lq + sf_eps);
        const double log10_sim = std::log10(sf_sim + sf_eps);
        const double log10_error = std::fabs(log10_sim - log10_lq);

        total_log10_error += log10_error;
        count++;
        
        cout << setw(10) << fixed << setprecision(2) << dose
             << setw(12) << setprecision(4) << sf_lq
             << setw(12) << sf_calib
             << setw(12) << sf_sim
             << setw(14) << setprecision(3) << log10_error << endl;
    }
    
    const double avg_log10_error = (count > 0) ? total_log10_error / static_cast<double>(count) : 0.0;
    const double avg_factor_error = std::pow(10.0, avg_log10_error);
    cout << string(60, '-') << endl;
    cout << "Average |Δlog10(SF)|: " << fixed << setprecision(3) << avg_log10_error
         << "  (≈×" << setprecision(2) << avg_factor_error << ")" << endl;
    
    //------------------------------------------------------------------------
    // Write results to files
    //------------------------------------------------------------------------
    ofstream fileCompare, fileSF, fileCheckpoint;
    fileCompare.open((outputDir + "/comparison.csv").c_str());
    fileSF.open((outputDir + "/cellSurvival.csv").c_str());
    fileCheckpoint.open((outputDir + "/checkpoint_diagnostics.csv").c_str());
    
    fileCompare << "Dose,SF_LQ,SF_Calibrated,SF_Simulated,SF_Simulated_SD,AbsLog10Error" << endl;
    fileSF << "Dose,ColonyNumber,SurvivalFraction,SurvivalFraction_SD" << endl;
    fileCheckpoint << "Dose,CheckpointHoldCount,CheckpointHoldPercent,S1_FromCheckpoint,S2_FromCheckpoint,S3_FromCheckpoint,AvgHoldAge_h,MaxHoldAge_h,SF_Simulated,SF_LQ,OverpredictionFactor" << endl;
    
    for (size_t n = 0; n < simDoses.size(); n++) {
        double dose = simDoses[n];
        double sf_lq = LQ_survival(dose, alpha_LQ, beta_LQ);
        double sf_calib = sf_calibrated[n];
        double sf_sim = cellSF_results[n].sf;
        const double log10_lq  = std::log10(sf_lq + sf_eps);
        const double log10_sim = std::log10(sf_sim + sf_eps);
        const double log10_error = std::fabs(log10_sim - log10_lq);
        
        fileCompare << dose << "," << sf_lq << "," << sf_calib << "," << sf_sim << "," << cellSF_sd[n] << "," << log10_error << endl;
        fileSF << cellSF_results[n].dose << "," << cellSF_results[n].colonyNumber << "," << cellSF_results[n].sf << "," << cellSF_sd[n] << endl;
        
        // Write checkpoint diagnostics
        fileCheckpoint << checkpointDiag[n].dose << ","
                       << checkpointDiag[n].checkpointHoldCount << ","
                       << checkpointDiag[n].checkpointHoldPercent << ","
                       << checkpointDiag[n].S1_FromCheckpoint << ","
                       << checkpointDiag[n].S2_FromCheckpoint << ","
                       << checkpointDiag[n].S3_FromCheckpoint << ","
                       << checkpointDiag[n].avgHoldAge_h << ","
                       << checkpointDiag[n].maxHoldAge_h << ","
                       << checkpointDiag[n].sf_simulated << ","
                       << checkpointDiag[n].sf_lq << ","
                       << checkpointDiag[n].overpredictionFactor << endl;
    }
    
    fileCompare.close();
    fileSF.close();
    fileCheckpoint.close();
    
    // Save calibrated parameters
    ofstream fileParams;
    fileParams.open((outputDir + "/calibrated_params.txt").c_str());
    fileParams << "=== Calibration Validation Results ===" << endl;
    fileParams << "\nInput LQ Parameters:" << endl;
    fileParams << "  alpha_LQ = " << alpha_LQ << " Gy^-1" << endl;
    fileParams << "  beta_LQ = " << beta_LQ << " Gy^-2" << endl;
    fileParams << "\nCalibrated Reduced Parameters (fitted):" << endl;
    fileParams << "  e2 = " << result.params.e2 << endl;
    fileParams << "  e3 = " << result.params.e3 << endl;
    fileParams << "  a = " << result.params.a << endl;
    fileParams << "  k_error = " << result.params.k_error << " per (DSB*hour)" << endl;
    fileParams << "\nOptimized Timescales:" << endl;
    fileParams << "  T21 = " << result.params.T21 << " hours (from optimization)" << endl;
    fileParams << "  T23 = " << result.params.T23 << " hours (from optimization)" << endl;
    fileParams << "\nCalibrated Physical Parameters (sigma=" << sigma_chosen << "):" << endl;
    fileParams << "  E1 = " << calibratedParams.E1 << endl;
    fileParams << "  E2 = " << calibratedParams.E2 << endl;
    fileParams << "  E3 = " << calibratedParams.E3 << endl;
    fileParams << "  sigma = " << calibratedParams.sigma << endl;
    fileParams << "  alpha = " << calibratedParams.alpha << endl;
    fileParams << "  kappa = " << calibratedParams.kappa << " DSB/Gy" << endl;
    fileParams << "  k_error = " << calibratedParams.k_error << " per (DSB*hour)" << endl;
    fileParams << "\nProbabilistic Checkpoint Parameters (Saturating T_hold Model):" << endl;
    fileParams << "  mu_min_hold = 1.386 (log-mean for min hold, median = 4h)" << endl;
    fileParams << "  mu_max_hold = 2.773 (log-mean for max hold, median = 16h)" << endl;
    fileParams << "  w_E_hold = 10.0 (saturation width)" << endl;
    fileParams << "  sigmaT_hold = 0.5 (lognormal std dev)" << endl;
    fileParams << "\nLog-space Error (all doses):" << endl;
    fileParams << "  sf_eps = " << sf_eps << " (assay resolution floor)" << endl;
    fileParams << "  Average |Δlog10(SF)| = " << avg_log10_error << " (≈×" << avg_factor_error << ")" << endl;
    fileParams.close();
    
    cout << "\n--- Validation Complete ---" << endl;
    cout << "  Wall time: " << fixed << setprecision(2) << (wall_end - wall_start) << " seconds" << endl;
    cout << "  Output files:" << endl;
    cout << "    - " << outputDir << "/comparison.csv" << endl;
    cout << "    - " << outputDir << "/cellSurvival.csv" << endl;
    cout << "    - " << outputDir << "/checkpoint_diagnostics.csv" << endl;
    cout << "    - " << outputDir << "/calibrated_params.txt" << endl;
    
    return 0;
}
