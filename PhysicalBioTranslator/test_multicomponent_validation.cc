/**
 * @file test_multicomponent_validation.cc
 * @brief Multi-component energy model (Er + Ep) validation against PIDE 60Co data.
 *
 * Workflow:
 * 1. Load PIDE 60Co photon survival data
 * 2. Analytically calibrate 6 parameters per cell line (unchanged calibrator)
 * 3. Run simulation with Er/Ep split controlled by (Nc, omega_p, lambda_p)
 * 4. Output per-cell-line SF comparison CSV and aggregate error metrics
 *
 * New CLI flags: --nc, --omega-p, --lambda-p
 * These are passed to CellStateModel::SetMultiComponentParams().
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

#include "RepairMediatedCalibrator_with_k_error.hh"
#include "RepairMediatedCalibrator.hh"
#include "DataTypes.hh"
#include "PIDEDataReader.hh"

using namespace std;
using namespace CellStateCalibration;
using namespace PIDE;

struct PerCellLineParams {
    double Nc;
    double omega_p;
    double lambda_p;
};

map<string, PerCellLineParams> loadParamsFile(const string& path) {
    map<string, PerCellLineParams> result;
    ifstream f(path);
    if (!f.is_open()) {
        cerr << "ERROR: Cannot open params file: " << path << endl;
        return result;
    }
    string line;
    getline(f, line); // skip header
    while (getline(f, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        string cellLine;
        char comma;
        PerCellLineParams p;
        getline(iss, cellLine, ',');
        iss >> p.Nc >> comma >> p.omega_p >> comma >> p.lambda_p;
        if (!iss.fail()) {
            result[cellLine] = p;
            cout << "  Loaded params for " << cellLine
                 << ": Nc=" << p.Nc << " wp=" << p.omega_p << " lp=" << p.lambda_p << endl;
        }
    }
    return result;
}

double LQ_survival(double dose, double alpha_LQ, double beta_LQ) {
    return exp(-alpha_LQ * dose - beta_LQ * dose * dose);
}

struct CellValidationResult {
    string cellLine;
    string photonRadiation;
    double alpha_LQ;
    double beta_LQ;
    double alpha_beta_ratio;
    double e2, e3, a, k_error;
    double T21, T23;
    bool calibrationConverged;
    vector<double> doses;
    vector<double> sf_lq;
    vector<double> sf_calibrated;
    vector<double> sf_sim;
    vector<double> sf_sim_sd;
    double avg_log10_error;
    double max_log10_error;
};

void createDirectoryIfNotExists(const char* path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0755);
    }
}

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

CellValidationResult runValidationForCellType(
    const string& cellLine,
    const string& photonRadiation,
    double alpha_LQ,
    double beta_LQ,
    int cellNum,
    int numReplicates,
    double repairRateScale,
    const string& outputDir,
    unsigned long long baseSeed,
    double param_Nc,
    double param_omega_p,
    double param_lambda_p,
    bool poissonDSBEnabled,
    bool useFittedTimescales,
    int checkpointMode,
    double lambdaCommit,
    bool noCheckpoint,
    bool noMisrepair,
    bool force6paramCalib)
{
    CellValidationResult result;
    result.cellLine = cellLine;
    result.photonRadiation = photonRadiation;
    result.alpha_LQ = alpha_LQ;
    result.beta_LQ = beta_LQ;
    result.alpha_beta_ratio = (beta_LQ > 0) ? (alpha_LQ / beta_LQ) : 0.0;

    cout << "\n  Validating: " << cellLine
         << " (α=" << alpha_LQ << ", β=" << beta_LQ << ")" << endl;

    // Generate synthetic LQ data for calibration
    vector<double> fakeDoses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
    vector<DataPoint> fakeData;
    for (double dose : fakeDoses) {
        double sf = LQ_survival(dose, alpha_LQ, beta_LQ);
        fakeData.push_back(DataPoint(dose, sf, 0.0));
    }

    const double T_cellCycle = 10.0;

    FitConfigS2 cfg;
    cfg.kappa = 40.0;
    cfg.dose_max_fit = 6.0;
    cfg.T_assay_h = 27.8;
    cfg.setTimescalesToCellCycle(T_cellCycle);

    ParamsS2 fitParams;
    bool fitConverged = false;
    vector<double> sf_calib_predict;

    bool use5paramCalib = noMisrepair && !force6paramCalib;

    if (use5paramCalib) {
        // 5-param calibration (no k_error)
        RepairMediatedCalibrator calibrator5(cfg);
        calibrator5.setData(fakeData);
        calibrator5.setDefaultUncertainties(0.1);

        ParamsS2 x0;
        x0.e2 = 2.0; x0.e3 = 10.0; x0.a = 0.05;
        x0.T21 = T_cellCycle; x0.T23 = T_cellCycle; x0.k_error = 0.0;

        ParamsS2 step;
        step.e2 = 0.5; step.e3 = 2.0; step.a = 0.01;
        step.T21 = 1.0; step.T23 = 1.0; step.k_error = 0.0;

        RepairMediatedFitResult fitResult5 = calibrator5.fit(x0, step);
        fitParams = fitResult5.params;
        fitParams.k_error = 0.0;
        fitConverged = fitResult5.converged;

        vector<double> simDosesTmp = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
        sf_calib_predict = calibrator5.predictSurvival(fitResult5.params, simDosesTmp);
    } else {
        // 6-param calibration (with k_error)
        RepairMediatedMisrepairCalibrator calibrator6(cfg);
        calibrator6.setData(fakeData);
        calibrator6.setDefaultUncertainties(0.1);

        ParamsS2 x0;
        x0.e2 = 2.0; x0.e3 = 10.0; x0.a = 0.05;
        x0.T21 = T_cellCycle; x0.T23 = T_cellCycle; x0.k_error = 1e-4;

        ParamsS2 step;
        step.e2 = 0.5; step.e3 = 2.0; step.a = 0.01;
        step.T21 = 1.0; step.T23 = 1.0; step.k_error = 5e-5;

        RepairMediatedMisrepairFitResult fitResult6 = calibrator6.fit(x0, step);
        fitParams = fitResult6.params;
        fitConverged = fitResult6.converged;

        vector<double> simDosesTmp = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
        sf_calib_predict = calibrator6.predictSurvival(fitResult6.params, simDosesTmp);
    }

    result.e2 = fitParams.e2;
    result.e3 = fitParams.e3;
    result.a = fitParams.a;
    result.k_error = fitParams.k_error;
    result.T21 = fitParams.T21;
    result.T23 = fitParams.T23;
    result.calibrationConverged = fitConverged;

    double sigma_chosen = 10.0;
    CellStateModelParamsS2 calibratedParams = CellStateModelParamsS2::fromReduced(
        fitParams, sigma_chosen, cfg.kappa);
    calibratedParams.k_error = fitParams.k_error;

    // Simulation parameters
    double xDim = 5.0, yDim = 5.0, zDim = 0.0;
    double gridSize = 0.05;
    double deltaT_cellPhaseUpdate = 60;
    double deltaT_cellStateUpdate = 60;
    double T = 10;
    int totalTimeStepNum = 10000;

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
    double lambda1 = 3.31 * repairRateScale;
    double lambda2 = 0.14 * repairRateScale;

    double tau1 = 1.0 / lambda1;
    double tau2 = 1.0 / lambda2;
    double effective_repair_time = f1 * tau1 + f2 * tau2;
    double T21_to_use = useFittedTimescales ? fitParams.T21 : effective_repair_time;
    double T23_to_use = useFittedTimescales ? fitParams.T23 : effective_repair_time;

    double DSB_per_Gy = cfg.kappa;

    vector<double> simDoses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};

    CellLayoutInitializer layout;
    layout.SetRandomSeed(baseSeed);
    layout.SetCellHomeParamter(gridSize, gridSize, gridSize);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);

    result.doses = simDoses;
    result.sf_lq.resize(simDoses.size());
    result.sf_calibrated.resize(simDoses.size());
    result.sf_sim.resize(simDoses.size());
    result.sf_sim_sd.resize(simDoses.size());

    for (size_t i = 0; i < simDoses.size(); i++) {
        result.sf_calibrated[i] = sf_calib_predict[i];
    }

    const double sf_eps = std::max(1e-12, 0.5 / (static_cast<double>(cellNum) * static_cast<double>(numReplicates)));

    #pragma omp parallel for schedule(dynamic)
    for (size_t n = 0; n < simDoses.size(); n++)
    {
        double dose = simDoses[n];
        double meanDSB = dose * DSB_per_Gy;

        std::vector<int> initialDSB(cellNum, 0);
        if (poissonDSBEnabled && meanDSB > 0) {
            std::mt19937_64 dsbRng(baseSeed + n * 999999ULL);
            std::poisson_distribution<int> poissonDist(meanDSB);
            for (int cid = 0; cid < cellNum; cid++) {
                initialDSB[cid] = poissonDist(dsbRng);
            }
        } else {
            for (int cid = 0; cid < cellNum; cid++) {
                initialDSB[cid] = static_cast<int>(meanDSB);
            }
        }

        RunningStats sfStats;

        for (int rep = 0; rep < numReplicates; rep++)
        {
            double clock_cellPhaseUpdate = 0;
            double clock_cellStateUpdate = 0;
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
            statePara.E1 = E1; statePara.E2 = E2; statePara.E3 = E3;
            statePara.sigma = sigma_model; statePara.alpha = alpha_model;
            statePara.beta = beta_model;
            statePara.To12 = T_cellCycle; statePara.To13 = T_cellCycle;
            statePara.To21 = T21_to_use; statePara.To23 = T23_to_use;

            CellDNADamageRepairParameter dnaRepairPara;
            dnaRepairPara.f1 = f1; dnaRepairPara.f2 = f2;
            dnaRepairPara.lambda1 = lambda1; dnaRepairPara.lambda2 = lambda2;

            CellBystanderSignalParameter bystanderPara;
            bystanderPara.Rt = 5000; bystanderPara.mu = 200;

            testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);

            CellStateModel myCellState;
            myCellState.SetRandomSeed(baseSeed + n * 100000ULL + rep * 1000ULL);
            myCellState.TissueGeometryInitialization(xDim, yDim, zDim, gridSize);
            myCellState.CellStateModelParameterSetup(testCell);
            if (!noMisrepair) {
                myCellState.SetMisrepairRate("Epithelia", calibratedParams.k_error);
            }
            myCellState.SetMultiComponentParams("Epithelia", param_Nc, param_omega_p, param_lambda_p);
            myCellState.SetUpContactInhibition(true);
            myCellState.SetCheckpointEnabled(!noCheckpoint);
            myCellState.SetSoftSaturationEnabled(true);
            myCellState.SetCheckpointMode(checkpointMode);
            myCellState.SetCommitmentHazardRate(lambdaCommit);

            for (int i = 0; i < cellNum; i++) {
                double cX = layout.GetCellPositionX(i) + xDim / 2.0;
                double cY = layout.GetCellPositionY(i) + yDim / 2.0;
                double cZ = layout.GetCellPositionZ(i) + zDim / 2.0;
                myCellState.CellPositionInitialization(i, cX, cY, cZ);
                myCellState.CellTypeInitialiation(i, testCell);
                myCellState.CellPhaseInitializationRandom(i);
                myCellState.CellStateInitialization(i, "S1");
            }

            for (int i = 0; i < totalTimeStepNum; i++)
            {
                clock_cellStateUpdate += T;
                clock_cellPhaseUpdate += T;

                map<int, string> phaseMap = myCellState.GetCellPhase();

                if (clock_cellPhaseUpdate >= deltaT_cellPhaseUpdate) {
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

            for (const auto& cell : finalStateMap) {
                if (cell.second != "S1") continue;
                auto itA = cellAncestryIDMap.find(cell.first);
                if (itA == cellAncestryIDMap.end()) continue;
                cellColonySizeMap[itA->second]++;
            }

            int colonyNumber = 0;
            for (const auto& colony : cellColonySizeMap) {
                if (colony.second >= 4) colonyNumber++;
            }

            sfStats.add(static_cast<double>(colonyNumber) / static_cast<double>(cellNum));
        }

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
        if (log10_error > result.max_log10_error)
            result.max_log10_error = log10_error;
    }
    result.avg_log10_error = total_log10_error / simDoses.size();

    cout << "    Avg |Δlog10(SF)|: " << fixed << setprecision(3) << result.avg_log10_error
         << "  (e2=" << setprecision(2) << result.e2
         << " e3=" << result.e3
         << " a=" << result.a << ")" << endl;

    // Save per-cell-line CSV
    string safeCellLine = sanitizeFilename(cellLine);
    string cellOutputFile = outputDir + "/mc_validation_" + safeCellLine + ".csv";
    ofstream cellFile(cellOutputFile);
    cellFile << "Dose,SF_LQ,SF_Calibrated,SF_Simulated,SF_Simulated_SD,AbsLog10Error" << endl;
    for (size_t n = 0; n < simDoses.size(); n++) {
        double log10_error = std::fabs(std::log10(result.sf_sim[n] + sf_eps) -
                                       std::log10(result.sf_lq[n] + sf_eps));
        cellFile << simDoses[n] << "," << result.sf_lq[n] << "," << result.sf_calibrated[n] << ","
                 << result.sf_sim[n] << "," << result.sf_sim_sd[n] << "," << log10_error << endl;
    }
    cellFile.close();

    return result;
}

int main(int argc, char** argv)
{
    int numReplicates = 3;
    int cellNum = 200;
    int maxCellTypes = 10;
    double repairRateScale = 0.1;
    string pideDir = "../PIDE3.4";
    string outputDir = "./multicomponent_validation_output";

    double param_Nc = 30.0;
    double param_omega_p = 1.5;
    double param_lambda_p = 0.03;
    bool poissonDSBEnabled = false;
    bool useFittedTimescales = false;
    int checkpointMode = 0;
    double lambdaCommit = 0.0;
    string paramsFile = "";
    string cellLineFilter = "";
    bool noCheckpoint = false;
    bool noMisrepair = false;
    bool force6paramCalib = false;

    for (int i = 1; i < argc; i++) {
        string arg = string(argv[i]);
        if (arg == "--pide-dir" && i + 1 < argc) { pideDir = argv[++i]; }
        else if (arg == "--output-dir" && i + 1 < argc) { outputDir = argv[++i]; }
        else if (arg == "--replicates" && i + 1 < argc) { numReplicates = atoi(argv[++i]); }
        else if (arg == "--cells" && i + 1 < argc) { cellNum = atoi(argv[++i]); }
        else if (arg == "--max-cell-types" && i + 1 < argc) { maxCellTypes = atoi(argv[++i]); }
        else if (arg == "--repair-scale" && i + 1 < argc) { repairRateScale = atof(argv[++i]); }
        else if (arg == "--nc" && i + 1 < argc) { param_Nc = atof(argv[++i]); }
        else if (arg == "--omega-p" && i + 1 < argc) { param_omega_p = atof(argv[++i]); }
        else if (arg == "--lambda-p" && i + 1 < argc) { param_lambda_p = atof(argv[++i]); }
        else if (arg == "--params-file" && i + 1 < argc) { paramsFile = argv[++i]; }
        else if (arg == "--cell-line" && i + 1 < argc) { cellLineFilter = argv[++i]; }
        else if (arg == "--poisson-dsb") { poissonDSBEnabled = true; }
        else if (arg == "--use-fitted-timescales") { useFittedTimescales = true; }
        else if (arg == "--no-checkpoint") { noCheckpoint = true; }
        else if (arg == "--no-misrepair") { noMisrepair = true; }
        else if (arg == "--force-6param-calib") { force6paramCalib = true; }
        else if (arg == "--checkpoint-mode" && i + 1 < argc) { checkpointMode = atoi(argv[++i]); }
        else if (arg == "--lambda-commit" && i + 1 < argc) { lambdaCommit = atof(argv[++i]); }
        else if (arg == "--help" || arg == "-h") {
            cout << "Usage: " << argv[0] << " [OPTIONS]" << endl;
            cout << "Options:" << endl;
            cout << "  --pide-dir PATH      PIDE3.4 directory (default: ../PIDE3.4)" << endl;
            cout << "  --output-dir PATH    Output directory" << endl;
            cout << "  --replicates N       Replicates per dose (default: 3)" << endl;
            cout << "  --cells N            Cells per simulation (default: 200)" << endl;
            cout << "  --max-cell-types N   Maximum cell types (default: 10)" << endl;
            cout << "  --repair-scale V     Repair rate scale (default: 0.1)" << endl;
            cout << "  --nc V               Half-saturation DSB count Nc (default: 30)" << endl;
            cout << "  --omega-p V          Persistent energy weight omega_p (default: 1.5)" << endl;
            cout << "  --lambda-p V         Persistent decay rate 1/hr (default: 0.03)" << endl;
            cout << "  --params-file PATH   CSV with per-cell-line Nc,omega_p,lambda_p" << endl;
            cout << "  --cell-line NAME     Run only this cell line (e.g. CHO-10B)" << endl;
            cout << "  --no-checkpoint      Disable damage checkpoint (diversion + gating)" << endl;
            cout << "  --no-misrepair       Disable misrepair channel (k_error=0)" << endl;
            cout << "  --force-6param-calib Force 6-param calibrator even with --no-misrepair" << endl;
            cout << "  --poisson-dsb        Enable Poisson DSB sampling per cell" << endl;
            cout << "  --use-fitted-timescales  Use calibrated T21/T23 (default: use effective repair time)" << endl;
            cout << "  --checkpoint-mode N  Checkpoint mode (0=gating, 1=attempt-based, default: 0)" << endl;
            cout << "  --lambda-commit V    Commitment hazard rate (default: 0.0)" << endl;
            return 0;
        }
    }

    map<string, PerCellLineParams> perCellParams;
    if (!paramsFile.empty()) {
        cout << "Loading per-cell-line parameters from: " << paramsFile << endl;
        perCellParams = loadParamsFile(paramsFile);
        if (perCellParams.empty()) {
            cerr << "WARNING: No parameters loaded from params file!" << endl;
        }
    }

    cout << "============================================" << endl;
    cout << "  Multi-Component Energy Model Validation" << endl;
    cout << "  Er + Ep split with PIDE 60Co data" << endl;
    cout << "============================================" << endl;
    cout << "Configuration:" << endl;
    cout << "  PIDE directory:   " << pideDir << endl;
    cout << "  Output directory:  " << outputDir << endl;
    cout << "  Replicates/dose:  " << numReplicates << endl;
    cout << "  Cells/sim:        " << cellNum << endl;
    cout << "  Repair scale:     " << repairRateScale << endl;
    cout << "  Nc:               " << param_Nc << endl;
    cout << "  omega_p:          " << param_omega_p << endl;
    cout << "  lambda_p:         " << param_lambda_p << " (half-life "
         << fixed << setprecision(1) << (param_lambda_p > 0 ? 0.693 / param_lambda_p : 1e9) << " h)" << endl;
    cout << "  Poisson DSB:      " << (poissonDSBEnabled ? "enabled" : "disabled") << endl;
    cout << "  Fitted timescales:" << (useFittedTimescales ? "enabled" : "disabled (use effective repair time)") << endl;
    cout << "  Checkpoint:       " << (noCheckpoint ? "DISABLED" : "enabled") << endl;
    cout << "  Misrepair:        " << (noMisrepair ? "DISABLED (k_error=0)" : "enabled") << endl;
    bool use5paramCalibMain = noMisrepair && !force6paramCalib;
    cout << "  Calibrator:       " << (use5paramCalibMain ? "RepairMediatedCalibrator (5-param, no k_error)" : "RepairMediatedMisrepairCalibrator (6-param)") << endl;
    cout << "  Checkpoint mode:  " << checkpointMode << endl;
    cout << "  Lambda commit:    " << lambdaCommit << endl;

    createDirectoryIfNotExists(outputDir.c_str());

    // Load PIDE data
    cout << "\n--- Loading PIDE Database ---" << endl;
    PIDEDataReader reader;
    if (!reader.loadData(pideDir)) {
        cerr << "ERROR: Failed to load PIDE data from: " << pideDir << endl;
        return 1;
    }
    reader.printSummary();

    // Find 60Co experiments
    vector<ExperimentEntry> co60Exps;
    map<string, ExperimentEntry> uniqueCellLines;

    for (size_t i = 0; i < reader.getNumExperiments(); i++) {
        const ExperimentEntry* exp = reader.getExperimentByID(i + 1);
        if (!exp) continue;

        string photonRad = exp->photonRadiation;
        bool is60Co = (photonRad.find("60Co") != string::npos ||
                       photonRad.find("Co60") != string::npos ||
                       photonRad.find("Co-60") != string::npos);
        if (!is60Co) continue;

        bool hasValidLQ = exp->photonLQ_paper.isValid || exp->photonLQ_fit.isValid;
        if (!hasValidLQ) continue;

        if (uniqueCellLines.find(exp->cellLine) == uniqueCellLines.end()) {
            uniqueCellLines[exp->cellLine] = *exp;
            co60Exps.push_back(*exp);
        }
    }

    cout << "Found " << co60Exps.size() << " unique 60Co cell lines" << endl;
    if (co60Exps.empty()) {
        cerr << "ERROR: No 60Co experiments found!" << endl;
        return 1;
    }

    if (!cellLineFilter.empty()) {
        vector<ExperimentEntry> filtered;
        for (const auto& e : co60Exps) {
            if (e.cellLine == cellLineFilter) {
                filtered.push_back(e);
                break;
            }
        }
        if (filtered.empty()) {
            cerr << "ERROR: Cell line '" << cellLineFilter << "' not found in PIDE 60Co data!" << endl;
            cerr << "Available cell lines:" << endl;
            for (const auto& e : co60Exps) cerr << "  " << e.cellLine << endl;
            return 1;
        }
        co60Exps = filtered;
        cout << "Filtered to cell line: " << cellLineFilter << endl;
    }

    if (co60Exps.size() > static_cast<size_t>(maxCellTypes)) {
        co60Exps.resize(maxCellTypes);
    }

    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    cout << "\n--- Running Validation (OpenMP threads: " << numThreads << ") ---" << endl;

    double wall_start = omp_get_wtime();

    vector<CellValidationResult> allResults;
    unsigned long long baseSeed = 12345ULL;

    for (size_t i = 0; i < co60Exps.size(); i++) {
        const auto& exp = co60Exps[i];
        const LQParameters& lq = exp.photonLQ_paper.isValid ? exp.photonLQ_paper : exp.photonLQ_fit;

        double use_Nc = param_Nc;
        double use_omega_p = param_omega_p;
        double use_lambda_p = param_lambda_p;

        if (!perCellParams.empty()) {
            auto it = perCellParams.find(exp.cellLine);
            if (it != perCellParams.end()) {
                use_Nc = it->second.Nc;
                use_omega_p = it->second.omega_p;
                use_lambda_p = it->second.lambda_p;
                cout << "  Using per-cell-line params for " << exp.cellLine
                     << ": Nc=" << use_Nc << " wp=" << use_omega_p << " lp=" << use_lambda_p << endl;
            } else {
                cout << "  No per-cell-line params for " << exp.cellLine
                     << ", using global defaults" << endl;
            }
        }

        CellValidationResult result = runValidationForCellType(
            exp.cellLine, exp.photonRadiation,
            lq.alpha, lq.beta,
            cellNum, numReplicates, repairRateScale,
            outputDir, baseSeed + i * 1000000ULL,
            use_Nc, use_omega_p, use_lambda_p,
            poissonDSBEnabled, useFittedTimescales,
            checkpointMode, lambdaCommit,
            noCheckpoint, noMisrepair, force6paramCalib);

        allResults.push_back(result);
    }

    double wall_end = omp_get_wtime();

    // Save summary CSV
    string summaryFile = outputDir + "/mc_validation_summary.csv";
    ofstream summaryCSV(summaryFile);
    summaryCSV << "CellLine,PhotonRadiation,Alpha_LQ,Beta_LQ,AlphaBetaRatio,"
               << "e2,e3,a,k_error,T21,T23,CalibConverged,"
               << "AvgLog10Error,MaxLog10Error,"
               << "Nc,omega_p,lambda_p" << endl;

    for (const auto& r : allResults) {
        double out_Nc = param_Nc, out_wp = param_omega_p, out_lp = param_lambda_p;
        if (!perCellParams.empty()) {
            auto it = perCellParams.find(r.cellLine);
            if (it != perCellParams.end()) {
                out_Nc = it->second.Nc;
                out_wp = it->second.omega_p;
                out_lp = it->second.lambda_p;
            }
        }
        summaryCSV << r.cellLine << "," << r.photonRadiation << ","
                   << r.alpha_LQ << "," << r.beta_LQ << "," << r.alpha_beta_ratio << ","
                   << r.e2 << "," << r.e3 << "," << r.a << "," << r.k_error << ","
                   << r.T21 << "," << r.T23 << "," << (r.calibrationConverged ? "yes" : "no") << ","
                   << r.avg_log10_error << "," << r.max_log10_error << ","
                   << out_Nc << "," << out_wp << "," << out_lp << endl;
    }
    summaryCSV.close();

    // Print final summary
    cout << "\n============================================" << endl;
    cout << "  VALIDATION COMPLETE" << endl;
    cout << "============================================" << endl;
    cout << "Nc=" << param_Nc << "  omega_p=" << param_omega_p
         << "  lambda_p=" << param_lambda_p << endl;
    cout << "Wall time: " << fixed << setprecision(1) << (wall_end - wall_start) << " s" << endl;

    cout << "\n" << setw(15) << "Cell Line" << setw(10) << "α"
         << setw(10) << "β" << setw(14) << "Avg|Δlog10|" << setw(14) << "Max|Δlog10|" << endl;
    cout << string(63, '-') << endl;

    double totalAvgError = 0.0;
    for (const auto& r : allResults) {
        cout << setw(15) << r.cellLine
             << setw(10) << fixed << setprecision(4) << r.alpha_LQ
             << setw(10) << r.beta_LQ
             << setw(14) << setprecision(3) << r.avg_log10_error
             << setw(14) << r.max_log10_error << endl;
        totalAvgError += r.avg_log10_error;
    }

    cout << string(63, '-') << endl;
    double overallAvg = totalAvgError / allResults.size();
    cout << "Overall avg |Δlog10(SF)|: " << fixed << setprecision(3) << overallAvg << endl;

    cout << "\nOutput: " << summaryFile << endl;
    return 0;
}
