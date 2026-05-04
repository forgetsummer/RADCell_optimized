/**
 * @file test_multicomponent_calibration.cc
 * @brief Joint calibration of all parameters including Nc, omega_p, lambda_p
 *
 * Replaces the grid-scan workflow by using the analytical multi-component
 * model to fit all 7 parameters jointly via Nelder-Mead.
 *
 * Workflow per PIDE cell line:
 * 1. Generate synthetic LQ survival data
 * 2. Call RepairMediatedMultiComponentCalibrator::fit() for all 7 parameters
 * 3. Run Monte Carlo simulation with calibrated parameters
 * 4. Compare analytical prediction, simulation, and LQ target
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
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

#include "RepairMediatedCalibrator_multicomponent.hh"
#include "DataTypes.hh"
#include "PIDEDataReader.hh"

using namespace std;
using namespace CellStateCalibration;
using namespace PIDE;

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
    double Nc, omega_p, lambda_p;
    bool calibrationConverged;
    int calibrationIterations;
    double calibrationNLL;
    vector<double> doses;
    vector<double> sf_lq;
    vector<double> sf_calibrated;
    vector<double> sf_sim;
    vector<double> sf_sim_sd;
    double avg_log10_error;
    double max_log10_error;
    double avg_log10_error_analytical;
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
    bool runSimulation)
{
    CellValidationResult result;
    result.cellLine = cellLine;
    result.photonRadiation = photonRadiation;
    result.alpha_LQ = alpha_LQ;
    result.beta_LQ = beta_LQ;
    result.alpha_beta_ratio = (beta_LQ > 0) ? (alpha_LQ / beta_LQ) : 0.0;

    cout << "\n  Validating: " << cellLine
         << " (alpha=" << alpha_LQ << ", beta=" << beta_LQ << ")" << endl;

    vector<double> fakeDoses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
    vector<DataPoint> fakeData;
    for (double dose : fakeDoses) {
        double sf = LQ_survival(dose, alpha_LQ, beta_LQ);
        fakeData.push_back(DataPoint(dose, sf, 0.0));
    }

    const double T_cellCycle = 10.0;
    const double lambda1_base = 3.31;
    const double lambda2_base = 0.14;

    // Compute effective repair time (same formula used in the simulation for To21/To23)
    const double lambda1_scaled = lambda1_base * repairRateScale;
    const double lambda2_scaled = lambda2_base * repairRateScale;
    const double f1_cfg = 0.62, f2_cfg = 0.38;
    const double effective_repair_time = f1_cfg / lambda1_scaled + f2_cfg / lambda2_scaled;

    FitConfigS2 cfg;
    cfg.kappa = 40.0;
    cfg.dose_max_fit = 6.0;
    cfg.T_assay_h = 27.8;
    // S2 timescales: fix to effective_repair_time (matching simulation's To21/To23)
    cfg.fixTimescales(effective_repair_time, effective_repair_time);
    cfg.f1 = f1_cfg;
    cfg.f2 = f2_cfg;
    cfg.lambda1_rep = lambda1_scaled;
    cfg.lambda2_rep = lambda2_scaled;
    cfg.fix_Nc = false;
    // S1 timescales: cell cycle time (matching simulation's To12/To13)
    cfg.To12 = T_cellCycle;
    cfg.To13 = T_cellCycle;

    RepairMediatedMultiComponentCalibrator calibrator(cfg);
    calibrator.setData(fakeData);
    calibrator.setDefaultUncertainties(0.1);

    ParamsS2 x0;
    x0.e2 = 2.0; x0.e3 = 10.0; x0.a = 0.05;
    x0.T21 = T_cellCycle; x0.T23 = T_cellCycle; x0.k_error = 1e-4;
    x0.Nc = 30.0; x0.omega_p = 1.5; x0.lambda_p = 0.03;

    ParamsS2 pstep;
    pstep.e2 = 0.5; pstep.e3 = 2.0; pstep.a = 0.01;
    pstep.T21 = 1.0; pstep.T23 = 1.0; pstep.k_error = 5e-5;
    pstep.Nc = 10.0; pstep.omega_p = 0.3; pstep.lambda_p = 0.01;

    MultiComponentFitResult fitResult = calibrator.fit(x0, pstep);

    result.e2 = fitResult.params.e2;
    result.e3 = fitResult.params.e3;
    result.a = fitResult.params.a;
    result.k_error = fitResult.params.k_error;
    result.T21 = fitResult.params.T21;
    result.T23 = fitResult.params.T23;
    result.Nc = fitResult.params.Nc;
    result.omega_p = fitResult.params.omega_p;
    result.lambda_p = fitResult.params.lambda_p;
    result.calibrationConverged = fitResult.converged;
    result.calibrationIterations = fitResult.iterations;
    result.calibrationNLL = fitResult.negLogLikelihood;

    double sigma_chosen = 10.0;
    CellStateModelParamsS2 calibratedParams = CellStateModelParamsS2::fromReduced(
        fitResult.params, sigma_chosen, cfg.kappa);

    cout << "    Fit: e2=" << fixed << setprecision(2) << result.e2
         << " e3=" << result.e3 << " a=" << result.a
         << " k_err=" << scientific << setprecision(2) << result.k_error
         << " Nc=" << fixed << setprecision(1) << result.Nc
         << " wp=" << setprecision(2) << result.omega_p
         << " lp=" << setprecision(3) << result.lambda_p
         << " (" << (fitResult.converged ? "converged" : "NOT converged")
         << ", " << fitResult.iterations << " iters)" << endl;

    // Compute analytical error
    vector<double> simDoses = fakeDoses;
    result.doses = simDoses;
    result.sf_lq.resize(simDoses.size());
    result.sf_calibrated.resize(simDoses.size());
    result.sf_sim.resize(simDoses.size(), 0.0);
    result.sf_sim_sd.resize(simDoses.size(), 0.0);

    vector<double> sf_calibrated_vec = calibrator.predictSurvival(fitResult.params, simDoses);
    const double sf_eps = std::max(1e-12, 0.5 / (static_cast<double>(cellNum) * static_cast<double>(numReplicates)));

    double total_analytical_error = 0.0;
    for (size_t i = 0; i < simDoses.size(); i++) {
        result.sf_lq[i] = LQ_survival(simDoses[i], alpha_LQ, beta_LQ);
        result.sf_calibrated[i] = sf_calibrated_vec[i];
        double log_err = std::fabs(std::log10(result.sf_calibrated[i] + sf_eps) -
                                   std::log10(result.sf_lq[i] + sf_eps));
        total_analytical_error += log_err;
    }
    result.avg_log10_error_analytical = total_analytical_error / simDoses.size();

    cout << "    Analytical avg |dlog10(SF)|: " << fixed << setprecision(3)
         << result.avg_log10_error_analytical << endl;

    // Always write per-cell-line CSV (analytical columns are always available)
    {
        string safeCellLine = sanitizeFilename(cellLine);
        string cellOutputFile = outputDir + "/mc_calib_" + safeCellLine + ".csv";
        ofstream cellFile(cellOutputFile);
        cellFile << "Dose,SF_LQ,SF_Calibrated,SF_Simulated,SF_Simulated_SD,AbsLog10Error" << endl;
        for (size_t i = 0; i < simDoses.size(); i++) {
            cellFile << simDoses[i] << "," << result.sf_lq[i] << "," << result.sf_calibrated[i]
                     << ",0,0,0" << endl;
        }
        cellFile.close();
    }

    if (!runSimulation) {
        result.avg_log10_error = result.avg_log10_error_analytical;
        result.max_log10_error = 0.0;
        return result;
    }

    // --- Run simulation with calibrated parameters ---
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
    double f1 = f1_cfg, f2 = f2_cfg;
    double lambda1 = lambda1_scaled;
    double lambda2 = lambda2_scaled;

    // Always use effective_repair_time for To21/To23 to match calibrator's assumption
    double T21_to_use = effective_repair_time;
    double T23_to_use = effective_repair_time;

    double DSB_per_Gy = cfg.kappa;

    CellLayoutInitializer layout;
    layout.SetRandomSeed(baseSeed);
    layout.SetCellHomeParamter(gridSize, gridSize, gridSize);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);

    #pragma omp parallel for schedule(dynamic)
    for (size_t n = 0; n < simDoses.size(); n++)
    {
        double dose = simDoses[n];
        double meanDSB = dose * DSB_per_Gy;

        std::vector<int> initialDSB(cellNum, static_cast<int>(meanDSB));

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
            myCellState.SetMisrepairRate("Epithelia", calibratedParams.k_error);
            myCellState.SetMultiComponentParams("Epithelia",
                calibratedParams.Nc, calibratedParams.omega_p, calibratedParams.lambda_p);
            myCellState.SetUpContactInhibition(true);
            myCellState.SetCheckpointEnabled(true);
            myCellState.SetSoftSaturationEnabled(true);

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

        result.sf_sim[n] = sfStats.mean();
        result.sf_sim_sd[n] = sfStats.sd();
    }

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

    cout << "    Sim avg |dlog10(SF)|: " << fixed << setprecision(3)
         << result.avg_log10_error << endl;

    // Per-cell-line CSV
    string safeCellLine = sanitizeFilename(cellLine);
    string cellOutputFile = outputDir + "/mc_calib_" + safeCellLine + ".csv";
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
    string outputDir = "./multicomponent_calibration_output";
    bool runSimulation = true;
    // (removed useFittedTimescales flag — always use effective_repair_time)
    vector<string> filterCellLines;

    for (int i = 1; i < argc; i++) {
        string arg = string(argv[i]);
        if (arg == "--pide-dir" && i + 1 < argc) { pideDir = argv[++i]; }
        else if (arg == "--output-dir" && i + 1 < argc) { outputDir = argv[++i]; }
        else if (arg == "--replicates" && i + 1 < argc) { numReplicates = atoi(argv[++i]); }
        else if (arg == "--cells" && i + 1 < argc) { cellNum = atoi(argv[++i]); }
        else if (arg == "--max-cell-types" && i + 1 < argc) { maxCellTypes = atoi(argv[++i]); }
        else if (arg == "--repair-scale" && i + 1 < argc) { repairRateScale = atof(argv[++i]); }
        else if (arg == "--analytical-only") { runSimulation = false; }
        // (removed --use-fitted-timescales / --use-effective-repair-time flags)
        else if (arg == "--cell-lines" && i + 1 < argc) {
            string list = argv[++i];
            stringstream ss(list);
            string token;
            while (getline(ss, token, ',')) {
                if (!token.empty()) filterCellLines.push_back(token);
            }
        }
        else if (arg == "--help" || arg == "-h") {
            cout << "Usage: " << argv[0] << " [OPTIONS]" << endl;
            cout << "Options:" << endl;
            cout << "  --pide-dir PATH      PIDE3.4 directory (default: ../PIDE3.4)" << endl;
            cout << "  --output-dir PATH    Output directory" << endl;
            cout << "  --replicates N       Replicates per dose (default: 3)" << endl;
            cout << "  --cells N            Cells per simulation (default: 200)" << endl;
            cout << "  --max-cell-types N   Maximum cell types (default: 10)" << endl;
            cout << "  --repair-scale V     Repair rate scale (default: 0.1)" << endl;
            cout << "  --analytical-only    Skip MC simulation, only test analytical fit" << endl;
            cout << "  --use-fitted-timescales   Use calibrated T21/T23 for simulation (default)" << endl;
            cout << "  --use-effective-repair-time  Use effective_repair_time for To21/To23" << endl;
            cout << "  --cell-lines A,B,C   Comma-separated list of cell lines to run" << endl;
            return 0;
        }
    }

    cout << "============================================" << endl;
    cout << "  Multi-Component Joint Calibration" << endl;
    cout << "  (Nc, omega_p, lambda_p fitted analytically)" << endl;
    cout << "============================================" << endl;
    cout << "Configuration:" << endl;
    cout << "  PIDE directory:   " << pideDir << endl;
    cout << "  Output directory:  " << outputDir << endl;
    cout << "  Replicates/dose:  " << numReplicates << endl;
    cout << "  Cells/sim:        " << cellNum << endl;
    cout << "  Repair scale:     " << repairRateScale << endl;
    cout << "  Run simulation:   " << (runSimulation ? "yes" : "no (analytical only)") << endl;
    cout << "  Timescale mode:   effective_repair_time (matching simulation To21/To23)" << endl;
    if (!filterCellLines.empty()) {
        cout << "  Cell line filter: ";
        for (size_t i = 0; i < filterCellLines.size(); i++) {
            if (i > 0) cout << ", ";
            cout << filterCellLines[i];
        }
        cout << endl;
    }

    createDirectoryIfNotExists(outputDir.c_str());

    cout << "\n--- Loading PIDE Database ---" << endl;
    PIDEDataReader reader;
    if (!reader.loadData(pideDir)) {
        cerr << "ERROR: Failed to load PIDE data from: " << pideDir << endl;
        return 1;
    }
    reader.printSummary();

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

    if (!filterCellLines.empty()) {
        vector<ExperimentEntry> filtered;
        for (const auto& exp : co60Exps) {
            for (const auto& target : filterCellLines) {
                if (exp.cellLine == target) {
                    filtered.push_back(exp);
                    break;
                }
            }
        }
        co60Exps = filtered;
        cout << "Filtered to " << co60Exps.size() << " cell lines" << endl;
    } else if (co60Exps.size() > static_cast<size_t>(maxCellTypes)) {
        co60Exps.resize(maxCellTypes);
    }

    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    cout << "\n--- Running Calibration (OpenMP threads: " << numThreads << ") ---" << endl;

    double wall_start = omp_get_wtime();

    vector<CellValidationResult> allResults;
    unsigned long long baseSeed = 12345ULL;

    for (size_t i = 0; i < co60Exps.size(); i++) {
        const auto& exp = co60Exps[i];
        const LQParameters& lq = exp.photonLQ_paper.isValid ? exp.photonLQ_paper : exp.photonLQ_fit;

        CellValidationResult result = runValidationForCellType(
            exp.cellLine, exp.photonRadiation,
            lq.alpha, lq.beta,
            cellNum, numReplicates, repairRateScale,
            outputDir, baseSeed + i * 1000000ULL,
            runSimulation);

        allResults.push_back(result);
    }

    double wall_end = omp_get_wtime();

    // Save summary CSV
    string summaryFile = outputDir + "/mc_calibration_summary.csv";
    ofstream summaryCSV(summaryFile);
    summaryCSV << "CellLine,PhotonRadiation,Alpha_LQ,Beta_LQ,AlphaBetaRatio,"
               << "e2,e3,a,k_error,T21,T23,Nc,omega_p,lambda_p,"
               << "CalibConverged,CalibIters,CalibNLL,"
               << "AvgLog10Error_Analytical,AvgLog10Error_Sim,MaxLog10Error_Sim" << endl;

    for (const auto& r : allResults) {
        summaryCSV << r.cellLine << "," << r.photonRadiation << ","
                   << r.alpha_LQ << "," << r.beta_LQ << "," << r.alpha_beta_ratio << ","
                   << r.e2 << "," << r.e3 << "," << r.a << "," << r.k_error << ","
                   << r.T21 << "," << r.T23 << ","
                   << r.Nc << "," << r.omega_p << "," << r.lambda_p << ","
                   << (r.calibrationConverged ? "yes" : "no") << ","
                   << r.calibrationIterations << "," << r.calibrationNLL << ","
                   << r.avg_log10_error_analytical << ","
                   << r.avg_log10_error << "," << r.max_log10_error << endl;
    }
    summaryCSV.close();

    // Print final summary
    cout << "\n============================================" << endl;
    cout << "  JOINT CALIBRATION COMPLETE" << endl;
    cout << "============================================" << endl;
    cout << "Wall time: " << fixed << setprecision(1) << (wall_end - wall_start) << " s" << endl;

    cout << "\n" << setw(15) << "Cell Line" << setw(8) << "alpha"
         << setw(8) << "beta" << setw(7) << "Nc"
         << setw(7) << "wp" << setw(8) << "lp"
         << setw(12) << "Anal.Err" << setw(12) << "Sim.Err" << endl;
    cout << string(77, '-') << endl;

    double totalAnalErr = 0.0;
    double totalSimErr = 0.0;
    for (const auto& r : allResults) {
        cout << setw(15) << r.cellLine
             << setw(8) << fixed << setprecision(3) << r.alpha_LQ
             << setw(8) << r.beta_LQ
             << setw(7) << setprecision(1) << r.Nc
             << setw(7) << setprecision(2) << r.omega_p
             << setw(8) << setprecision(3) << r.lambda_p
             << setw(12) << setprecision(3) << r.avg_log10_error_analytical
             << setw(12) << r.avg_log10_error << endl;
        totalAnalErr += r.avg_log10_error_analytical;
        totalSimErr += r.avg_log10_error;
    }

    cout << string(77, '-') << endl;
    double overallAnalAvg = totalAnalErr / allResults.size();
    double overallSimAvg = totalSimErr / allResults.size();
    cout << "Overall avg |dlog10(SF)| analytical: " << fixed << setprecision(3) << overallAnalAvg << endl;
    if (runSimulation) {
        cout << "Overall avg |dlog10(SF)| simulation: " << overallSimAvg << endl;
    }

    cout << "\nOutput: " << summaryFile << endl;
    return 0;
}
