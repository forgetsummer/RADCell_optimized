/**
 * @file test_simulation_calibration.cc
 * @brief Simulation-based calibration: optimizes parameters directly against
 *        the stochastic simulation output instead of an analytical model.
 *
 * Flow:
 * 1. For each PIDE cell line, run analytical calibration (warm start)
 * 2. Build a lambda that runs a coarse simulation and returns SF
 * 3. Pass it to SimulationBasedCalibrator to optimize all 6 parameters
 * 4. Validate with a detailed simulation and compare to analytical calibration
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

#include "RepairMediatedCalibrator_with_k_error.hh"
#include "SimulationBasedCalibrator.hh"
#include "DataTypes.hh"
#include "PIDEDataReader.hh"

using namespace std;
using namespace CellStateCalibration;
using namespace PIDE;

//============================================================================
// Helpers
//============================================================================
double LQ_survival(double dose, double alpha_LQ, double beta_LQ) {
    return exp(-alpha_LQ * dose - beta_LQ * dose * dose);
}

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
    int n = 0; double sum = 0.0, sumsq = 0.0;
    void add(double x) { n++; sum += x; sumsq += x * x; }
    double mean() const { return (n > 0) ? (sum / n) : 0.0; }
    double sd() const {
        if (n < 2) return 0.0;
        double m = mean();
        return std::sqrt(std::max(0.0, sumsq / n - m * m));
    }
};

//============================================================================
// Core simulation function: runs one simulation at a given dose and returns SF
//============================================================================
struct SimConfig {
    double sigma_chosen = 10.0;
    double kappa = 40.0;
    double repairRateScale = 0.1;
    double T_cellCycle = 10.0;
    double xDim = 5.0, yDim = 5.0, zDim = 0.0, d = 0.05;
    int totalTimeStepNum = 10000;
    double T_microstep = 10.0;
    double deltaT_update = 60.0;
};

double simulateSF(double dose, const ParamsS2& params,
                  int cellNum, int nReps,
                  unsigned long long baseSeed,
                  const SimConfig& sc)
{
    CellStateModelParamsS2 phys = CellStateModelParamsS2::fromReduced(
        params, sc.sigma_chosen, sc.kappa);
    phys.k_error = params.k_error;

    double f1 = 0.62, f2 = 0.38;
    double lambda1 = 3.31 * sc.repairRateScale;
    double lambda2 = 0.14 * sc.repairRateScale;

    double tG1 = 1.5, sigmaG1 = 0.25;
    double tS = 6.0,  sigmaS  = 0.25;
    double tG2 = 1.5, sigmaG2 = 0.25;
    double tM = 1.0,  sigmaM  = 0.25;
    double fG1 = 1.0, fS = 0.816326, fG2 = 1.015306122, fM = 1.015306122;
    double beta_model = 25.71;

    CellLayoutInitializer layout;
    layout.SetRandomSeed(baseSeed);
    layout.SetCellHomeParamter(sc.d, sc.d, sc.d);
    layout.RectangularSlab(sc.xDim, sc.yDim, sc.zDim, cellNum);

    double meanDSB = dose * sc.kappa;
    RunningStats sfStats;

    for (int rep = 0; rep < nReps; rep++) {
        Cell testCell;
        testCell.CellConstruct("Epithelia", "Cytoplasm Nucleus", "Sphere", "Blue Green");

        CellCycleParameter cyclePara;
        cyclePara.mTG1 = tG1; cyclePara.sigmaTG1 = sigmaG1;
        cyclePara.mTS = tS;   cyclePara.sigmaTS  = sigmaS;
        cyclePara.mTG2 = tG2; cyclePara.sigmaTG2 = sigmaG2;
        cyclePara.mTM = tM;   cyclePara.sigmaTM  = sigmaM;
        cyclePara.fG1 = fG1; cyclePara.fS = fS; cyclePara.fG2 = fG2; cyclePara.fM = fM;

        CellStateParameter statePara;
        statePara.E1 = phys.E1;
        statePara.E2 = phys.E2;
        statePara.E3 = phys.E3;
        statePara.sigma = phys.sigma;
        statePara.alpha = phys.alpha;
        statePara.beta = beta_model;
        statePara.To12 = sc.T_cellCycle;
        statePara.To13 = sc.T_cellCycle;
        statePara.To21 = params.T21;
        statePara.To23 = params.T23;

        CellDNADamageRepairParameter dnaRepairPara;
        dnaRepairPara.f1 = f1; dnaRepairPara.f2 = f2;
        dnaRepairPara.lambda1 = lambda1; dnaRepairPara.lambda2 = lambda2;

        CellBystanderSignalParameter bystanderPara;
        bystanderPara.Rt = 5000; bystanderPara.mu = 200;

        testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);

        CellStateModel myCellState;
        myCellState.SetRandomSeed(baseSeed + static_cast<unsigned long long>(rep) * 7777ULL);
        myCellState.TissueGeometryInitialization(sc.xDim, sc.yDim, sc.zDim, sc.d);
        myCellState.CellStateModelParameterSetup(testCell);
        myCellState.SetMisrepairRate("Epithelia", phys.k_error);
        myCellState.SetUpContactInhibition(true);
        myCellState.SetCheckpointEnabled(true);
        myCellState.SetSoftSaturationEnabled(true);

        for (int i = 0; i < cellNum; i++) {
            double cX = layout.GetCellPositionX(i) + sc.xDim / 2.0;
            double cY = layout.GetCellPositionY(i) + sc.yDim / 2.0;
            double cZ = layout.GetCellPositionZ(i) + sc.zDim / 2.0;
            myCellState.CellPositionInitialization(i, cX, cY, cZ);
            myCellState.CellTypeInitialiation(i, testCell);
            myCellState.CellPhaseInitializationRandom(i);
            myCellState.CellStateInitialization(i, "S1");
        }

        int initialDSB_val = static_cast<int>(meanDSB);
        double clock_phase = 0, clock_state = 0;
        int cycle_state = 0;

        for (int t = 0; t < sc.totalTimeStepNum; t++) {
            clock_phase += sc.T_microstep;
            clock_state += sc.T_microstep;

            if (clock_phase >= sc.deltaT_update) {
                map<int, string> phaseMap = myCellState.GetCellPhase();
                for (auto& cell : phaseMap)
                    myCellState.CellPhaseUpdate(cell.first, true, sc.deltaT_update, 1);
                clock_phase = 0;
            }

            if (clock_state >= sc.deltaT_update) {
                cycle_state++;
                map<int, string> stateMap = myCellState.GetCellState();
                for (auto& cell : stateMap) {
                    int dsb = (cycle_state == 1) ? initialDSB_val : 0;
                    myCellState.CellStateUpdate(cell.first, dsb, 0, sc.deltaT_update, 1);
                }
                clock_state = 0;
            }
        }

        map<int, string> finalState = myCellState.GetCellState();
        map<int, int> ancestry = myCellState.GetCellAncestryID();
        map<int, int> colonySizes;
        for (int i = 0; i < cellNum; i++) colonySizes[i] = 0;
        for (const auto& cell : finalState) {
            if (cell.second != "S1") continue;
            auto itA = ancestry.find(cell.first);
            if (itA != ancestry.end()) colonySizes[itA->second]++;
        }
        int colonies = 0;
        for (const auto& c : colonySizes)
            if (c.second >= 4) colonies++;

        sfStats.add(static_cast<double>(colonies) / static_cast<double>(cellNum));
    }

    return sfStats.mean();
}

//============================================================================
// Run a detailed validation for a given set of parameters
//============================================================================
struct ValidationResult {
    vector<double> doses, sf_lq, sf_sim, sf_sim_sd;
    double avg_log10_err = 0, max_log10_err = 0;
};

ValidationResult runDetailedValidation(
    const ParamsS2& params, double alpha_LQ, double beta_LQ,
    int cellNum, int nReps, unsigned long long baseSeed,
    const SimConfig& sc)
{
    ValidationResult vr;
    vr.doses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
    int nDoses = (int)vr.doses.size();
    vr.sf_lq.resize(nDoses);
    vr.sf_sim.resize(nDoses);
    vr.sf_sim_sd.resize(nDoses);

    const double sf_eps = std::max(1e-12, 0.5 / (double(cellNum) * double(nReps)));

    // Flatten dose x rep into parallel tasks for better thread utilization
    int nTasks = nDoses * nReps;
    vector<double> taskSF(nTasks, 0.0);

    #pragma omp parallel for schedule(dynamic)
    for (int flat = 0; flat < nTasks; flat++) {
        int di = flat / nReps;
        int ri = flat % nReps;
        taskSF[flat] = simulateSF(vr.doses[di], params, cellNum, 1,
                                  baseSeed + di * 100000ULL + ri * 1000ULL, sc);
    }

    // Aggregate per dose
    for (int di = 0; di < nDoses; di++) {
        RunningStats st;
        for (int ri = 0; ri < nReps; ri++)
            st.add(taskSF[di * nReps + ri]);
        vr.sf_lq[di] = LQ_survival(vr.doses[di], alpha_LQ, beta_LQ);
        vr.sf_sim[di] = st.mean();
        vr.sf_sim_sd[di] = st.sd();
    }

    double total = 0;
    for (int i = 0; i < nDoses; i++) {
        double err = fabs(log10(vr.sf_sim[i] + sf_eps) - log10(vr.sf_lq[i] + sf_eps));
        total += err;
        if (err > vr.max_log10_err) vr.max_log10_err = err;
    }
    vr.avg_log10_err = total / nDoses;
    return vr;
}

//============================================================================
// Parse a previously saved summary CSV to load calibrated parameters
//============================================================================
struct SavedCellParams {
    string cellLine;
    double alpha_LQ, beta_LQ;
    ParamsS2 anaParams;
    ParamsS2 simParams;
};

vector<SavedCellParams> loadParamsFromCSV(const string& csvPath) {
    vector<SavedCellParams> out;
    ifstream f(csvPath);
    if (!f.is_open()) return out;
    string header;
    getline(f, header);
    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string tok;
        vector<string> toks;
        while (getline(ss, tok, ',')) toks.push_back(tok);
        if (toks.size() < 21) continue;
        SavedCellParams p;
        p.cellLine = toks[0];
        p.alpha_LQ = stod(toks[1]);
        p.beta_LQ = stod(toks[2]);
        p.anaParams.e2 = stod(toks[3]);
        p.anaParams.e3 = stod(toks[4]);
        p.anaParams.a = stod(toks[5]);
        p.anaParams.k_error = stod(toks[6]);
        p.anaParams.T21 = stod(toks[7]);
        p.anaParams.T23 = stod(toks[8]);
        p.simParams.e2 = stod(toks[11]);
        p.simParams.e3 = stod(toks[12]);
        p.simParams.a = stod(toks[13]);
        p.simParams.k_error = stod(toks[14]);
        p.simParams.T21 = stod(toks[15]);
        p.simParams.T23 = stod(toks[16]);
        out.push_back(p);
    }
    return out;
}

//============================================================================
// Main
//============================================================================
int main(int argc, char** argv)
{
    string pideDir = "../PIDE3.4";
    string outputDir = "./sim_calibration_output";
    int maxCellTypes = 5;
    int coarseCells = 100;
    int coarseReps = 2;
    int nmMaxIter = 300;
    int detailedCells = 500;
    int detailedReps = 3;
    bool validationOnly = false;
    string paramsCSV = "";

    for (int i = 1; i < argc; i++) {
        string arg(argv[i]);
        if (arg == "--pide-dir" && i+1 < argc) pideDir = argv[++i];
        else if (arg == "--output-dir" && i+1 < argc) outputDir = argv[++i];
        else if (arg == "--max-cell-types" && i+1 < argc) maxCellTypes = atoi(argv[++i]);
        else if (arg == "--coarse-cells" && i+1 < argc) coarseCells = atoi(argv[++i]);
        else if (arg == "--coarse-reps" && i+1 < argc) coarseReps = atoi(argv[++i]);
        else if (arg == "--nm-max-iter" && i+1 < argc) nmMaxIter = atoi(argv[++i]);
        else if (arg == "--detailed-cells" && i+1 < argc) detailedCells = atoi(argv[++i]);
        else if (arg == "--detailed-reps" && i+1 < argc) detailedReps = atoi(argv[++i]);
        else if (arg == "--validation-only" && i+1 < argc) {
            validationOnly = true;
            paramsCSV = argv[++i];
        }
        else if (arg == "--help" || arg == "-h") {
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
                 << "  --pide-dir PATH         PIDE3.4 directory\n"
                 << "  --output-dir PATH       Output directory\n"
                 << "  --max-cell-types N      Max cell lines (default: 5)\n"
                 << "  --coarse-cells N        Cells for calibration sim (default: 100)\n"
                 << "  --coarse-reps N         Reps for calibration sim (default: 2)\n"
                 << "  --nm-max-iter N         Nelder-Mead iterations (default: 300)\n"
                 << "  --detailed-cells N      Cells for validation (default: 500)\n"
                 << "  --detailed-reps N       Reps for validation (default: 3)\n"
                 << "  --validation-only CSV   Skip calibration, load params from CSV\n";
            return 0;
        }
    }

    //========================================================================
    // Validation-only mode: load previously calibrated params and re-run sims
    //========================================================================
    if (validationOnly) {
        cout << "============================================\n"
             << "  Validation-Only Mode\n"
             << "============================================\n"
             << "  Params CSV: " << paramsCSV << "\n"
             << "  Detailed sim: " << detailedCells << " cells, " << detailedReps << " reps\n"
             << "============================================\n";

        auto saved = loadParamsFromCSV(paramsCSV);
        if (saved.empty()) {
            cerr << "ERROR: Could not load params from " << paramsCSV << "\n";
            return 1;
        }

        createDirectoryIfNotExists(outputDir.c_str());
        double wall_start = omp_get_wtime();
        SimConfig sc;
        unsigned long long masterSeed = 42ULL;

        double tau1 = 1.0 / (3.31 * sc.repairRateScale);
        double tau2 = 1.0 / (0.14 * sc.repairRateScale);
        double eff_repair = 0.62 * tau1 + 0.38 * tau2;

        ofstream summaryCSV(outputDir + "/sim_calib_summary.csv");
        summaryCSV << "CellLine,Alpha_LQ,Beta_LQ,"
                   << "Ana_e2,Ana_e3,Ana_a,Ana_kerr,Ana_T21,Ana_T23,Ana_AvgErr,Ana_MaxErr,"
                   << "Sim_e2,Sim_e3,Sim_a,Sim_kerr,Sim_T21,Sim_T23,Sim_NLL,Sim_Iter,Sim_Conv,Sim_AvgErr,Sim_MaxErr\n";

        ofstream paramsFile(outputDir + "/sim_calib_params.txt");
        paramsFile << "=== Validation-Only Results (" << detailedCells
                   << " cells, " << detailedReps << " reps) ===\n\n";

        for (size_t ci = 0; ci < saved.size(); ci++) {
            const auto& sp = saved[ci];
            unsigned long long cellSeed = masterSeed + ci * 10000000ULL;

            cout << "\n========================================\n"
                 << "  [" << ci+1 << "/" << saved.size() << "] " << sp.cellLine << "\n"
                 << "  alpha=" << sp.alpha_LQ << ", beta=" << sp.beta_LQ << "\n"
                 << "========================================\n";

            // Analytical params need T21/T23 set to effective repair time
            ParamsS2 anaForSim = sp.anaParams;
            anaForSim.T21 = eff_repair;
            anaForSim.T23 = eff_repair;

            cout << "  Running detailed validation (analytical params)...\n";
            ValidationResult vrAna = runDetailedValidation(
                anaForSim, sp.alpha_LQ, sp.beta_LQ, detailedCells, detailedReps,
                cellSeed + 5000000ULL, sc);

            cout << "  Running detailed validation (sim-calibrated params)...\n";
            ValidationResult vrSim = runDetailedValidation(
                sp.simParams, sp.alpha_LQ, sp.beta_LQ, detailedCells, detailedReps,
                cellSeed + 5000000ULL, sc);

            cout << "  Analytical: avg|dlog10SF|=" << vrAna.avg_log10_err
                 << ", max=" << vrAna.max_log10_err << "\n";
            cout << "  Sim-calib:  avg|dlog10SF|=" << vrSim.avg_log10_err
                 << ", max=" << vrSim.max_log10_err << "\n";

            string safeName = sanitizeFilename(sp.cellLine);
            ofstream cellCSV(outputDir + "/sim_calib_" + safeName + ".csv");
            cellCSV << "Dose,SF_LQ,SF_Analytical_Sim,SF_Analytical_Sim_SD,SF_SimCalib_Sim,SF_SimCalib_Sim_SD\n";
            for (size_t j = 0; j < vrAna.doses.size(); j++) {
                cellCSV << vrAna.doses[j] << ","
                        << vrAna.sf_lq[j] << ","
                        << vrAna.sf_sim[j] << "," << vrAna.sf_sim_sd[j] << ","
                        << vrSim.sf_sim[j] << "," << vrSim.sf_sim_sd[j] << "\n";
            }
            cellCSV.close();

            summaryCSV << sp.cellLine << "," << sp.alpha_LQ << "," << sp.beta_LQ << ","
                       << sp.anaParams.e2 << "," << sp.anaParams.e3 << ","
                       << sp.anaParams.a << "," << sp.anaParams.k_error << ","
                       << sp.anaParams.T21 << "," << sp.anaParams.T23 << ","
                       << vrAna.avg_log10_err << "," << vrAna.max_log10_err << ","
                       << sp.simParams.e2 << "," << sp.simParams.e3 << ","
                       << sp.simParams.a << "," << sp.simParams.k_error << ","
                       << sp.simParams.T21 << "," << sp.simParams.T23 << ","
                       << "0,0,reused,"
                       << vrSim.avg_log10_err << "," << vrSim.max_log10_err << "\n";

            paramsFile << "--- " << sp.cellLine << " ---\n"
                       << "  LQ: alpha=" << sp.alpha_LQ << ", beta=" << sp.beta_LQ << "\n"
                       << "  Analytical: e2=" << sp.anaParams.e2
                       << ", e3=" << sp.anaParams.e3
                       << ", a=" << sp.anaParams.a
                       << ", k_err=" << sp.anaParams.k_error
                       << ", T21=" << sp.anaParams.T21
                       << ", T23=" << sp.anaParams.T23 << "\n"
                       << "    -> avg|dlog10SF|=" << vrAna.avg_log10_err
                       << ", max=" << vrAna.max_log10_err << "\n"
                       << "  Sim-calib: e2=" << sp.simParams.e2
                       << ", e3=" << sp.simParams.e3
                       << ", a=" << sp.simParams.a
                       << ", k_err=" << sp.simParams.k_error
                       << ", T21=" << sp.simParams.T21
                       << ", T23=" << sp.simParams.T23 << "\n"
                       << "    -> avg|dlog10SF|=" << vrSim.avg_log10_err
                       << ", max=" << vrSim.max_log10_err << "\n\n";
        }

        summaryCSV.close();
        paramsFile.close();
        double wall_end = omp_get_wtime();
        cout << "\n============================================\n"
             << "  Total wall time: " << fixed << setprecision(1) << (wall_end - wall_start) << " s\n"
             << "  Results in: " << outputDir << "\n"
             << "============================================\n";
        return 0;
    }

    //========================================================================
    // Full calibration + validation mode
    //========================================================================
    cout << "============================================\n"
         << "  Simulation-Based Calibration\n"
         << "============================================\n"
         << "  Coarse sim: " << coarseCells << " cells, " << coarseReps << " reps\n"
         << "  NM max iter: " << nmMaxIter << "\n"
         << "  Detailed sim: " << detailedCells << " cells, " << detailedReps << " reps\n"
         << "  Max cell types: " << maxCellTypes << "\n"
         << "============================================\n";

    createDirectoryIfNotExists(outputDir.c_str());
    double wall_start = omp_get_wtime();

    // Load PIDE data
    PIDEDataReader reader;
    if (!reader.loadData(pideDir)) {
        cerr << "ERROR: Failed to load PIDE data from: " << pideDir << endl;
        return 1;
    }
    reader.printSummary();

    // Filter for 60Co experiments with valid LQ parameters (unique cell lines)
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

    if (co60Exps.empty()) {
        cerr << "ERROR: No 60Co experiments found!\n";
        return 1;
    }
    if ((int)co60Exps.size() > maxCellTypes)
        co60Exps.resize(maxCellTypes);

    cout << "\nCell lines to calibrate: " << co60Exps.size() << "\n";

    SimConfig sc;
    unsigned long long masterSeed = 42ULL;

    // Summary storage
    ofstream summaryCSV(outputDir + "/sim_calib_summary.csv");
    summaryCSV << "CellLine,Alpha_LQ,Beta_LQ,"
               << "Ana_e2,Ana_e3,Ana_a,Ana_kerr,Ana_T21,Ana_T23,Ana_AvgErr,Ana_MaxErr,"
               << "Sim_e2,Sim_e3,Sim_a,Sim_kerr,Sim_T21,Sim_T23,Sim_NLL,Sim_Iter,Sim_Conv,Sim_AvgErr,Sim_MaxErr\n";

    ofstream paramsFile(outputDir + "/sim_calib_params.txt");
    paramsFile << "=== Simulation-Based Calibration Results ===\n\n";

    for (size_t ci = 0; ci < co60Exps.size(); ci++) {
        const auto& exp = co60Exps[ci];
        const LQParameters& lq = exp.photonLQ_paper.isValid ? exp.photonLQ_paper : exp.photonLQ_fit;
        double alpha_LQ = lq.alpha;
        double beta_LQ = lq.beta;

        cout << "\n========================================\n"
             << "  [" << ci+1 << "/" << co60Exps.size() << "] " << exp.cellLine << "\n"
             << "  alpha=" << alpha_LQ << ", beta=" << beta_LQ << "\n"
             << "========================================\n";

        // Generate calibration data from LQ
        vector<double> calibDoses = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
        vector<DataPoint> calibData;
        for (double dose : calibDoses)
            calibData.push_back(DataPoint(dose, LQ_survival(dose, alpha_LQ, beta_LQ), 0.0));

        // Step 1: Analytical calibration (warm start)
        FitConfigS2 cfg;
        cfg.kappa = 40.0;
        cfg.dose_max_fit = 6.0;
        cfg.T_assay_h = 27.8;
        cfg.setTimescalesToCellCycle(sc.T_cellCycle);

        RepairMediatedMisrepairCalibrator analyticalCalib(cfg);
        analyticalCalib.setData(calibData);
        analyticalCalib.setDefaultUncertainties(0.1);

        ParamsS2 x0_ana(2.0, 10.0, 0.05, sc.T_cellCycle, sc.T_cellCycle, 1e-4);
        ParamsS2 step_ana(0.5, 2.0, 0.01, 1.0, 1.0, 5e-5);
        RepairMediatedMisrepairFitResult anaResult = analyticalCalib.fit(x0_ana, step_ana);

        cout << "  Analytical: e2=" << anaResult.params.e2
             << ", e3=" << anaResult.params.e3
             << ", a=" << anaResult.params.a
             << ", k_err=" << anaResult.params.k_error
             << " (" << (anaResult.converged ? "conv" : "FAIL") << ")\n";

        // Step 2: Simulation-based calibration
        cout << "  Starting simulation-based calibration...\n";

        unsigned long long cellSeed = masterSeed + ci * 10000000ULL;

        SimulationBasedCalibrator simCalib(cfg);
        simCalib.setData(calibData);
        simCalib.setDefaultUncertainties(0.1);

        simCalib.setBatchSurvivalFunction(
            [&sc, coarseCells, coarseReps, cellSeed](
                const vector<double>& doses, const ParamsS2& p) -> vector<double>
            {
                vector<double> results(doses.size());
                #pragma omp parallel for schedule(dynamic)
                for (size_t i = 0; i < doses.size(); i++) {
                    results[i] = simulateSF(doses[i], p, coarseCells, coarseReps,
                                            cellSeed + i * 77777ULL, sc);
                }
                return results;
            }
        );

        ParamsS2 x0_sim = anaResult.params;
        ParamsS2 step_sim(0.3, 1.0, 0.005, 2.0, 2.0, 3e-5);
        SimulationBasedFitResult simResult = simCalib.fit(x0_sim, step_sim, nmMaxIter);

        cout << "  Sim-calib: e2=" << simResult.params.e2
             << ", e3=" << simResult.params.e3
             << ", a=" << simResult.params.a
             << ", k_err=" << simResult.params.k_error
             << ", T21=" << simResult.params.T21
             << ", T23=" << simResult.params.T23
             << " (" << (simResult.converged ? "conv" : "not conv")
             << ", " << simResult.iterations << " iter, NLL=" << simResult.negLogLikelihood << ")\n";

        // Step 3: Detailed validation with both parameter sets
        cout << "  Running detailed validation (analytical params)...\n";
        ParamsS2 anaForSim = anaResult.params;
        double tau1 = 1.0 / (3.31 * sc.repairRateScale);
        double tau2 = 1.0 / (0.14 * sc.repairRateScale);
        double eff_repair = 0.62 * tau1 + 0.38 * tau2;
        anaForSim.T21 = eff_repair;
        anaForSim.T23 = eff_repair;
        ValidationResult vrAna = runDetailedValidation(
            anaForSim, alpha_LQ, beta_LQ, detailedCells, detailedReps,
            cellSeed + 5000000ULL, sc);

        cout << "  Running detailed validation (sim-calibrated params)...\n";
        ValidationResult vrSim = runDetailedValidation(
            simResult.params, alpha_LQ, beta_LQ, detailedCells, detailedReps,
            cellSeed + 5000000ULL, sc);

        cout << "  Analytical: avg|dlog10SF|=" << vrAna.avg_log10_err
             << ", max=" << vrAna.max_log10_err << "\n";
        cout << "  Sim-calib:  avg|dlog10SF|=" << vrSim.avg_log10_err
             << ", max=" << vrSim.max_log10_err << "\n";

        // Step 4: Write per-cell-line CSV
        string safeName = sanitizeFilename(exp.cellLine);
        ofstream cellCSV(outputDir + "/sim_calib_" + safeName + ".csv");
        cellCSV << "Dose,SF_LQ,SF_Analytical_Sim,SF_Analytical_Sim_SD,SF_SimCalib_Sim,SF_SimCalib_Sim_SD\n";
        for (size_t j = 0; j < vrAna.doses.size(); j++) {
            cellCSV << vrAna.doses[j] << ","
                    << vrAna.sf_lq[j] << ","
                    << vrAna.sf_sim[j] << "," << vrAna.sf_sim_sd[j] << ","
                    << vrSim.sf_sim[j] << "," << vrSim.sf_sim_sd[j] << "\n";
        }
        cellCSV.close();

        // Summary row
        summaryCSV << exp.cellLine << "," << alpha_LQ << "," << beta_LQ << ","
                   << anaResult.params.e2 << "," << anaResult.params.e3 << ","
                   << anaResult.params.a << "," << anaResult.params.k_error << ","
                   << anaResult.params.T21 << "," << anaResult.params.T23 << ","
                   << vrAna.avg_log10_err << "," << vrAna.max_log10_err << ","
                   << simResult.params.e2 << "," << simResult.params.e3 << ","
                   << simResult.params.a << "," << simResult.params.k_error << ","
                   << simResult.params.T21 << "," << simResult.params.T23 << ","
                   << simResult.negLogLikelihood << "," << simResult.iterations << ","
                   << (simResult.converged ? "yes" : "no") << ","
                   << vrSim.avg_log10_err << "," << vrSim.max_log10_err << "\n";

        paramsFile << "--- " << exp.cellLine << " ---\n"
                   << "  LQ: alpha=" << alpha_LQ << ", beta=" << beta_LQ << "\n"
                   << "  Analytical: e2=" << anaResult.params.e2
                   << ", e3=" << anaResult.params.e3
                   << ", a=" << anaResult.params.a
                   << ", k_err=" << anaResult.params.k_error
                   << ", T21=" << anaResult.params.T21
                   << ", T23=" << anaResult.params.T23 << "\n"
                   << "    -> avg|dlog10SF|=" << vrAna.avg_log10_err
                   << ", max=" << vrAna.max_log10_err << "\n"
                   << "  Sim-calib: e2=" << simResult.params.e2
                   << ", e3=" << simResult.params.e3
                   << ", a=" << simResult.params.a
                   << ", k_err=" << simResult.params.k_error
                   << ", T21=" << simResult.params.T21
                   << ", T23=" << simResult.params.T23 << "\n"
                   << "    -> avg|dlog10SF|=" << vrSim.avg_log10_err
                   << ", max=" << vrSim.max_log10_err << "\n\n";
    }

    summaryCSV.close();
    paramsFile.close();

    double wall_end = omp_get_wtime();
    cout << "\n============================================\n"
         << "  Total wall time: " << fixed << setprecision(1) << (wall_end - wall_start) << " s\n"
         << "  Results in: " << outputDir << "\n"
         << "============================================\n";

    return 0;
}
