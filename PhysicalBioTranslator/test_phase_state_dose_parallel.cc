/**
 * @file test_phase_state_dose_parallel.cc
 * @brief Parallel version of test_phase_state_dose using OpenMP
 * 
 * This test simulates cell population response to different radiation doses
 * with OpenMP parallelization across dose levels for faster execution.
 * 
 * Key features:
 * - OpenMP parallel for loop across dose levels
 * - Thread-local CellStateModel instances
 * - Thread-safe output collection
 * - Significant speedup on multi-core systems
 * 
 * @author Ruirui Liu (original), parallelized version
 * @date 2025
 * 
 * Usage:
 *   export OMP_NUM_THREADS=8  # Set number of threads
 *   ./test_phase_state_dose_parallel
 * 
 * Output:
 * - cellState.csv: dose, S1, S2, S3 counts
 * - cellPhase.csv: dose, G1, S, G2, M, G0 counts
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <map>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <omp.h>  // OpenMP header

#include "Cell.hh"
#include "CellStateModel.hh"
#include "CellLayoutInitializer.hh"

using namespace std;

//============================================================================
// Class to read radiation transport results (like old version)
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

// Structure to store results for each dose level
struct DoseResult {
    double dose;
    int G1_num, S_num, G2_num, M_num, G0_num;
    int S1_num, S2_num, S3_num;
    int totalCells;
};

// Structure to store cell survival fraction results
struct CellSFResult {
    double dose;
    int colonyNumber;
    double sf;
};

// Helper function to create directory if it doesn't exist
void createDirectoryIfNotExists(const char* path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0755);
        cout << "Created directory: " << path << endl;
    }
}

int main(int argc, char** argv)
{
    cout << "============================================" << endl;
    cout << "  Test: Cell Phase & State vs Dose" << endl;
    cout << "  (Parallel Version with OpenMP)" << endl;
    cout << "============================================" << endl;

    //------------------------------------------------------------------------
    // Simulation parameters (matched to OLD VERSION test.cc)
    //------------------------------------------------------------------------
    double xDim = 5.0;      // unit in mm (OLD: 5mm, was 1mm)
    double yDim = 5.0;      // unit in mm (OLD: 5mm, was 1mm)
    double zDim = 0.0;      // unit in mm (2D monolayer)
    double d = 0.05;        // grid size, unit in mm (OLD: 0.05mm, was 0.01mm)
    
    // Reaction-diffusion parameters
    double D = 1E-5;        // Diffusion coefficient, unit in mm^2/s
    double r = 2E-7;        // Reaction rate constant (1/#.s)
    int Rt = 5000;          // Total receptors
    double mu = 200;        // Signal production rate
    
    // Time parameters
    double deltaT_diffusion = 1;        // unit in second
    double deltaT_diffusionUpdate = 60; // unit in second
    double deltaT_cellPhaseUpdate = 60; // unit in second
    double deltaT_cellStateUpdate = 60; // unit in second
    double T = 10;                      // time step for one MCS, unit in second
    int totalTimeStepNum = 10000;       // total time steps (~27.8 hours) (OLD: 10000, was 18000)
    int cellNum = 1000;                 // number of cells
    
    // Cell cycle parameters
    double tG1 = 1.5;       // mean G1 duration, unit in hour
    double sigmaG1 = 0.25;
    double tS = 6.0;        // mean S duration
    double sigmaS = 0.25;
    double tG2 = 1.5;       // mean G2 duration
    double sigmaG2 = 0.25;
    double tM = 1.0;        // mean M duration
    double sigmaM = 0.25;
    
    // Cell state parameters
    double alpha = 0.2497;  // unit in 1/DSB
    double beta = 25.71;    // unit in ml/pg
    double E1 = 0.0;
    double E2 = 18.15168;
    double E3 = 48.47;
    double sigma = 6.96;
    double fG1 = 1.0;
    double fG2 = 1.015306122;
    double fM = 1.015306122;
    double fS = 0.816326;
    double f1 = 0.62;
    double f2 = 0.38;
    double lambda1 = 3.31;  // unit in 1/hour
    double lambda2 = 0.14;  // unit in 1/hour
    
    // Total cell cycle time (observation window for state transitions)
    double T_cellCycle = tG1 + tS + tG2 + tM;  // unit in hour
    
    // Build dose vector as in test.cc
    int doseCaseNum = 20;
    int maxDose = 8;
    double step_dose = (double)maxDose / doseCaseNum;
    
    std::vector<double> prescribedDoseVec;
    for (int i = 0; i <= 10; i++) {
        prescribedDoseVec.push_back(0.1 * i);  // 0.0, 0.1, 0.2, ..., 1.0 Gy
    }
    for (int i = 3; i <= 20; i++) {
        prescribedDoseVec.push_back(i * step_dose);  // 1.2, 1.6, 2.0, ..., 8.0 Gy
    }
    
    // DSB per Gy (fallback if no MC data available)
    double DSB_per_Gy = 40.0;
    
    // Number of threads (can be overridden by OMP_NUM_THREADS env var)
    int numThreads = 10;
    omp_set_num_threads(numThreads);

    cout << "\n--- Simulation Parameters ---" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "  Grid size: " << d << " mm" << endl;
    cout << "  Cell number: " << cellNum << endl;
    cout << "  Total time steps per dose: " << totalTimeStepNum << endl;
    cout << "  Time step: " << T << " seconds" << endl;
    cout << "  Simulated time per dose: " << (totalTimeStepNum * T / 3600.0) << " hours" << endl;
    cout << "  Number of dose levels: " << prescribedDoseVec.size() << endl;
    cout << "  Dose range: " << prescribedDoseVec.front() << " to " << prescribedDoseVec.back() << " Gy" << endl;
    cout << "  DSB per Gy (fallback): " << DSB_per_Gy << endl;
    cout << "  Cell cycle time (T_cellCycle): " << T_cellCycle << " hours" << endl;
    cout << "  State transition window (To12=To13=To21=To23): " << T_cellCycle << " hours" << endl;
    cout << "  OpenMP threads: " << numThreads << endl;

    //------------------------------------------------------------------------
    // Create Cell Layout (shared across all threads - read only)
    //------------------------------------------------------------------------
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d, d, d);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);
    cout << "\nThe total cell number is " << layout.GetCellNumber() << endl;

    //------------------------------------------------------------------------
    // Load radiation transport data (like old version)
    // This provides cell-specific DSB values from Monte Carlo simulation
    //------------------------------------------------------------------------
    cout << "\n--- Loading Radiation Transport Data ---" << endl;
    RadiationTransportData radiationData;
    
    // File paths (adjust these to match your MC simulation output)
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
    createDirectoryIfNotExists("./dose_output");
    
    //------------------------------------------------------------------------
    // Prepare results storage (one entry per dose level)
    //------------------------------------------------------------------------
    std::vector<DoseResult> results(prescribedDoseVec.size());
    std::vector<CellSFResult> cellSF_results(prescribedDoseVec.size());

    //------------------------------------------------------------------------
    // Main parallel dose loop
    //------------------------------------------------------------------------
    cout << "\n--- Starting Parallel Dose-Response Simulation ---" << endl;
    double wall_start = omp_get_wtime();
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t n = 0; n < prescribedDoseVec.size(); n++)
    {
        // Get dose for this level
        double dose = prescribedDoseVec[n];
        
        // Thread-safe console output
        #pragma omp critical
        {
            cout << "Thread " << omp_get_thread_num() << " processing dose level " 
                 << n << ": " << dose << " Gy" << endl;
        }
        
        // Reset clocks for this dose level
        double clock_cellPhaseUpdate = 0;
        double clock_cellStateUpdate = 0;
        int cycle_cellPhaseUpdate = 0;
        int cycle_cellStateUpdate = 0;
        
        //--------------------------------------------------------------------
        // Create thread-local Cell object
        //--------------------------------------------------------------------
        Cell testCell;
        testCell.CellConstruct("Epithelia", "Cytoplasm Nucleus", "Sphere", "Blue Green");
        
        CellCycleParameter cyclePara;
        cyclePara.mTG1 = tG1;
        cyclePara.sigmaTG1 = sigmaG1;
        cyclePara.mTS = tS;
        cyclePara.sigmaTS = sigmaS;
        cyclePara.mTG2 = tG2;
        cyclePara.sigmaTG2 = sigmaG2;
        cyclePara.mTM = tM;
        cyclePara.sigmaTM = sigmaM;
        cyclePara.fG1 = fG1;
        cyclePara.fS = fS;
        cyclePara.fG2 = fG2;
        cyclePara.fM = fM;
        
        CellStateParameter statePara;
        statePara.E1 = E1;
        statePara.E2 = E2;
        statePara.E3 = E3;
        statePara.sigma = sigma;
        statePara.alpha = alpha;
        statePara.beta = beta;
        // Use same observation window width for all state transitions
        statePara.To12 = T_cellCycle;  // S1→S2: observation window = cell cycle time
        statePara.To13 = T_cellCycle;  // S1→S3: observation window = cell cycle time
        statePara.To21 = T_cellCycle;  // S2→S1: observation window = cell cycle time
        statePara.To23 = T_cellCycle;  // S2→S3: observation window = cell cycle time
        
        CellDNADamageRepairParameter dnaRepairPara;
        dnaRepairPara.f1 = f1;
        dnaRepairPara.f2 = f2;
        dnaRepairPara.lambda1 = lambda1;
        dnaRepairPara.lambda2 = lambda2;
        
        CellBystanderSignalParameter bystanderPara;
        bystanderPara.Rt = 5000;
        bystanderPara.mu = 200;
        
        testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);
        
        //--------------------------------------------------------------------
        // Initialize thread-local CellStateModel
        //--------------------------------------------------------------------
        CellStateModel myCellState;
        myCellState.TissueGeometryInitialization(xDim, yDim, zDim, d);
        myCellState.CellStateModelParameterSetup(testCell);
        myCellState.SetUpContactInhibition(true);
        
        // Initialize cell positions
        for (int i = 0; i < cellNum; i++)
        {
            double cX = layout.GetCellPositionX(i) + xDim / 2.0;
            double cY = layout.GetCellPositionY(i) + yDim / 2.0;
            double cZ = layout.GetCellPositionZ(i) + zDim / 2.0;
            myCellState.CellPositionInitialization(i, cX, cY, cZ);
        }
        
        // Initialize cell types
        for (int i = 0; i < cellNum; i++)
        {
            myCellState.CellTypeInitialiation(i, testCell);
        }
        
        // Initialize cell phases - RANDOM phases (like OLD version)
        // OLD version used: myCellState.CellPhaseInitializationRandom(i,1);
        for (int i = 0; i < cellNum; i++)
        {
            myCellState.CellPhaseInitializationRandom(i);
        }
        
        // Initialize cell states (all start in S1)
        for (int i = 0; i < cellNum; i++)
        {
            myCellState.CellStateInitialization(i, "S1");
        }

        //--------------------------------------------------------------------
        // Time loop for this dose level
        //--------------------------------------------------------------------
        for (int i = 0; i < totalTimeStepNum; i++)
        {
            clock_cellStateUpdate += T;
            clock_cellPhaseUpdate += T;
            
            // Get current phase map
            std::map<int, std::string> phaseMap = myCellState.GetCellPhase();
            
            // Cell phase update
            if (clock_cellPhaseUpdate >= deltaT_cellPhaseUpdate)
            {
                cycle_cellPhaseUpdate++;
                
                for (auto& cell : phaseMap)
                {
                    myCellState.CellPhaseUpdate(cell.first, true, deltaT_cellPhaseUpdate, 1);
                }
                clock_cellPhaseUpdate = 0;
            }
            
            // Cell state update
            if (clock_cellStateUpdate >= deltaT_cellStateUpdate)
            {
                cycle_cellStateUpdate++;
                
                std::map<int, std::string> stateMap = myCellState.GetCellState();
                
                for (auto& cell : stateMap)
                {
                    // Apply DSB only at the FIRST time step (cycle 1), no DSB afterwards
                    int currentDSB = 0;
                    if (cycle_cellStateUpdate == 1)
                    {
                        // Get cell-specific DSB from MC simulation data (like old version)
                        // Each cell gets DIFFERENT DSB based on its dose fraction
                        double DSB_cell = radiationData.GetAbsDSBOfCell(cell.first, dose);
                        currentDSB = static_cast<int>(DSB_cell);
                    }
                    myCellState.CellStateUpdate(cell.first, currentDSB, 0, deltaT_cellStateUpdate, 1);
                }
                clock_cellStateUpdate = 0;
            }
        }

        //--------------------------------------------------------------------
        // Collect final results for this dose level
        //--------------------------------------------------------------------
        std::map<int, std::string> finalPhaseMap = myCellState.GetCellPhase();
        std::map<int, std::string> finalStateMap = myCellState.GetCellState();
        
        // Count phases
        int G1_num = 0, S_num = 0, G2_num = 0, M_num = 0, G0_num = 0;
        for (auto& cell : finalPhaseMap)
        {
            if (cell.second == "G1") G1_num++;
            else if (cell.second == "S") S_num++;
            else if (cell.second == "G2") G2_num++;
            else if (cell.second == "M") M_num++;
            else if (cell.second == "G0") G0_num++;
        }
        
        // Count states
        int S1_num = 0, S2_num = 0, S3_num = 0;
        for (auto& cell : finalStateMap)
        {
            if (cell.second == "S1") S1_num++;
            else if (cell.second == "S2") S2_num++;
            else if (cell.second == "S3") S3_num++;
        }
        
        //--------------------------------------------------------------------
        // 3.2 Cell survival fraction calculation (colony formation assay)
        //--------------------------------------------------------------------
        std::map<int, int> cellAncestryIDMap = myCellState.GetCellAncestryID();
        
        // Initialize colony size map for each original cell
        std::map<int, int> cellColonySizeMap;
        for (int i = 0; i < cellNum; i++)
        {
            cellColonySizeMap[i] = 0;
        }
        
        // Count colony size for each original cell (only alive cells S1 or S21)
        for (auto& cell : finalStateMap)
        {
            // Only consider alive cells (S1 = healthy, S21 = recovering from stressed)
            if (cell.second == "S1" || cell.second == "S21")
            {
                int ancestryID = cellAncestryIDMap[cell.first];
                cellColonySizeMap[ancestryID] = cellColonySizeMap[ancestryID] + 1;
            }
        }
        
        // Count colonies that exceed the threshold
        int colonySizeThreshold = 2;  // Minimum colony size to count as surviving
        int colonyNumber = 0;
        
        for (auto& colony : cellColonySizeMap)
        {
            if (colony.second >= colonySizeThreshold)
            {
                colonyNumber = colonyNumber + 1;
            }
        }
        
        // Calculate survival fraction
        double survivalFraction = (double)colonyNumber / cellNum;
        
        // Store SF results (thread-safe since each thread writes to different index)
        cellSF_results[n].dose = dose;
        cellSF_results[n].colonyNumber = colonyNumber;
        cellSF_results[n].sf = survivalFraction;
        
        // Store results (thread-safe since each thread writes to different index)
        results[n].dose = dose;
        results[n].G1_num = G1_num;
        results[n].S_num = S_num;
        results[n].G2_num = G2_num;
        results[n].M_num = M_num;
        results[n].G0_num = G0_num;
        results[n].S1_num = S1_num;
        results[n].S2_num = S2_num;
        results[n].S3_num = S3_num;
        results[n].totalCells = finalPhaseMap.size();
        
        #pragma omp critical
        {
            cout << "  Dose " << dose << " Gy complete: S1=" << S1_num 
                 << ", S2=" << S2_num << ", S3=" << S3_num 
                 << ", Total=" << results[n].totalCells 
                 << ", SF=" << survivalFraction << endl;
        }
    }
    
    double wall_end = omp_get_wtime();
    double wall_time = wall_end - wall_start;
    
    //------------------------------------------------------------------------
    // Write results to files (sequential, after parallel section)
    //------------------------------------------------------------------------
    ofstream filePhase, fileState, fileSF;
    filePhase.open("./dose_output/cellPhase.csv");
    fileState.open("./dose_output/cellState.csv");
    fileSF.open("./dose_output/cellSurvival.csv");
    
    // Write headers
    filePhase << "Dose,G1,S,G2,M,G0,TotalCells" << endl;
    fileState << "Dose,S1,S2,S3,TotalCells" << endl;
    fileSF << "Dose,ColonyNumber,SurvivalFraction" << endl;
    
    // Write data (sorted by dose level)
    for (size_t n = 0; n < results.size(); n++)
    {
        filePhase << results[n].dose << "," << results[n].G1_num << "," 
                  << results[n].S_num << "," << results[n].G2_num << "," 
                  << results[n].M_num << "," << results[n].G0_num << "," 
                  << results[n].totalCells << endl;
        fileState << results[n].dose << "," << results[n].S1_num << "," 
                  << results[n].S2_num << "," << results[n].S3_num << "," 
                  << results[n].totalCells << endl;
        fileSF << cellSF_results[n].dose << "," << cellSF_results[n].colonyNumber << "," 
               << cellSF_results[n].sf << endl;
    }
    
    filePhase.close();
    fileState.close();
    fileSF.close();
    
    cout << "\n--- Simulation Complete ---" << endl;
    cout << "  Wall clock time: " << fixed << setprecision(2) << wall_time << " seconds" << endl;
    cout << "  Speedup estimate: ~" << setprecision(1) << (prescribedDoseVec.size() * 1.0 / (wall_time / (totalTimeStepNum * T / 3600.0 * 60))) << "x" << endl;
    cout << "  Output files:" << endl;
    cout << "    - ./dose_output/cellPhase.csv" << endl;
    cout << "    - ./dose_output/cellState.csv" << endl;
    cout << "    - ./dose_output/cellSurvival.csv" << endl;
    
    return 0;
}

