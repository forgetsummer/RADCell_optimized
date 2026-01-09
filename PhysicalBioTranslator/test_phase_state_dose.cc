/**
 * @file test_phase_state_dose.cc
 * @brief Test file for cell phase and state transitions under varying dose levels
 * 
 * This test simulates cell population response to different radiation doses.
 * For each dose level (0 to 99), it runs a full simulation and records the
 * final state distribution. This allows analysis of dose-response relationships.
 * 
 * Key features:
 * - Tests dose-response relationship for cell state transitions
 * - Outputs both phase and state information
 * - Uses the new 3-state model (S1, S2, S3)
 * - Adapted from old version's test_phase_state_dose.cc
 * 
 * @author Ruirui Liu (original), adapted for new version
 * @date 2025
 * 
 * Parameters:
 * - 100 cells in 1mm x 1mm 2D monolayer
 * - Dose levels: 0-99 (DSB = dose * 0.05 * 40 = 0 to 198 DSBs)
 * - Simulation time: ~16.8 hours per dose level
 * - Reaction rate constant: 2E-7 (1/#.s)
 * 
 * Output:
 * - cellState.csv: dose, S1, S2, S3 counts
 * - cellPhase.csv: dose, G1, S, G2, M, G0 counts
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <map>
#include <string>
#include <vector>
#include <sys/stat.h>

#include "Cell.hh"
#include "CellStateModel.hh"
#include "CellLayoutInitializer.hh"

using namespace std;

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
    cout << "  (New Version - 3 State Model)" << endl;
    cout << "============================================" << endl;

    //------------------------------------------------------------------------
    // Simulation parameters
    //------------------------------------------------------------------------
    double xDim = 1.0;      // unit in mm
    double yDim = 1.0;      // unit in mm
    double zDim = 0.0;      // unit in mm (2D monolayer)
    double d = 0.01;        // grid size, unit in mm
    
    // Reaction-diffusion parameters (updated reaction rate)
    double D = 1E-5;        // Diffusion coefficient, unit in mm^2/s
    double r = 2E-7;        // Reaction rate constant (1/#.s) - NEW VALUE
    int Rt = 5000;          // Total receptors
    double mu = 200;        // Signal production rate
    
    // Time parameters
    double deltaT_diffusion = 1;        // unit in second
    double deltaT_diffusionUpdate = 60; // unit in second
    double deltaT_cellPhaseUpdate = 60; // unit in second
    double deltaT_cellStateUpdate = 60; // unit in second
    double T = 10;                      // time step for one MCS, unit in second
    int totalTimeStepNum = 10000;       // total time steps (~27.8 hours) - reduced for faster testing
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
    
    // Build dose vector as in test.cc
    // Dose levels: 0.1, 0.2, ..., 1.0 Gy, then 1.2, 1.6, 2.0, ..., 8.0 Gy
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
    
    // DSB per Gy (typical value for gamma radiation)
    double DSB_per_Gy = 40.0;

    cout << "\n--- Simulation Parameters ---" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "  Grid size: " << d << " mm" << endl;
    cout << "  Cell number: " << cellNum << endl;
    cout << "  Total time steps per dose: " << totalTimeStepNum << endl;
    cout << "  Time step: " << T << " seconds" << endl;
    cout << "  Simulated time per dose: " << (totalTimeStepNum * T / 3600.0) << " hours" << endl;
    cout << "  Number of dose levels: " << prescribedDoseVec.size() << endl;
    cout << "  Dose range: " << prescribedDoseVec.front() << " to " << prescribedDoseVec.back() << " Gy" << endl;
    cout << "  DSB per Gy: " << DSB_per_Gy << endl;
    cout << "  Reaction rate constant: " << r << " (1/#.s)" << endl;

    //------------------------------------------------------------------------
    // Create Cell Layout
    //------------------------------------------------------------------------
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d, d, d);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);
    cout << "\nThe total cell number is " << layout.GetCellNumber() << endl;

    //------------------------------------------------------------------------
    // Create and configure Cell object (NEW API)
    //------------------------------------------------------------------------
    Cell testCell;
    testCell.CellConstruct("Epithelia", "Cytoplasm Nucleus", "Sphere", "Blue Green");
    
    // Set up cell cycle parameters
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
    
    // Set up cell state parameters
    // Use separate observation window widths for each state transition route
    CellStateParameter statePara;
    statePara.E1 = E1;
    statePara.E2 = E2;
    statePara.E3 = E3;
    statePara.sigma = sigma;
    statePara.alpha = alpha;
    statePara.beta = beta;
    statePara.To12 = 6.0;    // S1→S2: 6 hours
    statePara.To13 = 48.0;   // S1→S3: 48 hours
    statePara.To21 = 24.0;   // S2→S1: 24 hours
    statePara.To23 = 24.0;   // S2→S3: 24 hours
    
    // DNA damage repair parameters
    CellDNADamageRepairParameter dnaRepairPara;
    dnaRepairPara.f1 = f1;
    dnaRepairPara.f2 = f2;
    dnaRepairPara.lambda1 = lambda1;
    dnaRepairPara.lambda2 = lambda2;
    
    // Bystander signal parameters (not used in this test, but required)
    CellBystanderSignalParameter bystanderPara;
    bystanderPara.Rt = 5000;
    bystanderPara.mu = 200;
    
    testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);

    //------------------------------------------------------------------------
    // Create output directory and files
    //------------------------------------------------------------------------
    createDirectoryIfNotExists("./dose_output");
    
    ofstream filePhase, fileState;
    filePhase.open("./dose_output/cellPhase.csv");
    fileState.open("./dose_output/cellState.csv");
    
    // Write headers
    filePhase << "Dose,G1,S,G2,M,G0,TotalCells" << endl;
    fileState << "Dose,S1,S2,S3,TotalCells" << endl;

    //------------------------------------------------------------------------
    // Main dose loop
    //------------------------------------------------------------------------
    cout << "\n--- Starting Dose-Response Simulation ---" << endl;
    clock_t total_start = clock();
    
    for (size_t n = 0; n < prescribedDoseVec.size(); n++)
    {
        // Get dose and calculate DSB number
        double dose = prescribedDoseVec[n];
        int DSBNum = static_cast<int>(dose * DSB_per_Gy);
        
        cout << "PROCESSING DOSE LEVEL " << n << ": " << dose << " Gy (DSB = " << DSBNum << ")" << endl;
        
        // Reset clocks for this dose level
        double clock_cellPhaseUpdate = 0;
        double clock_cellStateUpdate = 0;
        int cycle_cellPhaseUpdate = 0;
        int cycle_cellStateUpdate = 0;
        
        //--------------------------------------------------------------------
        // Initialize CellStateModel for this dose level
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
        
        // Initialize cell phases - all start in G1
        for (int i = 0; i < cellNum; i++)
        {
            myCellState.CellPhaseInitialization(i, "G1");
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
                    int currentDSB = (cycle_cellStateUpdate == 1) ? DSBNum : 0;
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
        
        // Write to files (use dose instead of index n)
        int totalCells = finalPhaseMap.size();
        filePhase << dose << "," << G1_num << "," << S_num << "," << G2_num << "," 
                  << M_num << "," << G0_num << "," << totalCells << endl;
        fileState << dose << "," << S1_num << "," << S2_num << "," << S3_num << "," 
                  << totalCells << endl;
        
        // Print progress every 5 dose levels or at the end
        if (n % 5 == 0 || n == prescribedDoseVec.size() - 1)
        {
            cout << "  Dose " << dose << " Gy: Phases(G1=" << G1_num << ",S=" << S_num 
                 << ",G2=" << G2_num << ",M=" << M_num << ",G0=" << G0_num 
                 << ") States(S1=" << S1_num << ",S2=" << S2_num << ",S3=" << S3_num 
                 << ") Total=" << totalCells << endl;
        }
    }
    
    filePhase.close();
    fileState.close();
    
    clock_t total_end = clock();
    double total_time = double(total_end - total_start) / CLOCKS_PER_SEC;
    
    cout << "\n--- Simulation Complete ---" << endl;
    cout << "  Total time: " << fixed << setprecision(2) << total_time << " seconds" << endl;
    cout << "  Output files:" << endl;
    cout << "    - ./dose_output/cellPhase.csv" << endl;
    cout << "    - ./dose_output/cellState.csv" << endl;
    
    return 0;
}

