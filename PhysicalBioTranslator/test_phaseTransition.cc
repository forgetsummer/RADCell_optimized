//============================================================================
// Test Program: Cell Phase Transition
//============================================================================
//
// Purpose:
//   This test program focuses on cell phase transitions (G1 -> S -> G2 -> M -> division).
//   It is adapted from the original testCellPhaseFirstRun.cc to work with the new
//   CellStateModel design that uses Cell objects instead of raw parameters.
//
// What it tests:
//   - Cell phase initialization (G1, S, G2, M, G0)
//   - Phase duration sampling from Gaussian distribution
//   - Phase transitions when age exceeds duration
//   - Cell division (mitosis) and daughter cell creation
//   - Contact inhibition effects on division
//   - Population growth over time
//
// Author: Ruirui Liu
//         Department of Nuclear Engineering and Radiation Health Physics
//         Oregon State University
//
// References:
//   - PhD Thesis: https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/rf55zf01n
//   - RADCell Paper: https://iopscience.iop.org/article/10.1088/1361-6560/abd4f9/meta
//
// Contact: liuruirui.nova@gmail.com
//
// Usage:
//   Compile with: make test_phaseTransition
//   Run with: ./test_phaseTransition
//
// Output:
//   cellPhase.csv - Phase distribution over time
//
//============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <sys/stat.h>
#include "Cell.hh"
#include "CellStateModel.hh"
#include "CellLayoutInitializer.hh"

using namespace std;

// Helper function to create output directory if it doesn't exist
void createDirectoryIfNotExists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        mkdir(path.c_str(), 0755);
        cout << "Created directory: " << path << endl;
    }
}

int main(int argc, char** argv)
{
    cout << "============================================" << endl;
    cout << "  Cell Phase Transition Test               " << endl;
    cout << "  (Adapted from testCellPhaseFirstRun.cc)  " << endl;
    cout << "============================================" << endl;

    // Create output directory
    createDirectoryIfNotExists("./phase_output");

    //------------------------------------------------------------------------
    // Simulation Parameters (from original test file)
    //------------------------------------------------------------------------
    double xDim = 1.0;    // unit in mm
    double yDim = 1.0;    // unit in mm
    double zDim = 0.0;    // unit in mm (2D monolayer)
    double d = 0.01;      // grid size, unit in mm
    
    // Time parameters
    double T = 10.0;                    // time step for one MCS, unit in second
    double deltaT_cellPhaseUpdate = 60; // unit in second
    int totalTimeStepNum = 60480;       // total time steps
    int cellNum = 1000;                 // number of cells
    
    // Cell cycle parameters (from original test)
    double tG1 = 9.0;       // mean G1 duration, unit in hour
    double sigmaG1 = 1.8;   // G1 std dev
    double tS = 11.0;       // mean S duration
    double sigmaS = 2.2;
    double tG2 = 1.0;       // mean G2 duration
    double sigmaG2 = 0.2;
    double tM = 1.0;        // mean M duration
    double sigmaM = 0.2;
    
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

    cout << "\n--- Simulation Parameters ---" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "  Grid size: " << d << " mm" << endl;
    cout << "  Cell number: " << cellNum << endl;
    cout << "  Total time steps: " << totalTimeStepNum << endl;
    cout << "  Time step: " << T << " seconds" << endl;
    cout << "  Phase update interval: " << deltaT_cellPhaseUpdate << " seconds" << endl;
    
    cout << "\n--- Cell Cycle Parameters ---" << endl;
    cout << "  G1: " << tG1 << " +/- " << sigmaG1 << " hours" << endl;
    cout << "  S:  " << tS << " +/- " << sigmaS << " hours" << endl;
    cout << "  G2: " << tG2 << " +/- " << sigmaG2 << " hours" << endl;
    cout << "  M:  " << tM << " +/- " << sigmaM << " hours" << endl;
    cout << "  Expected cycle: ~" << (tG1 + tS + tG2 + tM) << " hours" << endl;

    //------------------------------------------------------------------------
    // Create Cell Layout using CellLayoutInitializer
    //------------------------------------------------------------------------
    cout << "\n--- Initializing Cell Layout ---" << endl;
    
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d, d, d);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);
    cout << "The total cell number is " << layout.GetCellNumber() << endl;

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
    CellStateParameter statePara;
    statePara.E1 = E1;
    statePara.E2 = E2;
    statePara.E3 = E3;
    statePara.sigma = sigma;
    statePara.alpha = alpha;
    statePara.beta = beta;
    statePara.To12 = 24.0;  // observation window for S1->S2
    statePara.To13 = 24.0;  // observation window for S1->S3
    statePara.To21 = 24.0;  // observation window for S2->S1
    statePara.To23 = 24.0;  // observation window for S2->S3
    
    // DNA damage repair parameters
    CellDNADamageRepairParameter dnaRepairPara;
    dnaRepairPara.f1 = f1;
    dnaRepairPara.f2 = f2;
    dnaRepairPara.lambda1 = lambda1;
    dnaRepairPara.lambda2 = lambda2;
    
    // Bystander signal parameters
    CellBystanderSignalParameter bystanderPara;
    bystanderPara.Rt = 5000;
    bystanderPara.mu = 200;
    
    testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);

    //------------------------------------------------------------------------
    // Initialize CellStateModel (NEW API)
    //------------------------------------------------------------------------
    cout << "\n--- Initializing CellStateModel ---" << endl;
    
    double clock_cellPhaseUpdate = 0;
    int cycle_cellPhaseUpdate = 0;
    
    CellStateModel myCellState;
    
    // NEW: Use Cell object to setup parameters
    myCellState.CellStateModelParameterSetup(testCell);
    myCellState.TissueGeometryInitialization(xDim, yDim, zDim, d);
    myCellState.SetUpContactInhibition(true);  // Enable contact inhibition
    
    //------------------------------------------------------------------------
    // Initialize cell positions
    //------------------------------------------------------------------------
    cout << "\n--- Initializing Cell Positions ---" << endl;
    
    for (int i = 0; i < cellNum; i++) {
        double cX = layout.GetCellPositionX(i) + xDim / 2.0;
        double cY = layout.GetCellPositionY(i) + yDim / 2.0;
        double cZ = layout.GetCellPositionZ(i) + zDim / 2.0;
        myCellState.CellPositionInitialization(i, cX, cY, cZ);
    }

    //------------------------------------------------------------------------
    // Initialize cell types (NEW API - required before phase initialization)
    //------------------------------------------------------------------------
    for (int i = 0; i < cellNum; i++) {
        myCellState.CellTypeInitialiation(i, testCell);
    }

    //------------------------------------------------------------------------
    // Initialize cell phases randomly
    //------------------------------------------------------------------------
    cout << "\n--- Initializing Cell Phases Randomly ---" << endl;
    
    for (int i = 0; i < cellNum; i++) {
        // NEW API: CellPhaseInitializationRandom only takes cellID
        myCellState.CellPhaseInitializationRandom(i);
    }

    //------------------------------------------------------------------------
    // Initialize cell states (all start healthy)
    //------------------------------------------------------------------------
    for (int i = 0; i < cellNum; i++) {
        myCellState.CellStateInitialization(i, "S1");
    }

    //------------------------------------------------------------------------
    // Open output file
    //------------------------------------------------------------------------
    ofstream file;
    file.open("./phase_output/cellPhase.csv");
    file << "Step,G0,G1,S,G2,M,TotalCells" << endl;

    //------------------------------------------------------------------------
    // Main simulation loop
    //------------------------------------------------------------------------
    cout << "\n--- Starting Simulation ---" << endl;
    cout << "Total simulation time: " << (totalTimeStepNum * T / 3600.0) << " hours" << endl;
    
    clock_t start_time = clock();
    
    for (int i = 0; i < totalTimeStepNum; i++) {
        clock_cellPhaseUpdate = clock_cellPhaseUpdate + T;
        
        // Get current phase map (dynamic cell indices)
        std::map<int, std::string> phaseMap;
        phaseMap = myCellState.GetCellPhase();
        
        if (clock_cellPhaseUpdate >= deltaT_cellPhaseUpdate) {
            cycle_cellPhaseUpdate = cycle_cellPhaseUpdate + 1;
            
            // Count cells in each phase
            int G0Num = 0;
            int G1Num = 0;
            int SNum = 0;
            int G2Num = 0;
            int mitosisNum = 0;
            
            for (std::map<int, std::string>::iterator mitr_cell = phaseMap.begin(); 
                 mitr_cell != phaseMap.end(); mitr_cell++) {
                if (mitr_cell->second == "G0") G0Num++;
                else if (mitr_cell->second == "G1") G1Num++;
                else if (mitr_cell->second == "S") SNum++;
                else if (mitr_cell->second == "G2") G2Num++;
                else if (mitr_cell->second == "M") mitosisNum++;
            }
            
            // Write to file
            file << cycle_cellPhaseUpdate << "," << G0Num << "," << G1Num << "," 
                 << SNum << "," << G2Num << "," << mitosisNum << "," << phaseMap.size() << endl;
            
            // Update all cells' phases
            // NEW API: CellPhaseUpdate(cellID, proliferationState, deltaT, frequency)
            for (std::map<int, std::string>::iterator mitr_cell = phaseMap.begin(); 
                 mitr_cell != phaseMap.end(); mitr_cell++) {
                myCellState.CellPhaseUpdate(mitr_cell->first, true, deltaT_cellPhaseUpdate, 1);
            }
            
            clock_cellPhaseUpdate = 0;
            
            // Print progress every 100 cycles
            if (cycle_cellPhaseUpdate % 100 == 0) {
                double elapsed_hr = cycle_cellPhaseUpdate * deltaT_cellPhaseUpdate / 3600.0;
                cout << "  Step " << cycle_cellPhaseUpdate << " (t=" << fixed << setprecision(1) 
                     << elapsed_hr << " hr): " << phaseMap.size() << " cells" << endl;
            }
        }
    }
    
    clock_t end_time = clock();
    double elapsed_seconds = double(end_time - start_time) / CLOCKS_PER_SEC;
    
    file.close();
    
    //------------------------------------------------------------------------
    // Final Summary
    //------------------------------------------------------------------------
    cout << "\n--- Final Summary ---" << endl;
    
    std::map<int, std::string> phaseMap = myCellState.GetCellPhase();
    std::map<int, double> ageMap = myCellState.GetCellAge();
    std::map<int, double> durationMap = myCellState.GetCellPhaseDuration();
    
    // Count final phase distribution
    int G0Num = 0, G1Num = 0, SNum = 0, G2Num = 0, MNum = 0;
    for (auto& p : phaseMap) {
        if (p.second == "G0") G0Num++;
        else if (p.second == "G1") G1Num++;
        else if (p.second == "S") SNum++;
        else if (p.second == "G2") G2Num++;
        else if (p.second == "M") MNum++;
    }
    
    cout << "  Initial cells: " << cellNum << endl;
    cout << "  Final cells: " << phaseMap.size() << endl;
    cout << "  Growth factor: " << fixed << setprecision(2) 
         << (double)phaseMap.size() / cellNum << "x" << endl;
    cout << "\n  Phase distribution:" << endl;
    cout << "    G0: " << G0Num << " (" << 100.0 * G0Num / phaseMap.size() << "%)" << endl;
    cout << "    G1: " << G1Num << " (" << 100.0 * G1Num / phaseMap.size() << "%)" << endl;
    cout << "    S:  " << SNum << " (" << 100.0 * SNum / phaseMap.size() << "%)" << endl;
    cout << "    G2: " << G2Num << " (" << 100.0 * G2Num / phaseMap.size() << "%)" << endl;
    cout << "    M:  " << MNum << " (" << 100.0 * MNum / phaseMap.size() << "%)" << endl;
    
    cout << "\n  Simulation time: " << fixed << setprecision(2) << elapsed_seconds << " seconds" << endl;
    cout << "  Total phase update cycles: " << cycle_cellPhaseUpdate << endl;
    cout << "  Simulated time: " << (cycle_cellPhaseUpdate * deltaT_cellPhaseUpdate / 3600.0) << " hours" << endl;
    
    cout << "\n============================================" << endl;
    cout << "  Test completed successfully!              " << endl;
    cout << "  Output: ./phase_output/cellPhase.csv      " << endl;
    cout << "============================================" << endl;

    return 0;
}
