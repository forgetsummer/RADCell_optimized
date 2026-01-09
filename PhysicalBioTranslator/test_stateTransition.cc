//============================================================================
// Test Program: Cell State Transition
//============================================================================
//
// Purpose:
//   This test program focuses on cell state transitions (S1 <-> S2 <-> S3).
//   It is adapted from the original test_stateTransition.cc to work with the
//   new CellStateModel design that uses Cell objects instead of raw parameters.
//
// Cell States (NEW VERSION - 3 states only):
//   - S1: Healthy state (normal proliferating cell)
//   - S2: Stressed/damaged state (can recover to S1 or progress to S3)
//   - S3: Apoptotic/dead state (terminal, absorbing state)
//
// Note: The old version had S21 (early stressed) and S22 (late stressed) 
//       substates. The new version simplifies this to just S2.
//
// What it tests:
//   - Cell state initialization (S1, S2, S3)
//   - State transitions based on DNA damage (DSB count)
//   - State transitions based on bystander signal concentration
//   - State energy calculation: E = alpha * DSB + beta * integralConcentration
//   - Probabilistic state jumps (instantaneous and delayed)
//   - State distribution over time
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
//   Compile with: make test_stateTransition
//   Run with: ./test_stateTransition
//
// Output:
//   state_output/cellState.csv - State distribution over time
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
    cout << "  Cell State Transition Test               " << endl;
    cout << "  (Adapted from test_stateTransition.cc)   " << endl;
    cout << "============================================" << endl;

    // Create output directory
    createDirectoryIfNotExists("./state_output");

    //------------------------------------------------------------------------
    // Simulation Parameters
    //------------------------------------------------------------------------
    double xDim = 1.0;    // unit in mm
    double yDim = 1.0;    // unit in mm
    double zDim = 0.0;    // unit in mm (2D monolayer)
    double d = 0.01;      // grid size, unit in mm
    
    // Time parameters
    double T = 10.0;                     // time step for one MCS, unit in second
    double deltaT_cellStateUpdate = 60;  // unit in second
    int totalTimeStepNum = 10000;        // total time steps (same as test.cc)
    int cellNum = 1;                     // SINGLE CELL for trajectory test
    
    // Cell cycle parameters (shorter cycle for faster testing)
    double tG1 = 1.5;       // mean G1 duration, unit in hour
    double sigmaG1 = 0.25;  // G1 std dev
    double tS = 6.0;        // mean S duration
    double sigmaS = 0.25;
    double tG2 = 1.5;       // mean G2 duration
    double sigmaG2 = 0.25;
    double tM = 1.0;        // mean M duration
    double sigmaM = 0.25;
    
    // Cell state parameters
    double alpha = 0.2497;  // unit in 1/DSB (translates DSB to energy)
    double beta = 25.71;    // unit in ml/pg (translates concentration to energy)
    double E1 = 0.0;        // energy threshold for S1
    double E2 = 18.15168;   // energy threshold for S2
    double E3 = 48.47;      // energy threshold for S3
    double sigma = 6.96;    // sigma for state transition probability
    
    // Phase-dependent sensitivity factors
    double fG1 = 1.0;
    double fG2 = 1.015306122;
    double fM = 1.015306122;
    double fS = 0.816326;
    
    // DNA damage repair parameters
    double f1 = 0.62;       // fraction of fast repair
    double f2 = 0.38;       // fraction of slow repair
    double lambda1 = 3.31;  // fast repair rate, unit in 1/hour
    double lambda2 = 0.14;  // slow repair rate, unit in 1/hour

    // Test parameters for state transition
    int DSBNum = 0;  // Number of DNA double-strand breaks (NO radiation damage)
    double massOfBystanderSignal = 0;  // NO bystander signal
    double volume = d * d * d;
    double concentrationInMass = 0;  // NO external perturbation

    cout << "\n--- Simulation Parameters ---" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "  Grid size: " << d << " mm" << endl;
    cout << "  Cell number: " << cellNum << endl;
    cout << "  Total time steps: " << totalTimeStepNum << endl;
    cout << "  Time step: " << T << " seconds" << endl;
    cout << "  State update interval: " << deltaT_cellStateUpdate << " seconds" << endl;
    cout << "  Simulated time: " << (totalTimeStepNum * T / 3600.0) << " hours" << endl;
    
    cout << "\n--- Cell State Parameters ---" << endl;
    cout << "  alpha (DSB->Energy): " << alpha << " 1/DSB" << endl;
    cout << "  beta (Conc->Energy): " << beta << " ml/pg" << endl;
    cout << "  E1 (S1 threshold): " << E1 << endl;
    cout << "  E2 (S2 threshold): " << E2 << endl;
    cout << "  E3 (S3 threshold): " << E3 << endl;
    cout << "  sigma: " << sigma << endl;
    cout << "  Observation window: " << (tG1 + tS + tG2 + tM) << " hours (= tG1+tS+tG2+tM)" << endl;
    
    cout << "\n--- Test Stimulus ---" << endl;
    cout << "  DSB count: " << DSBNum << endl;
    cout << "  Bystander signal concentration: " << concentrationInMass << " pg/ml" << endl;
    cout << "  Expected energy from DSB: " << (alpha * DSBNum) << endl;
    cout << "  Expected energy from signal: " << (beta * concentrationInMass) << endl;

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
    // Observation window = total cell cycle time (tG1 + tS + tG2 + tM)
    double observationWindow = tG1 + tS + tG2 + tM;  // = 10 hours for this test
    
    CellStateParameter statePara;
    statePara.E1 = E1;
    statePara.E2 = E2;
    statePara.E3 = E3;
    statePara.sigma = sigma;
    statePara.alpha = alpha;
    statePara.beta = beta;
    statePara.To12 = observationWindow;  // observation window for S1->S2
    statePara.To13 = observationWindow;  // observation window for S1->S3
    statePara.To21 = observationWindow;  // observation window for S2->S1
    statePara.To23 = observationWindow;  // observation window for S2->S3
    
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
    
    double clock_cellStateUpdate = 0;
    int cycle_cellStateUpdate = 0;
    
    CellStateModel myCellState;
    
    // NEW: Use Cell object to setup parameters
    myCellState.CellStateModelParameterSetup(testCell);
    myCellState.TissueGeometryInitialization(xDim, yDim, zDim, d);
    myCellState.SetUpContactInhibition(false);  // Disable contact inhibition for state test
    
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
    // Initialize cell types (NEW API - required before phase/state initialization)
    //------------------------------------------------------------------------
    for (int i = 0; i < cellNum; i++) {
        myCellState.CellTypeInitialiation(i, testCell);
    }

    //------------------------------------------------------------------------
    // Initialize cell phases (random for realistic state testing)
    //------------------------------------------------------------------------
    cout << "\n--- Initializing Cell Phases (Random) ---" << endl;
    
    for (int i = 0; i < cellNum; i++) {
        myCellState.CellPhaseInitializationRandom(i);
    }

    //------------------------------------------------------------------------
    // Initialize cell states (all start healthy in S1)
    //------------------------------------------------------------------------
    cout << "\n--- Initializing Cell States to S1 (Healthy) ---" << endl;
    
    for (int i = 0; i < cellNum; i++) {
        myCellState.CellStateInitialization(i, "S1");
    }

    //------------------------------------------------------------------------
    // Open output file
    //------------------------------------------------------------------------
    ofstream file;
    file.open("./state_output/cellState.csv");
    file << "Step,S1,S2,S3,TotalCells" << endl;

    //------------------------------------------------------------------------
    // Main simulation loop - State Transition Test
    //------------------------------------------------------------------------
    cout << "\n--- Starting State Transition Simulation ---" << endl;
    cout << "Total simulation time: " << (totalTimeStepNum * T / 3600.0) << " hours" << endl;
    cout << "Applying constant stimulus: DSB=" << DSBNum << ", Concentration=" << concentrationInMass << endl;
    
    clock_t start_time = clock();
    
    for (int i = 0; i < totalTimeStepNum; i++) {
        clock_cellStateUpdate = clock_cellStateUpdate + T;
        
        if (clock_cellStateUpdate >= deltaT_cellStateUpdate) {
            cycle_cellStateUpdate = cycle_cellStateUpdate + 1;
            
            // Get current state map
            std::map<int, std::string> stateMap = myCellState.GetCellState();
            std::map<int, double> ageMap = myCellState.GetCellStateAge();
            std::map<int, double> durationMap = myCellState.GetCellStateDuration();
            
            // Count cells in each state (NEW VERSION: S1, S2, S3 only)
            int S1_num = 0;
            int S2_num = 0;
            int S3_num = 0;
            
            for (std::map<int, std::string>::iterator mitr_cell = stateMap.begin();
                 mitr_cell != stateMap.end(); mitr_cell++) {
                if (mitr_cell->second == "S1") {
                    S1_num++;
                } else if (mitr_cell->second == "S2") {
                    S2_num++;
                } else if (mitr_cell->second == "S3") {
                    S3_num++;
                }
            }
            
            // Write to file
            file << cycle_cellStateUpdate << "," << S1_num << "," << S2_num << ","
                 << S3_num << "," << stateMap.size() << endl;
            
            // Update all cells' states
            // NEW API: CellStateUpdate(cellID, DSBNum, integralConcentration, deltaT, frequency)
            for (std::map<int, std::string>::iterator mitr_cell = stateMap.begin();
                 mitr_cell != stateMap.end(); mitr_cell++) {
                myCellState.CellStateUpdate(mitr_cell->first, DSBNum, concentrationInMass, 
                                           deltaT_cellStateUpdate, 1);
            }
            
            clock_cellStateUpdate = 0;
            
            // Print progress every 1000 cycles
            if (cycle_cellStateUpdate % 1000 == 0) {
                double elapsed_hr = cycle_cellStateUpdate * deltaT_cellStateUpdate / 3600.0;
                cout << "  Step " << cycle_cellStateUpdate << " (t=" << fixed << setprecision(1) 
                     << elapsed_hr << " hr): S1=" << S1_num << ", S2=" << S2_num 
                     << ", S3=" << S3_num << endl;
            }
            
            // Print detailed info every 5000 cycles
            if (cycle_cellStateUpdate % 5000 == 0 && cycle_cellStateUpdate > 0) {
                cout << "\n  --- Detailed State Info at Step " << cycle_cellStateUpdate << " ---" << endl;
                int count = 0;
                for (auto& cell : stateMap) {
                    if (count < 5) {  // Print first 5 cells
                        cout << "    Cell " << cell.first << ": State=" << cell.second 
                             << ", Age=" << fixed << setprecision(2) << ageMap[cell.first]
                             << ", Duration=" << durationMap[cell.first] << endl;
                        count++;
                    }
                }
                cout << endl;
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
    
    std::map<int, std::string> stateMap = myCellState.GetCellState();
    std::map<int, double> ageMap = myCellState.GetCellStateAge();
    std::map<int, double> durationMap = myCellState.GetCellStateDuration();
    
    // Count final state distribution (NEW VERSION: S1, S2, S3 only)
    int S1_num = 0, S2_num = 0, S3_num = 0;
    for (auto& s : stateMap) {
        if (s.second == "S1") S1_num++;
        else if (s.second == "S2") S2_num++;
        else if (s.second == "S3") S3_num++;
    }
    
    int totalCells = stateMap.size();
    
    cout << "  Initial cells: " << cellNum << endl;
    cout << "  Final cells: " << totalCells << endl;
    
    cout << "\n  State distribution:" << endl;
    cout << "    S1 (Healthy):  " << S1_num << " (" << fixed << setprecision(1) 
         << 100.0 * S1_num / totalCells << "%)" << endl;
    cout << "    S2 (Stressed): " << S2_num << " (" 
         << 100.0 * S2_num / totalCells << "%)" << endl;
    cout << "    S3 (Apoptotic): " << S3_num << " (" 
         << 100.0 * S3_num / totalCells << "%)" << endl;
    
    cout << "\n  Simulation time: " << fixed << setprecision(2) << elapsed_seconds << " seconds" << endl;
    cout << "  Total state update cycles: " << cycle_cellStateUpdate << endl;
    cout << "  Simulated time: " << (cycle_cellStateUpdate * deltaT_cellStateUpdate / 3600.0) << " hours" << endl;
    
    cout << "\n============================================" << endl;
    cout << "  Test completed successfully!              " << endl;
    cout << "  Output: ./state_output/cellState.csv      " << endl;
    cout << "============================================" << endl;

    return 0;
}

