//============================================================================
// Test Program: PhysicalBioTranslator Module - Cell State Model
//============================================================================
//
// Purpose:
//   This is a standalone test program to verify the Cell and CellStateModel
//   classes functionality. It tests the cell state transition model which
//   translates physical quantities (DNA damage, bystander signals) to
//   biological responses (cell state transitions, cell cycle progression).
//
// What it tests:
//   - Cell class construction and parameter setting
//   - CellStateModel initialization
//   - Cell phase initialization and transitions (G1, S, G2, M, G0)
//   - Cell state initialization and transitions (S1, S2, S3)
//   - Cell position tracking
//   - Cell cycle parameter setup
//   - Cell state parameter setup
//
// Author: Ruirui Liu
//         Department of Nuclear Engineering and Radiation Health Physics
//         Oregon State University
//
// References:
//   - PhD Thesis: https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/rf55zf01n
//   - RADCell Paper: https://iopscience.iop.org/article/10.1088/1361-6560/abd4f9/meta
//   - Cell State Transition Theory: https://arxiv.org/abs/2501.11875
//
// Contact: liuruirui.nova@gmail.com
//
// Usage:
//   Compile with: make test_cellStateModel
//   Run with: ./test_cellStateModel
//
// Output:
//   Console output showing cell state and phase transitions over time
//
//============================================================================

#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <sys/stat.h>
#include "Cell.hh"
#include "CellStateModel.hh"

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
    cout << "  PhysicalBioTranslator Module Test        " << endl;
    cout << "  Cell State Model Test                    " << endl;
    cout << "============================================" << endl;
    
    // Create output directory
    createDirectoryIfNotExists("./cell_state_output");
    
    //------------------------------------------------------------------------
    // 1. Create and configure a Cell
    //------------------------------------------------------------------------
    cout << "\n--- Test 1: Cell Construction ---" << endl;
    
    Cell testCell;
    testCell.CellConstruct("Epithelia", "Cytoplasm Nucleus", "Sphere", "Blue Green");
    
    cout << "Cell Type: " << testCell.GetCellType() << endl;
    cout << "Cell Shape: " << testCell.GetCellShape() << endl;
    cout << "Cell Organelles: ";
    for (const auto& org : testCell.GetCellOrganelle()) {
        cout << org << " ";
    }
    cout << endl;
    
    //------------------------------------------------------------------------
    // 2. Set up cell parameters
    //------------------------------------------------------------------------
    cout << "\n--- Test 2: Cell Parameter Setup ---" << endl;
    
    // Cell cycle parameters (typical values for epithelial cells)
    CellCycleParameter cyclePara;
    cyclePara.mTG1 = 10.0;      // mean G1 duration (hours)
    cyclePara.sigmaTG1 = 2.0;   // G1 std dev
    cyclePara.mTS = 8.0;        // mean S duration
    cyclePara.sigmaTS = 1.5;
    cyclePara.mTG2 = 4.0;       // mean G2 duration
    cyclePara.sigmaTG2 = 1.0;
    cyclePara.mTM = 1.0;        // mean M duration
    cyclePara.sigmaTM = 0.2;
    cyclePara.fG1 = 1.0;        // radiation sensitivity factors
    cyclePara.fS = 0.5;
    cyclePara.fG2 = 1.2;
    cyclePara.fM = 1.5;
    
    // Cell state parameters
    CellStateParameter statePara;
    statePara.E1 = 0.0;         // S1 state energy
    statePara.E2 = 1.5;         // S2 state energy
    statePara.E3 = 3.0;         // S3 state energy
    statePara.sigma = 0.5;      // energy fluctuation
    statePara.alpha = 0.1;      // DSB to energy conversion
    statePara.beta = 0.01;      // bystander signal to energy conversion
    statePara.To12 = 24.0;      // observation window S1->S2 (hours)
    statePara.To13 = 24.0;      // observation window S1->S3
    statePara.To21 = 24.0;      // observation window S2->S1
    statePara.To23 = 24.0;      // observation window S2->S3
    
    // DNA damage repair parameters
    CellDNADamageRepairParameter dnaRepairPara;
    dnaRepairPara.f1 = 0.8;     // fast repair fraction
    dnaRepairPara.f2 = 0.2;     // slow repair fraction
    dnaRepairPara.lambda1 = 2.0; // fast repair rate (1/hour)
    dnaRepairPara.lambda2 = 0.2; // slow repair rate (1/hour)
    
    // Bystander signal parameters
    CellBystanderSignalParameter bystanderPara;
    bystanderPara.Rt = 100;     // receptor threshold
    bystanderPara.mu = 0.001;   // secretion rate
    
    testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);
    
    cout << "Cell cycle parameters set:" << endl;
    cout << "  G1: " << cyclePara.mTG1 << " +/- " << cyclePara.sigmaTG1 << " hours" << endl;
    cout << "  S:  " << cyclePara.mTS << " +/- " << cyclePara.sigmaTS << " hours" << endl;
    cout << "  G2: " << cyclePara.mTG2 << " +/- " << cyclePara.sigmaTG2 << " hours" << endl;
    cout << "  M:  " << cyclePara.mTM << " +/- " << cyclePara.sigmaTM << " hours" << endl;
    
    cout << "Cell state parameters set:" << endl;
    cout << "  E1 (healthy): " << statePara.E1 << endl;
    cout << "  E2 (stressed): " << statePara.E2 << endl;
    cout << "  E3 (lethal): " << statePara.E3 << endl;
    cout << "  sigma: " << statePara.sigma << endl;
    
    //------------------------------------------------------------------------
    // 3. Initialize CellStateModel
    //------------------------------------------------------------------------
    cout << "\n--- Test 3: CellStateModel Initialization ---" << endl;
    
    CellStateModel csm;
    
    // Set up tissue geometry (1mm x 1mm x 0mm = 2D monolayer)
    double xDim = 1.0;  // mm
    double yDim = 1.0;  // mm
    double zDim = 0.0;  // mm (2D)
    double gridSize = 0.05; // mm (50 um cell size)
    
    csm.TissueGeometryInitialization(xDim, yDim, zDim, gridSize);
    csm.SetUpContactInhibition(true);
    
    cout << "Tissue geometry: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "Grid size: " << gridSize << " mm" << endl;
    
    // Set up model parameters
    csm.CellStateModelParameterSetup(testCell);
    
    //------------------------------------------------------------------------
    // 4. Create a small population of cells
    //------------------------------------------------------------------------
    cout << "\n--- Test 4: Cell Population Initialization ---" << endl;
    
    int numCells = 5;
    
    for (int i = 0; i < numCells; i++) {
        // Initialize cell type
        csm.CellTypeInitialiation(i, testCell);
        
        // Initialize cell position (spread across the grid)
        double cx = 0.1 + i * 0.15;  // mm
        double cy = 0.5;             // mm
        double cz = 0.0;             // mm
        csm.CellPositionInitialization(i, cx, cy, cz);
        
        // Initialize cell phase randomly
        csm.CellPhaseInitializationRandom(i);
        
        // Initialize cell state (all start healthy)
        csm.CellStateInitialization(i, "S1");
    }
    
    cout << "Initialized " << numCells << " cells" << endl;
    
    // Print initial state
    cout << "\nInitial cell states:" << endl;
    cout << setw(8) << "CellID" << setw(10) << "Phase" << setw(10) << "State" << endl;
    cout << string(28, '-') << endl;
    
    auto phases = csm.GetCellPhase();
    auto states = csm.GetCellState();
    
    for (int i = 0; i < numCells; i++) {
        cout << setw(8) << i << setw(10) << phases[i] << setw(10) << states[i] << endl;
    }
    
    //------------------------------------------------------------------------
    // 5. Simulate cell state and phase transitions
    //------------------------------------------------------------------------
    cout << "\n--- Test 5: Cell State/Phase Transition Simulation ---" << endl;
    
    double deltaT = 60.0;      // time step in seconds
    int frequency = 1;         // update frequency
    int numSteps = 100;        // number of simulation steps (= 100 minutes)
    
    // Open output file
    ofstream outFile("./cell_state_output/cell_transitions.csv");
    outFile << "Time_min,CellID,Phase,State,PosX,PosY" << endl;
    
    cout << "Simulating " << numSteps << " time steps (" << numSteps * deltaT / 60.0 << " minutes)..." << endl;
    
    for (int step = 0; step <= numSteps; step++) {
        double time_min = step * deltaT / 60.0;
        
        // Get current states
        phases = csm.GetCellPhase();
        states = csm.GetCellState();
        auto posX = csm.GetCellPositionX();
        auto posY = csm.GetCellPositionY();
        
        // Write to file
        for (auto& p : phases) {
            int cellID = p.first;
            outFile << time_min << "," << cellID << "," << p.second << "," 
                    << states[cellID] << "," << posX[cellID] << "," << posY[cellID] << endl;
        }
        
        // Update each cell (no radiation, no bystander signal for this test)
        for (auto& p : phases) {
            int cellID = p.first;
            csm.CellPhaseUpdate(cellID, true, deltaT, frequency);
            csm.CellStateUpdate(cellID, 0, 0.0, deltaT, frequency);  // DSB=0, no bystander
        }
        
        // Print progress every 20 steps
        if (step % 20 == 0) {
            cout << "  Step " << step << "/" << numSteps << " (t=" << time_min << " min)" << endl;
        }
    }
    
    outFile.close();
    
    //------------------------------------------------------------------------
    // 6. Print final state
    //------------------------------------------------------------------------
    cout << "\n--- Test 6: Final Cell States ---" << endl;
    
    phases = csm.GetCellPhase();
    states = csm.GetCellState();
    auto ages = csm.GetCellAge();
    auto durations = csm.GetCellPhaseDuration();
    
    cout << setw(8) << "CellID" << setw(10) << "Phase" << setw(10) << "State" 
         << setw(12) << "Age(hr)" << setw(12) << "Duration" << endl;
    cout << string(52, '-') << endl;
    
    for (auto& p : phases) {
        int cellID = p.first;
        cout << setw(8) << cellID 
             << setw(10) << p.second 
             << setw(10) << states[cellID]
             << setw(12) << fixed << setprecision(3) << ages[cellID]
             << setw(12) << durations[cellID] << endl;
    }
    
    //------------------------------------------------------------------------
    // 7. Test with radiation damage
    //------------------------------------------------------------------------
    cout << "\n--- Test 7: Radiation Damage Response ---" << endl;
    
    // Simulate radiation exposure (add DSBs to cells)
    cout << "Applying radiation damage (10 DSBs per cell)..." << endl;
    
    for (auto& p : phases) {
        int cellID = p.first;
        // Apply 10 DSBs, no bystander signal
        csm.CellStateUpdate(cellID, 10, 0.0, deltaT, frequency);
    }
    
    // Check state changes
    states = csm.GetCellState();
    cout << "\nCell states after radiation:" << endl;
    cout << setw(8) << "CellID" << setw(10) << "State" << endl;
    cout << string(18, '-') << endl;
    
    for (auto& s : states) {
        cout << setw(8) << s.first << setw(10) << s.second << endl;
    }
    
    //------------------------------------------------------------------------
    // Summary
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  Test completed successfully!              " << endl;
    cout << "  Output saved to: ./cell_state_output/     " << endl;
    cout << "============================================" << endl;
    
    return 0;
}

