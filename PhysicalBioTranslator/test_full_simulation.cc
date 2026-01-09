/**
 * @file test_full_simulation.cc
 * @brief Full cell biology simulation test - adapted from test.cc
 * 
 * This test is adapted from the original test.cc but uses the new CellStateModel API.
 * It simulates:
 * - Cell phase transitions (G1, S, G2, M, G0)
 * - Cell state transitions (S1, S2, S3)
 * - Dose-response relationship
 * - Colony formation survival
 * 
 * This is a standalone version that doesn't require Geant4 radiation transport.
 * DSB numbers are calculated directly from dose (DSB = dose * DSB_per_Gy).
 * 
 * Key differences from original test.cc:
 * - Uses new Cell-based parameter setup API
 * - No Geant4 radiation transport (DSB calculated from dose directly)
 * - Uses new 3-state model (S1, S2, S3) instead of 4-substate (S1, S21, S22, S3)
 * 
 * @author Adapted from test.cc by Ruirui Liu
 * @date January 2026
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <omp.h>

#include "Cell.hh"
#include "CellStateModel.hh"
#include "CellLayoutInitializer.hh"

using namespace std;

struct cellSFObject
{
    double dose;
    int colonyNumber;
    double sf;
};

// Helper function to create directory
void createDirectoryIfNotExists(const string& path) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) {
        mkdir(path.c_str(), 0755);
    }
}

int main(int argc, char** argv)
{
    clock_t t1, t2;
    t1 = clock();

    cout << "============================================" << endl;
    cout << "  Full Cell Biology Simulation Test" << endl;
    cout << "  (Adapted from test.cc - No Geant4)" << endl;
    cout << "============================================" << endl;

//*******************************************************************************
//**************************** 1, Start Simulation parameters
//*******************************************************************************
    double xDim = 5;    // unit in mm
    double yDim = 5;    // unit in mm
    double zDim = 0;    // unit in mm
    double d = 0.05;    // unit in mm, dimension of cell home
    double d_diffusion = 0.05;  // unit in mm, spatial step for diffusion
    double D = 1E-5;    // unit in mm^2/s
    double r = 2E-7;    // unit in 1/#.s
    int Rt = 5000;      // unit in #
    double mu = 2000 * 10;  // unit in #
    double deltaT_diffusion = 1;        // unit in second
    double deltaT_diffusionUpdate = 60; // unit in second
    double deltaT_cellPhaseUpdate = 60; // unit in second
    double deltaT_cellStateUpdate = 60; // unit in second
    int cellType = 1;
    double tG1 = 1.5;   // unit in hour
    double sigmaG1 = 0.25;  // unit in hour
    double tS = 6;
    double sigmaS = 0.25;
    double tG2 = 1.5;
    double sigmaG2 = 0.25;
    double tM = 1;
    double sigmaM = 0.25;
    double alpha = 0.2497; // unit in 1/DSB

    double beta = 25.71;    // unit in ml/pg
    double E1 = 0;
    double E2 = 18.15168;
    double E3 = 48.47;
    double sigma = 6.96;
    double fG1 = 1;
    double fG2 = 1.015306122;
    double fM = 1.015306122;
    double fS = 0.816326;
    double f1 = 0.62;
    double f2 = 0.38;
    double lambda1 = 3.31;  // unit in 1/hour
    double lambda2 = 0.14;  // unit in 1/hour
    double T = 10;          // time step for one MCS, unit in second
    int totalTimeStepNum = 10000;
    int cellNum = 1000;     // unit in number

    bool considerRIBE = false;

    // Calculate cell state model parameters using the same method as original test.cc
    CellStateModel cstate;
    double p_sp = 0.001;
    double p_mis_sp = 0.026;
    double avgN = 194.05;
    CellStateModelParameter cpr = cstate.GetCellStateModelParameters(p_sp, p_mis_sp, avgN);
    
    cout << "The p_sp = " << p_sp << endl;
    cout << "The p_mis_sp = " << p_mis_sp << endl;
    cout << "The alpha = " << cpr.alpha << endl;
    cout << "The beta = " << cpr.beta << endl;
    cout << "The E1 = " << cpr.E1 << endl;
    cout << "The E2 = " << cpr.E2 << endl;
    cout << "The E3 = " << cpr.E3 << endl;
    cout << "The sigma = " << cpr.sigma << endl;

    // Use calculated parameters 
    // NOTE: Original test.cc multiplies beta by 5, but this makes cells too sensitive
    // We use the original beta value to match the working test_phase_state_dose_parallel.cc
    E2 = cpr.E2;
    alpha = cpr.alpha;
    // beta = cpr.beta * 5;  // Original test.cc - too sensitive
    // Keep beta = 25.71 as defined above to match working test

    string DSBUpperOrLowerBoundType = "random";

//*********************************************************************************
//**********************End simulation parameters
//*************************************************************************************

//*************************************************************************************
//**********************2, Dose levels (from test.cc)
//*************************************************************************************
    int doseCaseNum = 20;
    int maxDose = 8;
    double step_dose = (double)maxDose / doseCaseNum;
    cout << step_dose << endl;

    vector<double> prescribedDoseVec;
    for (int i = 0; i <= 10; i++) {
        prescribedDoseVec.push_back(0.1 * i);
    }
    for (int i = 3; i <= 20; i++) {
        prescribedDoseVec.push_back(i * step_dose);
    }

    // Cell layout initialization
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d, d, d);  // cell home size is in unit of mm
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);  // unit here is also mm
    cout << "The total cell number is " << layout.GetCellNumber() << endl;

    cout << "\n--- Simulation Parameters ---" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "  Cell home size: " << d << " mm" << endl;
    cout << "  Cell number: " << cellNum << endl;
    cout << "  Total time steps: " << totalTimeStepNum << endl;
    cout << "  Time step: " << T << " seconds" << endl;
    cout << "  Simulated time: " << (totalTimeStepNum * T / 3600.0) << " hours" << endl;
    cout << "  Number of dose levels: " << prescribedDoseVec.size() << endl;
    cout << "  Alpha: " << alpha << endl;
    cout << "  Beta: " << beta << endl;

    // Create output directory
    string outputDir = "./full_simulation_output/";
    createDirectoryIfNotExists(outputDir);

    string cellSurvivalFileName = outputDir + "cellSurvival.csv";
    ofstream file_sf;
    file_sf.open(cellSurvivalFileName.c_str());

//*******************************************************************************************************************
//********************************3, Start Cell biology simulation
//*******************************************************************************************************************

//*****************************************************************************************************************
//*********************************3.1 Cell system CellPhaseInitialization
//*****************************************************************************************************************
    omp_set_num_threads(10);  // set up the thread number for parallel processing
    vector<cellSFObject> cellSF_private(prescribedDoseVec.size());

    #pragma omp parallel for
    for (int n = 0; n < (int)prescribedDoseVec.size(); n++)  // simulation for different prescribed dose
    {
        ofstream file_phase;
        string cellPhaseFileName;
        stringstream ss;
        ss << prescribedDoseVec[n];
        string cellDose = ss.str();
        cellPhaseFileName = outputDir + "cellPhase_dose_" + cellDose + ".csv";
        file_phase.open(cellPhaseFileName.c_str());

        ofstream file_state;
        string cellStateFileName;
        cellStateFileName = outputDir + "cellState_dose_" + cellDose + ".csv";
        file_state.open(cellStateFileName.c_str());

        double clock_cellPhaseUpdate = 0;
        double clock_cellStateUpdate = 0;
        double clock_diffusionUpdate = 0;
        int cycle_cellPhaseUpdate = 0;
        int cycle_cellStateUpdate = 0;
        int cycle_diffusionUpdate = 0;

        #pragma omp critical
        {
            cout << "PROCESSING DOSE STEP: " << n << " (Dose = " << prescribedDoseVec[n] << " Gy)" << endl;
        }

        CellStateModel myCellState;

        // ===== NEW API: Set up cell parameters using Cell object =====
        Cell testCell;
        testCell.CellConstruct("Epithelia", "Cytoplasma Nucleus", "Sphere", "Blue Green");

        CellCycleParameter cyclePara;
        cyclePara.mTG1 = tG1;
        cyclePara.sigmaTG1 = sigmaG1;
        cyclePara.mTS = tS;
        cyclePara.sigmaTS = sigmaS;
        cyclePara.mTG2 = tG2;
        cyclePara.sigmaTG2 = sigmaG2;
        cyclePara.mTM = tM;
        cyclePara.sigmaTM = sigmaM;

        CellStateParameter statePara;
        statePara.alpha = alpha;
        statePara.beta = beta;
        // Use observation windows from working test_phase_state_dose_parallel.cc
        statePara.To12 = 2.0;    // S1→S2: 2 hours (fast stress response)
        statePara.To13 = 8.0;    // S1→S3: 8 hours (fast direct death)
        statePara.To21 = 12.0;   // S2→S1: 12 hours (recovery)
        statePara.To23 = 4.0;    // S2→S3: 4 hours (fast stressed cell death)

        CellDNADamageRepairParameter dnaRepairPara;
        dnaRepairPara.f1 = f1;
        dnaRepairPara.f2 = f2;
        dnaRepairPara.lambda1 = lambda1;
        dnaRepairPara.lambda2 = lambda2;

        CellBystanderSignalParameter bystanderPara;
        bystanderPara.Rt = Rt;
        bystanderPara.mu = mu;

        testCell.SetCellParametersForCellStateModeling(cyclePara, statePara, dnaRepairPara, bystanderPara);

        // ===== NEW API: Initialize CellStateModel with Cell object =====
        myCellState.TissueGeometryInitialization(xDim, yDim, zDim, d);
        myCellState.CellStateModelParameterSetup(testCell);

        // Initialize cell positions
        for (int i = 0; i < layout.GetCellNumber(); i++) {
            double cX = layout.GetCellPositionX(i) + xDim / 2.0;
            double cY = layout.GetCellPositionY(i) + yDim / 2.0;
            double cZ = layout.GetCellPositionZ(i) + zDim / 2.0;
            myCellState.CellPositionInitialization(i, cX, cY, cZ);
        }

        // ===== NEW API: Initialize cell types =====
        for (int i = 0; i < layout.GetCellNumber(); i++) {
            myCellState.CellTypeInitialiation(i, testCell);
        }

        // ===== NEW API: Initialize cell phases (G1 to match working test) =====
        for (int i = 0; i < layout.GetCellNumber(); i++) {
            myCellState.CellPhaseInitialization(i, "G1");  // Start all cells in G1 phase
        }

        // Initialize cell states
        for (int i = 0; i < layout.GetCellNumber(); i++) {
            myCellState.CellStateInitialization(i, "S1");
        }

        // Time loop
        for (int i = 0; i < totalTimeStepNum; i++) {
            clock_cellStateUpdate = clock_cellStateUpdate + T;
            clock_cellPhaseUpdate = clock_cellPhaseUpdate + T;
            clock_diffusionUpdate = clock_diffusionUpdate + T;

            map<int, string> phaseMap;
            phaseMap = myCellState.GetCellPhase();

            // Cell phase update
            if (clock_cellPhaseUpdate >= deltaT_cellPhaseUpdate) {
                cycle_cellPhaseUpdate = cycle_cellPhaseUpdate + 1;

                map<int, double> ageMap;
                map<int, double> durationMap;
                ageMap = myCellState.GetCellAge();
                durationMap = myCellState.GetCellPhaseDuration();

                int mitosisNum = 0;
                int G0Num = 0;
                int G1Num = 0;
                int G2Num = 0;
                int SNum = 0;

                for (map<int, string>::iterator mitr_cell = phaseMap.begin(); mitr_cell != phaseMap.end(); mitr_cell++) {
                    if (mitr_cell->second == "G0") {
                        G0Num = G0Num + 1;
                    }
                    if (mitr_cell->second == "M") {
                        mitosisNum = mitosisNum + 1;
                    }
                    if (mitr_cell->second == "G1") {
                        G1Num = G1Num + 1;
                    }
                    if (mitr_cell->second == "G2") {
                        G2Num = G2Num + 1;
                    }
                    if (mitr_cell->second == "S") {
                        SNum = SNum + 1;
                    }
                }
                file_phase << cycle_cellPhaseUpdate << "," << G0Num << "," << G1Num << "," << SNum << "," << G2Num << "," << mitosisNum << "," << phaseMap.size() << endl;

                // ===== NEW API: CellPhaseUpdate without cellType argument =====
                for (map<int, string>::iterator mitr_cell = phaseMap.begin(); mitr_cell != phaseMap.end(); mitr_cell++) {
                    myCellState.CellPhaseUpdate(mitr_cell->first, true, deltaT_cellPhaseUpdate, 1);
                }
                clock_cellPhaseUpdate = 0;
            }

            // Cell state update
            if (clock_cellStateUpdate >= deltaT_cellStateUpdate) {
                cycle_cellStateUpdate = cycle_cellStateUpdate + 1;

                map<int, string> stateMap;
                map<int, double> durationMap;
                map<int, double> ageMap;
                stateMap = myCellState.GetCellState();
                durationMap = myCellState.GetCellStateDuration();
                ageMap = myCellState.GetCellStateAge();

                map<int, double> cellPositionX;
                map<int, double> cellPositionY;
                map<int, double> cellPositionZ;
                cellPositionX = myCellState.GetCellPositionX();
                cellPositionY = myCellState.GetCellPositionY();
                cellPositionZ = myCellState.GetCellPositionZ();

                // Count cells in each state (new 3-state model: S1, S2, S3)
                int S1_num = 0;
                int S2_num = 0;
                int S3_num = 0;

                for (map<int, string>::iterator mitr_cell = stateMap.begin(); mitr_cell != stateMap.end(); mitr_cell++) {
                    if (stateMap[mitr_cell->first] == "S1") {
                        S1_num = S1_num + 1;
                    }
                    if (stateMap[mitr_cell->first] == "S2") {
                        S2_num = S2_num + 1;
                    }
                    if (stateMap[mitr_cell->first] == "S3") {
                        S3_num = S3_num + 1;
                    }
                }

                file_state << cycle_cellStateUpdate << "," << S1_num << "," << S2_num << "," << S3_num << "," << stateMap.size() << endl;

                // ===== NEW API: CellStateUpdate without cellType argument =====
                for (map<int, string>::iterator mitr_cell = stateMap.begin(); mitr_cell != stateMap.end(); mitr_cell++) {
                    double cX = cellPositionX.at(mitr_cell->first);
                    double cY = cellPositionY.at(mitr_cell->first);
                    double cZ = cellPositionZ.at(mitr_cell->first);

                    double concentrationInMass = 0;
                    // Note: RIBE not implemented in this standalone version

                    if (cycle_cellStateUpdate == 1) {
                        // Calculate DSB from dose (approximation without MC transport)
                        // Using 80 DSB per Gy to match working test_phase_state_dose_parallel.cc
                        double DSB_cell = prescribedDoseVec[n] * 80;
                        myCellState.CellStateUpdate(mitr_cell->first, (int)DSB_cell, 0, deltaT_cellStateUpdate, 1);
                    } else {
                        myCellState.CellStateUpdate(mitr_cell->first, 0, concentrationInMass, deltaT_cellStateUpdate, 1);
                    }
                }

                clock_cellStateUpdate = 0;
            }
        }

//*********************************************************************************************************
//******************************************3.2 Cell cell surival fraction calculation
//********************************************************************************************************
        map<int, string> stateMap;
        map<int, double> durationMap;
        map<int, double> ageMap;
        stateMap = myCellState.GetCellState();
        ageMap = myCellState.GetCellStateAge();
        durationMap = myCellState.GetCellStateDuration();
        map<int, int> cellAncestryIDMap;
        cellAncestryIDMap = myCellState.GetCellAncestryID();
        map<int, int> cellColonySizeMap;

        for (int i = 0; i < cellNum; i++) {
            cellColonySizeMap[i] = 0;
        }

        // Count alive cells for colony formation (S1 only in new model)
        for (map<int, string>::iterator mitr_cell = stateMap.begin(); mitr_cell != stateMap.end(); mitr_cell++) {
            if (mitr_cell->second == "S1") {  // Only S1 is alive in new 3-state model
                cellColonySizeMap[cellAncestryIDMap[mitr_cell->first]] = cellColonySizeMap[cellAncestryIDMap[mitr_cell->first]] + 1;
            }
        }

        int colonySizeThreshold = 2;
        int colonyNumber = 0;

        for (map<int, int>::iterator mitr_cell = cellColonySizeMap.begin(); mitr_cell != cellColonySizeMap.end(); mitr_cell++) {
            if (mitr_cell->second >= colonySizeThreshold) {
                colonyNumber = colonyNumber + 1;
            }
        }

        cellSFObject runSFRes;
        runSFRes.dose = prescribedDoseVec[n];
        runSFRes.colonyNumber = stateMap.size();
        runSFRes.sf = (double)colonyNumber / cellNum;
        cellSF_private[n] = runSFRes;

        #pragma omp critical
        {
            cout << "  Dose " << prescribedDoseVec[n] << " Gy complete: Total cells=" << stateMap.size()
                 << ", Colonies=" << colonyNumber << ", SF=" << runSFRes.sf << endl;
        }

        file_phase.close();
        file_state.close();
    }

    cout << "finish dose loop" << endl;

    // Write survival results
    for (int i = 0; i < (int)cellSF_private.size(); i++) {
        file_sf << cellSF_private[i].dose << "," << cellSF_private[i].colonyNumber << "," << cellSF_private[i].sf << endl;
    }
    file_sf.close();

    t2 = clock();
    float diff((float)t2 - (float)t1);
    cout << "The processing time for the simulation is " << diff / CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}
