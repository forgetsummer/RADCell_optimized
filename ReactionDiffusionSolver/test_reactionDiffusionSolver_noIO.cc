//============================================================================
// Test Program: DiffusionReactionSolver - Minimal I/O Version
//============================================================================
//
// Purpose:
//   Test the DiffusionReactionSolver performance with minimal file I/O.
//   Only writes the final time step results to file.
//
// This test demonstrates the true computational performance of the
// diffusion-reaction solver without the file I/O bottleneck.
//
// Author: Ruirui Liu
//         Department of Nuclear Engineering and Radiation Health Physics
//         Oregon State University
//
// Usage:
//   Compile with: make test_reactionDiffusionSolver_noIO
//   Run with: ./test_reactionDiffusionSolver_noIO
//
//============================================================================

#include <fstream>
#include <iostream>
#include <ctime>
#include <sys/stat.h>

#include "ReactionDiffusionSimulation.hh"
#include "DiffusionReactionSolver.hh"
#include "CellLayoutInitializer.hh"

using namespace std;

bool createDirectoryIfNotExists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        if (mkdir(path.c_str(), 0755) == 0) {
            std::cout << "Created output directory: " << path << std::endl;
            return true;
        } else {
            std::cerr << "Error: Failed to create directory: " << path << std::endl;
            return false;
        }
    } else if (info.st_mode & S_IFDIR) {
        return true;
    } else {
        std::cerr << "Error: Path exists but is not a directory: " << path << std::endl;
        return false;
    }
}

int main(int argc, char** argv)
{
    cout << "============================================" << endl;
    cout << "  DiffusionReactionSolver - Minimal I/O    " << endl;
    cout << "============================================" << endl;
    
    //------------------------------------------------------------------------
    // Setup simulation parameters
    //------------------------------------------------------------------------
    double xDim = 1;            // Domain size in x direction (mm)
    double yDim = 1;            // Domain size in y direction (mm)
    double zDim = 1;            // Domain size in z direction (mm)
    double d = 0.01;            // Mesh spacing (mm)
    double D = 1E-5;            // Diffusion coefficient (mm^2/s)
    double r = 0.63E-17;        // Reaction rate constant (1/#.s)
    int Rt = 5000;              // Receptor threshold (#)
    double mu = 200;            // Secretion rate (#/s)
    int cellNum = 100;          // Number of cells
    double deltaT_diffusion = 1;       // Time step for diffusion (s)
    double deltaT_diffusionUpdate = 1; // Time step for update (s)
    int numTimeSteps = 100;            // Number of time steps to simulate
    
    cout << "\nSimulation Parameters:" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "  Grid points: " << (int)(xDim/d) << " x " << (int)(yDim/d) << " x " << (int)(zDim/d) 
         << " = " << (int)(xDim/d) * (int)(yDim/d) * (int)(zDim/d) << " total" << endl;
    cout << "  Mesh spacing: " << d << " mm" << endl;
    cout << "  Diffusion coefficient: " << D << " mm^2/s" << endl;
    cout << "  Time steps: " << numTimeSteps << endl;
    cout << "  Number of cells: " << cellNum << endl;
    
    //------------------------------------------------------------------------
    // Initialize diffusion-reaction solver
    //------------------------------------------------------------------------
    DiffusionReactionSolver diffusionSolver;
    diffusionSolver.DiffusionReactionInitialization(xDim, yDim, zDim, d, D, r, Rt, deltaT_diffusion);
    diffusionSolver.SetVerbose(0);
    
    //------------------------------------------------------------------------
    // Initialize cell layout
    //------------------------------------------------------------------------
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d, d, d);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);
    cout << "\nCell Layout Initialized:" << endl;
    cout << "  Total cells: " << layout.GetCellNumber() << endl;
    
    //------------------------------------------------------------------------
    // Setup cell positions and states
    //------------------------------------------------------------------------
    for (int i = 0; i < cellNum; i++)
    {
        double cX = layout.GetCellPositionX(i) + xDim / 2.0;
        double cY = layout.GetCellPositionY(i) + yDim / 2.0;
        double cZ = layout.GetCellPositionZ(i) + zDim / 2.0;
        
        // Alternate cell states: 1 (healthy) and 2 (signaling/arrest)
        int cellState = (i % 2 == 0) ? 1 : 2;
        diffusionSolver.CellStateUpdate(i, cX, cY, cZ, cellState, mu);
    }
    
    //------------------------------------------------------------------------
    // Run step-by-step diffusion-reaction simulation (NO FILE I/O)
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  Running simulation (NO intermediate I/O)" << endl;
    cout << "============================================" << endl;
    cout << "Monitoring concentration at point (0.5, 0.6, 0.5):" << endl;
    
    clock_t t1 = clock();
    
    for (int i = 0; i < numTimeSteps; i++)
    {
        diffusionSolver.DiffusionReactionCalculation(deltaT_diffusionUpdate, 1);
        
        // Query concentration at a sample point (very fast, negligible time)
        double concentration = diffusionSolver.GetConcentration(0.5, 0.6, 0.5);
        
        // Print progress every 20 steps
        if (i % 20 == 0) {
            cout << "  Step " << i << "/" << numTimeSteps << ": concentration = " << concentration << endl;
        }
        
        // NO FILE WRITING during simulation loop!
    }
    
    clock_t t2 = clock();
    double computation_time = (double)(t2 - t1) / CLOCKS_PER_SEC;
    
    //------------------------------------------------------------------------
    // Write ONLY final results
    //------------------------------------------------------------------------
    string outputPath = "./concentration_final/";
    if (!createDirectoryIfNotExists(outputPath)) {
        cerr << "Failed to create output directory. Exiting." << endl;
        return 1;
    }
    
    cout << "\nWriting FINAL concentration profile to: " << outputPath << endl;
    
    clock_t t3 = clock();
    diffusionSolver.WriteConcentrationToFile(outputPath, numTimeSteps - 1);
    clock_t t4 = clock();
    double io_time = (double)(t4 - t3) / CLOCKS_PER_SEC;
    
    //------------------------------------------------------------------------
    // Summary
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  PERFORMANCE SUMMARY" << endl;
    cout << "============================================" << endl;
    cout << "\nTime Breakdown:" << endl;
    cout << "  ┌─────────────────────────────────────────────────────────┐" << endl;
    printf("  │ %-35s %10.2f s         │\n", "Computation (diffusion):", computation_time);
    printf("  │ %-35s %10.2f s         │\n", "File I/O (final step only):", io_time);
    cout << "  ├─────────────────────────────────────────────────────────┤" << endl;
    printf("  │ %-35s %10.2f s         │\n", "TOTAL:", computation_time + io_time);
    cout << "  └─────────────────────────────────────────────────────────┘" << endl;
    
    cout << "\nComparison with full I/O version:" << endl;
    cout << "  Full I/O version:    ~430 seconds (100 file writes)" << endl;
    printf("  This version:        %.2f seconds (1 file write)\n", computation_time + io_time);
    printf("  Speedup:             %.1fx faster\n", 430.0 / (computation_time + io_time));
    
    cout << "\n============================================" << endl;
    cout << "  Test completed successfully!              " << endl;
    cout << "============================================" << endl;

    return 0;
}

