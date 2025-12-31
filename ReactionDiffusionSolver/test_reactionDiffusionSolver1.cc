//============================================================================
// Test Program: DiffusionReactionSolver Step-by-Step Simulation
//============================================================================
//
// Purpose:
//   This is a standalone test program to verify the DiffusionReactionSolver
//   class, which provides a higher-level interface for running reaction-
//   diffusion simulations in a step-by-step manner. This is useful for
//   coupling with other simulation components (e.g., radiation transport).
//
// What it tests:
//   - DiffusionReactionSolver class initialization
//   - Step-by-step diffusion-reaction calculation
//   - Real-time concentration queries at specific points
//   - Incremental concentration profile output
//   - CellLayoutInitializer for cell positioning
//
// Key difference from test_reactionDiffusionKernel:
//   - Uses DiffusionReactionSolver (high-level, step-by-step interface)
//   - Allows querying concentration at arbitrary points during simulation
//   - Writes output at each time step (for animation/visualization)
//
// Author: Ruirui Liu
//         Department of Nuclear Engineering and Radiation Health Physics
//         Oregon State University
//         December 10, 2015
//
// References:
//   - PhD Thesis: https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/rf55zf01n
//   - RADCell Paper: https://iopscience.iop.org/article/10.1088/1361-6560/abd4f9/meta

//
// Contact: liuruirui.nova@gmail.com
//
// Usage:
//   Compile with: make test_reactionDiffusionSolver1
//   Run with: ./test_reactionDiffusionSolver1
//
// Output:
//   Concentration profile files written to ./concentration/ directory
//   Console output showing concentration at a sample point over time
//
//============================================================================

#include <fstream>
#include <iostream>
#include <ctime>
#include <sys/stat.h>  // For mkdir and stat

#include "ReactionDiffusionSimulation.hh"
#include "DiffusionReactionSolver.hh"
#include "CellLayoutInitializer.hh"

using namespace std;

// Helper function to create directory if it doesn't exist
bool createDirectoryIfNotExists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        // Directory doesn't exist, create it
        if (mkdir(path.c_str(), 0755) == 0) {
            std::cout << "Created output directory: " << path << std::endl;
            return true;
        } else {
            std::cerr << "Error: Failed to create directory: " << path << std::endl;
            return false;
        }
    } else if (info.st_mode & S_IFDIR) {
        // Directory already exists
        return true;
    } else {
        // Path exists but is not a directory
        std::cerr << "Error: Path exists but is not a directory: " << path << std::endl;
        return false;
    }
}

int main(int argc, char** argv)
{
    cout << "============================================" << endl;
    cout << "  DiffusionReactionSolver Test Program     " << endl;
    cout << "============================================" << endl;
    
    //------------------------------------------------------------------------
    // Setup simulation parameters
    //------------------------------------------------------------------------
    double xDim = 1;            // Domain size in x direction (mm)
    double yDim = 1;            // Domain size in y direction (mm)
    double zDim = 1;            // Domain size in z direction (mm)
    double d = 0.01;            // Mesh spacing (mm)
    double D = 1E-5;            // Diffusion coefficient (mm^2/s)
    double r = 2E-7;        // Reaction rate constant (1/#.s)
    int Rt = 5000;              // Receptor threshold (#)
    double mu = 200;            // Secretion rate (#/s)
    int cellNum = 100;          // Number of cells
    double deltaT_diffusion = 1;       // Time step for diffusion (s)
    double deltaT_diffusionUpdate = 1; // Time step for update (s)
    int numTimeSteps = 100;           // Number of time steps to simulate
    
    cout << "\nSimulation Parameters:" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
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
    // Setup output directory
    //------------------------------------------------------------------------
    string outputPath = "./concentration/";
    if (!createDirectoryIfNotExists(outputPath)) {
        cerr << "Failed to create output directory. Exiting." << endl;
        return 1;
    }
    
    //------------------------------------------------------------------------
    // Run step-by-step diffusion-reaction simulation
    //------------------------------------------------------------------------
    cout << "\nRunning step-by-step diffusion-reaction simulation..." << endl;
    cout << "Monitoring concentration at point (0.5, 0.6, 0.5):" << endl;
    
    clock_t t1 = clock();
    
    for (int i = 0; i < numTimeSteps; i++)
    {
        diffusionSolver.DiffusionReactionCalculation(deltaT_diffusionUpdate, 1);
        
        // Query concentration at a sample point
        double concentration = diffusionSolver.GetConcentration(0.5, 0.6, 0.5);
        
        // Print progress every 100 steps
        if (i % 100 == 0) {
            cout << "  Step " << i << ": concentration = " << concentration << endl;
        }
        
        // Write concentration profile to file
        diffusionSolver.WriteConcentrationToFile(outputPath, i);
    }
    
    clock_t t2 = clock();
    double elapsed = (double)(t2 - t1) / CLOCKS_PER_SEC;
    
    //------------------------------------------------------------------------
    // Summary
    //------------------------------------------------------------------------
    cout << "\nSimulation completed in " << elapsed << " seconds" << endl;
    cout << "Output files written to: " << outputPath << endl;
    
    cout << "\n============================================" << endl;
    cout << "  Test completed successfully!              " << endl;
    cout << "============================================" << endl;

    return 0;
}

