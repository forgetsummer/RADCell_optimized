//============================================================================
// Test Program: Reaction-Diffusion Simulation Kernel
//============================================================================
//
// Purpose:
//   This is a standalone test program to verify the reaction-diffusion 
//   simulation component of the RADCell project. It tests the bystander 
//   signal diffusion and reaction kinetics without running the full 
//   Geant4 radiation transport simulation.
//
// What it tests:
//   - ReactionDiffusionSimulation class
//   - CellLayoutInitializer class  
//   - Mesh geometry setup
//   - Cell state initialization with secretion rates
//   - Diffusion-reaction calculation over time
//   - Concentration profile output
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
//   Compile with: make test_reactionDiffusion
//   Run with: ./test_reactionDiffusion
//
// Output:
//   Concentration profile files written to the specified output directory
//
//============================================================================

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <sstream>
#include <sys/stat.h>  // For mkdir and stat

#include "CellLayoutInitializer.hh"
#include "ReactionDiffusionSimulation.hh"

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
        std::cout << "Output directory already exists: " << path << std::endl;
        return true;
    } else {
        // Path exists but is not a directory
        std::cerr << "Error: Path exists but is not a directory: " << path << std::endl;
        return false;
    }
}

using namespace std;

int main(int argc, char** argv)
{
    cout << "============================================" << endl;
    cout << "  Reaction-Diffusion Kernel Test Program   " << endl;
    cout << "============================================" << endl;
    
    //------------------------------------------------------------------------
    // Setup simulation parameters
    //------------------------------------------------------------------------
    double xDim = M_PI;     // Domain size in x direction (mm)
    double yDim = M_PI;     // Domain size in y direction (mm)
    double zDim = 0;        // Domain size in z direction (mm), 0 for 2D simulation
    double d = 0.01;        // Mesh spacing (mm)
    double D = 1E-6;        // Diffusion coefficient (mm^2/s)
    double r = 0.63E-17;    // Reaction rate constant (1/#.s)
    int Rt = 5000;          // Threshold parameter
    double mu = 200;        // Secretion rate (#/s)
    double deltaT = 1;      // Time step (s)
    double T = 100;         // Total simulation time (s)
    int cellNum = 100;      // Number of cells
    
    cout << "\nSimulation Parameters:" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "  Mesh spacing: " << d << " mm" << endl;
    cout << "  Diffusion coefficient: " << D << " mm^2/s" << endl;
    cout << "  Time step: " << deltaT << " s" << endl;
    cout << "  Total time: " << T << " s" << endl;
    cout << "  Number of cells: " << cellNum << endl;
    
    //------------------------------------------------------------------------
    // Initialize reaction-diffusion simulation
    //------------------------------------------------------------------------
    ReactionDiffusionSimulation react;
    react.MeshGeometry(xDim, yDim, zDim, d);
    
    //------------------------------------------------------------------------
    // Initialize cell layout
    //------------------------------------------------------------------------
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d, d, d);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);
    
    cout << "\nCell Layout Initialized:" << endl;
    cout << "  Total cells: " << cellNum << endl;
    
    //------------------------------------------------------------------------
    // Setup cell positions and states
    //------------------------------------------------------------------------
    std::map<int, double> cX;
    std::map<int, double> cY;
    std::map<int, double> cZ;
    std::map<int, int> cellState;
    std::map<int, double> secretionRate;

    for (int i = 0; i < cellNum; i++)
    {
        cX[i] = layout.GetCellPositionX(i) + xDim / 2.0;
        cY[i] = layout.GetCellPositionY(i) + yDim / 2.0;
        cZ[i] = layout.GetCellPositionZ(i) + zDim / 2.0;
        secretionRate[i] = mu;
        
        // Alternate cell states: 1 (healthy) and 2 (signaling)
        if (i % 2 == 0)
        {
            cellState[i] = 1;
        }
        else
        {
            cellState[i] = 2;
        }
    }
    
    react.CellStateInitialization(cX, cY, cZ, cellState, secretionRate);
    
    //------------------------------------------------------------------------
    // Setup initial concentration field
    //------------------------------------------------------------------------
    std::vector<std::vector<std::vector<double>>> initialC;
    
    int N_X = react.GetMeshDIMX();
    int N_Y = react.GetMeshDIMY();
    int N_Z = react.GetMeshDIMZ();
    
    cout << "\nMesh Dimensions:" << endl;
    cout << "  N_X = " << N_X << ", N_Y = " << N_Y << ", N_Z = " << N_Z << endl;

    for (int i = 0; i < N_X; i++)
    {
        initialC.push_back(std::vector<std::vector<double>>());
        for (int j = 0; j < N_Y; j++)
        {
            initialC[i].push_back(std::vector<double>());
            for (int k = 0; k < N_Z; k++)
            {
                initialC[i][j].push_back(0); // Initial concentration = 0
            }
        }
    }

    //------------------------------------------------------------------------
    // Run diffusion-reaction simulation
    //------------------------------------------------------------------------
    cout << "\nRunning diffusion-reaction simulation..." << endl;
    clock_t t1 = clock();
    
    react.DiffusionReactionCalculation(initialC, D, r, Rt, deltaT, T);
    
    clock_t t2 = clock();
    double elapsed = (double)(t2 - t1) / CLOCKS_PER_SEC;
    cout << "Simulation completed in " << elapsed << " seconds" << endl;
    
    //------------------------------------------------------------------------
    // Write output
    //------------------------------------------------------------------------
    string outputPath = "./concentration/";
    
    // Create output directory if it doesn't exist
    if (!createDirectoryIfNotExists(outputPath)) {
        cerr << "Failed to create output directory. Exiting." << endl;
        return 1;
    }
    
    cout << "\nWriting concentration profile to: " << outputPath << endl;
    react.WriteConcentrationProfile(outputPath);
    
    cout << "\n============================================" << endl;
    cout << "  Test completed successfully!              " << endl;
    cout << "============================================" << endl;

    return 0;
}
