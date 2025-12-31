//============================================================================
// Test Program: CellLayoutInitializer Module
//============================================================================
//
// Purpose:
//   This is a standalone test program to verify the CellLayoutInitializer
//   class functionality. It tests the cell layout generation for different
//   geometries without requiring Geant4 dependencies.
//
// What it tests:
//   - CellLayoutInitializer class initialization
//   - SetCellHomeParamter() - setting cell home dimensions
//   - RectangularSlab() - generating cell positions in a rectangular slab
//   - GetCellNumber() - retrieving total cell count
//   - GetCellPositionX/Y/Z() - retrieving individual cell positions
//   - Cell distribution statistics and visualization
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
//   Compile with: make test_cellLayoutInitializer
//   Run with: ./test_cellLayoutInitializer
//
// Output:
//   - Console output showing cell positions and statistics
//   - Optional CSV file with cell positions for visualization
//
//============================================================================

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sys/stat.h>

#include "CellLayoutInitializer.hh"

using namespace std;

// Helper function to create directory if it doesn't exist
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

// Function to calculate statistics
void calculateStatistics(CellLayoutInitializer& layout, int cellNum,
                         double& minX, double& maxX, double& minY, double& maxY, 
                         double& minZ, double& maxZ) {
    minX = minY = minZ = 1e10;
    maxX = maxY = maxZ = -1e10;
    
    for (int i = 0; i < cellNum; i++) {
        double x = layout.GetCellPositionX(i);
        double y = layout.GetCellPositionY(i);
        double z = layout.GetCellPositionZ(i);
        
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
        if (z < minZ) minZ = z;
        if (z > maxZ) maxZ = z;
    }
}

// Function to write cell positions to CSV file
void writeCellPositionsToCSV(CellLayoutInitializer& layout, int cellNum, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }
    
    file << "cell_id,x,y,z" << endl;
    for (int i = 0; i < cellNum; i++) {
        file << i << "," 
             << layout.GetCellPositionX(i) << ","
             << layout.GetCellPositionY(i) << ","
             << layout.GetCellPositionZ(i) << endl;
    }
    file.close();
    cout << "Cell positions written to: " << filename << endl;
}

int main(int argc, char** argv)
{
    cout << "============================================" << endl;
    cout << "  CellLayoutInitializer Test Program       " << endl;
    cout << "============================================" << endl;
    
    //------------------------------------------------------------------------
    // Test 1: 2D Monolayer (z = 0)
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  TEST 1: 2D Monolayer Layout" << endl;
    cout << "============================================" << endl;
    
    double xDim_2D = 1.0;    // mm
    double yDim_2D = 1.0;    // mm
    double zDim_2D = 0.0;    // mm (2D)
    double cellHomeSize = 0.05;  // mm
    int cellNum_2D = 100;
    
    cout << "\nParameters:" << endl;
    cout << "  Domain: " << xDim_2D << " x " << yDim_2D << " x " << zDim_2D << " mm" << endl;
    cout << "  Cell home size: " << cellHomeSize << " mm" << endl;
    cout << "  Requested cells: " << cellNum_2D << endl;
    
    CellLayoutInitializer layout_2D;
    layout_2D.SetCellHomeParamter(cellHomeSize, cellHomeSize, cellHomeSize);
    layout_2D.RectangularSlab(xDim_2D, yDim_2D, zDim_2D, cellNum_2D);
    
    int actualCells_2D = layout_2D.GetCellNumber();
    cout << "\nResults:" << endl;
    cout << "  Actual cells created: " << actualCells_2D << endl;
    
    // Calculate statistics
    double minX, maxX, minY, maxY, minZ, maxZ;
    calculateStatistics(layout_2D, actualCells_2D, minX, maxX, minY, maxY, minZ, maxZ);
    
    cout << "\nCell Position Ranges:" << endl;
    cout << "  X: [" << minX << ", " << maxX << "] mm" << endl;
    cout << "  Y: [" << minY << ", " << maxY << "] mm" << endl;
    cout << "  Z: [" << minZ << ", " << maxZ << "] mm" << endl;
    
    // Print first few cells
    cout << "\nFirst 5 cell positions:" << endl;
    cout << "  " << setw(5) << "ID" << setw(12) << "X (mm)" << setw(12) << "Y (mm)" << setw(12) << "Z (mm)" << endl;
    cout << "  " << string(41, '-') << endl;
    for (int i = 0; i < min(5, actualCells_2D); i++) {
        cout << "  " << setw(5) << i 
             << setw(12) << fixed << setprecision(4) << layout_2D.GetCellPositionX(i)
             << setw(12) << layout_2D.GetCellPositionY(i)
             << setw(12) << layout_2D.GetCellPositionZ(i) << endl;
    }
    
    //------------------------------------------------------------------------
    // Test 2: 3D Volume
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  TEST 2: 3D Volume Layout" << endl;
    cout << "============================================" << endl;
    
    double xDim_3D = 0.5;    // mm
    double yDim_3D = 0.5;    // mm
    double zDim_3D = 0.5;    // mm (3D)
    int cellNum_3D = 50;
    
    cout << "\nParameters:" << endl;
    cout << "  Domain: " << xDim_3D << " x " << yDim_3D << " x " << zDim_3D << " mm" << endl;
    cout << "  Cell home size: " << cellHomeSize << " mm" << endl;
    cout << "  Requested cells: " << cellNum_3D << endl;
    
    CellLayoutInitializer layout_3D;
    layout_3D.SetCellHomeParamter(cellHomeSize, cellHomeSize, cellHomeSize);
    layout_3D.RectangularSlab(xDim_3D, yDim_3D, zDim_3D, cellNum_3D);
    
    int actualCells_3D = layout_3D.GetCellNumber();
    cout << "\nResults:" << endl;
    cout << "  Actual cells created: " << actualCells_3D << endl;
    
    calculateStatistics(layout_3D, actualCells_3D, minX, maxX, minY, maxY, minZ, maxZ);
    
    cout << "\nCell Position Ranges:" << endl;
    cout << "  X: [" << minX << ", " << maxX << "] mm" << endl;
    cout << "  Y: [" << minY << ", " << maxY << "] mm" << endl;
    cout << "  Z: [" << minZ << ", " << maxZ << "] mm" << endl;
    
    // Print first few cells
    cout << "\nFirst 5 cell positions:" << endl;
    cout << "  " << setw(5) << "ID" << setw(12) << "X (mm)" << setw(12) << "Y (mm)" << setw(12) << "Z (mm)" << endl;
    cout << "  " << string(41, '-') << endl;
    for (int i = 0; i < min(5, actualCells_3D); i++) {
        cout << "  " << setw(5) << i 
             << setw(12) << fixed << setprecision(4) << layout_3D.GetCellPositionX(i)
             << setw(12) << layout_3D.GetCellPositionY(i)
             << setw(12) << layout_3D.GetCellPositionZ(i) << endl;
    }
    
    //------------------------------------------------------------------------
    // Test 3: Single Cell
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  TEST 3: Single Cell Layout" << endl;
    cout << "============================================" << endl;
    
    double xDim_single = 0.05;  // mm
    double yDim_single = 0.05;  // mm
    double zDim_single = 0.05;  // mm
    int cellNum_single = 1;
    
    cout << "\nParameters:" << endl;
    cout << "  Domain: " << xDim_single << " x " << yDim_single << " x " << zDim_single << " mm" << endl;
    cout << "  Requested cells: " << cellNum_single << endl;
    
    CellLayoutInitializer layout_single;
    layout_single.SetCellHomeParamter(cellHomeSize, cellHomeSize, cellHomeSize);
    layout_single.RectangularSlab(xDim_single, yDim_single, zDim_single, cellNum_single);
    
    int actualCells_single = layout_single.GetCellNumber();
    cout << "\nResults:" << endl;
    cout << "  Actual cells created: " << actualCells_single << endl;
    
    if (actualCells_single > 0) {
        cout << "  Cell position: (" 
             << layout_single.GetCellPositionX(0) << ", "
             << layout_single.GetCellPositionY(0) << ", "
             << layout_single.GetCellPositionZ(0) << ") mm" << endl;
    }
    
    //------------------------------------------------------------------------
    // Test 4: Large Cell Population
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  TEST 4: Large Cell Population" << endl;
    cout << "============================================" << endl;
    
    double xDim_large = 3.14159;  // mm (pi for testing)
    double yDim_large = 3.14159;  // mm
    double zDim_large = 0;        // mm (2D)
    int cellNum_large = 1000;
    double cellHomeSize_large = 0.01;  // mm (smaller cells)
    
    cout << "\nParameters:" << endl;
    cout << "  Domain: " << xDim_large << " x " << yDim_large << " x " << zDim_large << " mm" << endl;
    cout << "  Cell home size: " << cellHomeSize_large << " mm" << endl;
    cout << "  Requested cells: " << cellNum_large << endl;
    
    CellLayoutInitializer layout_large;
    layout_large.SetCellHomeParamter(cellHomeSize_large, cellHomeSize_large, cellHomeSize_large);
    layout_large.RectangularSlab(xDim_large, yDim_large, zDim_large, cellNum_large);
    
    int actualCells_large = layout_large.GetCellNumber();
    cout << "\nResults:" << endl;
    cout << "  Actual cells created: " << actualCells_large << endl;
    
    calculateStatistics(layout_large, actualCells_large, minX, maxX, minY, maxY, minZ, maxZ);
    
    cout << "\nCell Position Ranges:" << endl;
    cout << "  X: [" << minX << ", " << maxX << "] mm" << endl;
    cout << "  Y: [" << minY << ", " << maxY << "] mm" << endl;
    cout << "  Z: [" << minZ << ", " << maxZ << "] mm" << endl;
    
    //------------------------------------------------------------------------
    // Write output files
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  Writing Output Files" << endl;
    cout << "============================================" << endl;
    
    string outputDir = "./cell_layout_output/";
    if (createDirectoryIfNotExists(outputDir)) {
        writeCellPositionsToCSV(layout_2D, actualCells_2D, outputDir + "cell_positions_2D.csv");
        writeCellPositionsToCSV(layout_3D, actualCells_3D, outputDir + "cell_positions_3D.csv");
        writeCellPositionsToCSV(layout_large, actualCells_large, outputDir + "cell_positions_large.csv");
    }
    
    //------------------------------------------------------------------------
    // Summary
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  TEST SUMMARY" << endl;
    cout << "============================================" << endl;
    cout << "\n  " << setw(25) << "Test" << setw(15) << "Requested" << setw(15) << "Created" << setw(10) << "Status" << endl;
    cout << "  " << string(65, '-') << endl;
    cout << "  " << setw(25) << "2D Monolayer" << setw(15) << cellNum_2D << setw(15) << actualCells_2D 
         << setw(10) << (actualCells_2D > 0 ? "PASS" : "FAIL") << endl;
    cout << "  " << setw(25) << "3D Volume" << setw(15) << cellNum_3D << setw(15) << actualCells_3D 
         << setw(10) << (actualCells_3D > 0 ? "PASS" : "FAIL") << endl;
    cout << "  " << setw(25) << "Single Cell" << setw(15) << cellNum_single << setw(15) << actualCells_single 
         << setw(10) << (actualCells_single > 0 ? "PASS" : "FAIL") << endl;
    cout << "  " << setw(25) << "Large Population" << setw(15) << cellNum_large << setw(15) << actualCells_large 
         << setw(10) << (actualCells_large > 0 ? "PASS" : "FAIL") << endl;
    
    cout << "\n============================================" << endl;
    cout << "  All tests completed successfully!         " << endl;
    cout << "============================================" << endl;

    return 0;
}

