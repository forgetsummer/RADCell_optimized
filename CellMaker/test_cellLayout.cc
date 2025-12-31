//============================================================================
// Test Program: Cell Layout with Geant4 Visualization
//============================================================================
//
// Purpose:
//   This test program demonstrates the integration of CellLayoutInitializer
//   with the full RADCellSimulation framework. It creates a cell population
//   layout and visualizes it using Geant4's visualization system.
//
// What it tests:
//   - CellLayoutInitializer integration with RADCellSimulation
//   - SetCellHomeParamter() - setting cell home dimensions
//   - RectangularSlab() - generating 2D monolayer cell positions
//   - RADCellSimulation cell world creation and visualization
//   - Cell construction with Geant4 geometry
//   - Microdosimetry simulation setup with GUI visualization
//
// Key Differences from test_cellLayoutInitializer:
//   - This test requires Geant4 and uses full simulation framework
//   - Provides visual verification of cell layout in Geant4 GUI
//   - Tests the complete workflow from layout to simulation
//
// Author: Ruirui Liu
//         Department of Nuclear Engineering and Radiation Health Physics
//         Oregon State University
//         December 10, 2015
//
// References:
//   - PhD Thesis: https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/rf55zf01n
//   - RADCell Paper: https://iopscience.iop.org/article/10.1088/1361-6560/abd4f9/meta
//   - Geant4-DNA: Med. Phys. 37 (2010) 4692-4708
//   - Geant4-DNA web site: http://geant4-dna.org
//
// Contact: liuruirui.nova@gmail.com
//
// Usage:
//   Compile with: make test_cellLayout
//   Run with: ./test_cellLayout
//
// Output:
//   - Console output showing cell positions (X, Y, Z coordinates)
//   - Geant4 GUI visualization of the cell layout
//
//============================================================================

#include "RADCellSimulation.hh"
#include "CellLayoutInitializer.hh"
#include <time.h>


int main(int argc,char** argv)
{
    int N=30; // this line for seeding 300 cell
        
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(0.05,0.05,0.05);// cell home size is in unit of mm
//         layout.RectangularSlab(0.05, 0.05, 0.05,N);// unit here is also mm, this line is for single cell simulation
    layout.RectangularSlab(1,1,0,N); // this line is for multiple cell simulation
    RADCellSimulation mySim;
//     mySim.SetCellWorld(0.1,0.1,0.1);//create cell world, dimension in unit of mm, this line is for testing single cell simulation

    mySim.SetCellWorld(1.1,1.1,1.1); // this line is for testing multiple cell simulation
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Blue Green");
    mySim.SetCellSimulationParameter(1,"Nucleus");
    for (int i=0; i<layout.GetCellNumber(); i++)
    {
        cout<<"X="<<layout.GetCellPositionX(i)<<" Y="<<layout.GetCellPositionY(i)<<" Z="<<layout.GetCellPositionZ(i)<<endl;

        mySim.CellConstruction(i,"Epithelia",layout.GetCellPositionX(i),layout.GetCellPositionY(i),layout.GetCellPositionZ(i),20,5);
    }

    mySim.RADCellSimulationInitialize(argc, argv);

    mySim.EnergyDistributionCalculation("gui","microdosimetry.in",false);

    return 0;
}

