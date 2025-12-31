//============================================================================
// Test Program: DiffusionReactionSolver Profiling
//============================================================================
//
// Purpose:
//   Profile the DiffusionReactionSolver to identify performance bottlenecks.
//
//============================================================================

#include <fstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <sys/stat.h>

#include "ReactionDiffusionSimulation.hh"
#include "DiffusionReactionSolver.hh"
#include "CellLayoutInitializer.hh"

using namespace std;
using namespace std::chrono;

// Timing helper class
class Timer {
public:
    void start() { start_time = high_resolution_clock::now(); }
    void stop() { 
        end_time = high_resolution_clock::now(); 
        total_time += duration_cast<microseconds>(end_time - start_time).count();
        call_count++;
    }
    double getTotalSeconds() const { return total_time / 1e6; }
    int getCallCount() const { return call_count; }
    double getAverageMs() const { return (call_count > 0) ? (total_time / 1000.0 / call_count) : 0; }
    void reset() { total_time = 0; call_count = 0; }
private:
    high_resolution_clock::time_point start_time, end_time;
    double total_time = 0;  // in microseconds
    int call_count = 0;
};

// Global timers
Timer timer_total;
Timer timer_diffusion_calc;
Timer timer_file_write;
Timer timer_get_concentration;

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
    cout << "  DiffusionReactionSolver PROFILING        " << endl;
    cout << "============================================" << endl;
    
    //------------------------------------------------------------------------
    // Setup simulation parameters
    //------------------------------------------------------------------------
    double xDim = 1;
    double yDim = 1;
    double zDim = 1;
    double d = 0.01;
    double D = 1E-5;
    double r = 0.63E-17;
    int Rt = 5000;
    double mu = 200;
    int cellNum = 100;
    double deltaT_diffusion = 1;
    double deltaT_diffusionUpdate = 1;
    int numTimeSteps = 100;
    
    cout << "\nSimulation Parameters:" << endl;
    cout << "  Domain: " << xDim << " x " << yDim << " x " << zDim << " mm" << endl;
    cout << "  Grid points: " << (int)(xDim/d) << " x " << (int)(yDim/d) << " x " << (int)(zDim/d) 
         << " = " << (int)(xDim/d) * (int)(yDim/d) * (int)(zDim/d) << " total" << endl;
    cout << "  Time steps: " << numTimeSteps << endl;
    cout << "  Number of cells: " << cellNum << endl;
    
    //------------------------------------------------------------------------
    // Initialize
    //------------------------------------------------------------------------
    DiffusionReactionSolver diffusionSolver;
    diffusionSolver.DiffusionReactionInitialization(xDim, yDim, zDim, d, D, r, Rt, deltaT_diffusion);
    diffusionSolver.SetVerbose(0);
    
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d, d, d);
    layout.RectangularSlab(xDim, yDim, zDim, cellNum);
    cout << "\nCell Layout Initialized: " << layout.GetCellNumber() << " cells" << endl;
    
    for (int i = 0; i < cellNum; i++)
    {
        double cX = layout.GetCellPositionX(i) + xDim / 2.0;
        double cY = layout.GetCellPositionY(i) + yDim / 2.0;
        double cZ = layout.GetCellPositionZ(i) + zDim / 2.0;
        int cellState = (i % 2 == 0) ? 1 : 2;
        diffusionSolver.CellStateUpdate(i, cX, cY, cZ, cellState, mu);
    }
    
    string outputPath = "./concentration_profile/";
    createDirectoryIfNotExists(outputPath);
    
    //------------------------------------------------------------------------
    // Run with profiling
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  PROFILING RUN" << endl;
    cout << "============================================" << endl;
    
    timer_total.start();
    
    for (int i = 0; i < numTimeSteps; i++)
    {
        // Profile: Diffusion calculation
        timer_diffusion_calc.start();
        diffusionSolver.DiffusionReactionCalculation(deltaT_diffusionUpdate, 1);
        timer_diffusion_calc.stop();
        
        // Profile: Get concentration query
        timer_get_concentration.start();
        double concentration = diffusionSolver.GetConcentration(0.5, 0.6, 0.5);
        timer_get_concentration.stop();
        
        // Profile: File writing
        timer_file_write.start();
        diffusionSolver.WriteConcentrationToFile(outputPath, i);
        timer_file_write.stop();
        
        if (i % 20 == 0) {
            cout << "  Step " << i << "/" << numTimeSteps << " completed" << endl;
        }
    }
    
    timer_total.stop();
    
    //------------------------------------------------------------------------
    // Print profiling results
    //------------------------------------------------------------------------
    cout << "\n============================================" << endl;
    cout << "  PROFILING RESULTS" << endl;
    cout << "============================================" << endl;
    
    double total_time = timer_total.getTotalSeconds();
    double diffusion_time = timer_diffusion_calc.getTotalSeconds();
    double file_time = timer_file_write.getTotalSeconds();
    double query_time = timer_get_concentration.getTotalSeconds();
    double other_time = total_time - diffusion_time - file_time - query_time;
    
    cout << "\nTime Breakdown:" << endl;
    cout << "  ┌─────────────────────────────────────────────────────────┐" << endl;
    printf("  │ %-30s %10.2f s  (%5.1f%%) │\n", "Diffusion Calculation:", diffusion_time, 100.0*diffusion_time/total_time);
    printf("  │ %-30s %10.2f s  (%5.1f%%) │\n", "File Writing:", file_time, 100.0*file_time/total_time);
    printf("  │ %-30s %10.2f s  (%5.1f%%) │\n", "Concentration Query:", query_time, 100.0*query_time/total_time);
    printf("  │ %-30s %10.2f s  (%5.1f%%) │\n", "Other (overhead):", other_time, 100.0*other_time/total_time);
    cout << "  ├─────────────────────────────────────────────────────────┤" << endl;
    printf("  │ %-30s %10.2f s  (100.0%%) │\n", "TOTAL:", total_time);
    cout << "  └─────────────────────────────────────────────────────────┘" << endl;
    
    cout << "\nPer-call Statistics:" << endl;
    printf("  Diffusion calc: %.2f ms/call (%d calls)\n", timer_diffusion_calc.getAverageMs(), timer_diffusion_calc.getCallCount());
    printf("  File write:     %.2f ms/call (%d calls)\n", timer_file_write.getAverageMs(), timer_file_write.getCallCount());
    printf("  Conc. query:    %.4f ms/call (%d calls)\n", timer_get_concentration.getAverageMs(), timer_get_concentration.getCallCount());
    
    cout << "\n============================================" << endl;
    cout << "  Profiling Complete!" << endl;
    cout << "============================================" << endl;

    return 0;
}

