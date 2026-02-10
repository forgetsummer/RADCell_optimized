/**
 * @file test_pide_reader.cc
 * @brief Test program for PIDE database reader
 * 
 * This program demonstrates how to use the PIDEDataReader class to:
 * 1. Load PIDE database
 * 2. Query experiments by cell line, ion, LET range
 * 3. Access raw survival data
 * 4. Compare simulation results with experimental data
 * 5. Fit LQ parameters to data
 */

#include "PIDEDataReader.hh"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace PIDE;

void printExperimentDetails(const ExperimentEntry& exp) {
    std::cout << "  ExpID: " << exp.expID 
              << ", Cell: " << exp.cellLine
              << ", Ion: " << exp.ion
              << ", LET: " << std::fixed << std::setprecision(2) << exp.LET_keV_um << " keV/μm"
              << ", Energy: " << exp.energy_MeV_u << " MeV/u"
              << std::endl;
    
    if (exp.ionLQ_paper.isValid) {
        std::cout << "    LQ (paper): α=" << std::setprecision(4) << exp.ionLQ_paper.alpha 
                  << " Gy⁻¹, β=" << exp.ionLQ_paper.beta << " Gy⁻²"
                  << ", D10=" << std::setprecision(2) << exp.ionLQ_paper.D10() << " Gy"
                  << std::endl;
    }
    if (exp.ionLQ_fit.isValid) {
        std::cout << "    LQ (fit):   α=" << std::setprecision(4) << exp.ionLQ_fit.alpha 
                  << " Gy⁻¹, β=" << exp.ionLQ_fit.beta << " Gy⁻²"
                  << ", D10=" << std::setprecision(2) << exp.ionLQ_fit.D10() << " Gy"
                  << std::endl;
    }
}

void printRawData(const RawSurvivalCurve& curve) {
    if (!curve.hasData) {
        std::cout << "  No raw data available" << std::endl;
        return;
    }
    
    std::cout << "  Raw data points (" << curve.dataPoints.size() << " points):" << std::endl;
    std::cout << "    Dose (Gy)   SF" << std::endl;
    for (const auto& pt : curve.dataPoints) {
        std::cout << "    " << std::fixed << std::setprecision(2) << std::setw(8) << pt.dose_Gy 
                  << "  " << std::setprecision(5) << pt.survivalFraction << std::endl;
    }
    
    auto range = curve.getDoseRange();
    std::cout << "  Dose range: " << range.first << " - " << range.second << " Gy" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "========================================" << std::endl;
    std::cout << "PIDE Database Reader Test" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Determine PIDE directory
    std::string pideDir = "../PIDE3.4";
    if (argc > 1) {
        pideDir = argv[1];
    }
    
    std::cout << "Loading PIDE data from: " << pideDir << std::endl;
    
    // Create reader and load data
    PIDEDataReader reader;
    if (!reader.loadData(pideDir)) {
        std::cerr << "Failed to load PIDE data!" << std::endl;
        std::cerr << "Make sure to run convert_pide_xlsx.py first to generate the CSV file." << std::endl;
        return 1;
    }
    
    // Print summary
    reader.printSummary();
    
    // Example 1: Query by cell line
    std::cout << "\n--- Example 1: Query V79 cell line experiments ---" << std::endl;
    auto v79Exps = reader.getExperimentsByCellLine("V79");
    std::cout << "Found " << v79Exps.size() << " experiments with V79 cells" << std::endl;
    if (!v79Exps.empty()) {
        std::cout << "First 3 experiments:" << std::endl;
        for (size_t i = 0; i < 3 && i < v79Exps.size(); ++i) {
            printExperimentDetails(v79Exps[i]);
        }
    }
    
    // Example 2: Query by ion type
    std::cout << "\n--- Example 2: Query Carbon-12 ion experiments ---" << std::endl;
    auto carbonExps = reader.getExperimentsByIon("12C");
    std::cout << "Found " << carbonExps.size() << " experiments with 12C ions" << std::endl;
    if (!carbonExps.empty()) {
        std::cout << "First 3 experiments:" << std::endl;
        for (size_t i = 0; i < 3 && i < carbonExps.size(); ++i) {
            printExperimentDetails(carbonExps[i]);
        }
    }
    
    // Example 3: Query by LET range (high LET: 50-200 keV/μm)
    std::cout << "\n--- Example 3: Query high-LET experiments (50-200 keV/μm) ---" << std::endl;
    auto highLETExps = reader.getExperimentsByLETRange(50.0, 200.0);
    std::cout << "Found " << highLETExps.size() << " experiments in LET range 50-200 keV/μm" << std::endl;
    if (!highLETExps.empty()) {
        std::cout << "First 3 experiments:" << std::endl;
        for (size_t i = 0; i < 3 && i < highLETExps.size(); ++i) {
            printExperimentDetails(highLETExps[i]);
        }
    }
    
    // Example 4: Get raw survival data
    std::cout << "\n--- Example 4: Raw survival data ---" << std::endl;
    if (!v79Exps.empty()) {
        int expID = v79Exps[0].expID;
        std::cout << "Ion raw data for experiment ID " << expID << ":" << std::endl;
        auto ionRaw = reader.getIonRawData(expID);
        printRawData(ionRaw);
        
        std::cout << "\nPhoton raw data for publication " << v79Exps[0].publicationID 
                  << ", exp " << v79Exps[0].photonExpNum << ":" << std::endl;
        auto photonRaw = reader.getPhotonRawData(v79Exps[0].publicationID, v79Exps[0].photonExpNum);
        printRawData(photonRaw);
    }
    
    // Example 5: Fit LQ parameters to simulation data
    std::cout << "\n--- Example 5: LQ Parameter Fitting ---" << std::endl;
    std::vector<double> simDoses = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<double> simSF = {1.0, 0.8, 0.5, 0.3, 0.15, 0.07, 0.03, 0.01, 0.003};
    
    auto fittedLQ = PIDEDataReader::fitLQParameters(simDoses, simSF);
    std::cout << "Simulated survival data:" << std::endl;
    for (size_t i = 0; i < simDoses.size(); ++i) {
        std::cout << "  D=" << simDoses[i] << " Gy, SF=" << simSF[i] << std::endl;
    }
    std::cout << "\nFitted LQ parameters:" << std::endl;
    std::cout << "  α = " << std::setprecision(4) << fittedLQ.alpha << " Gy⁻¹" << std::endl;
    std::cout << "  β = " << fittedLQ.beta << " Gy⁻²" << std::endl;
    std::cout << "  D10 = " << std::setprecision(2) << fittedLQ.D10() << " Gy" << std::endl;
    std::cout << "  D37 = " << fittedLQ.D37() << " Gy" << std::endl;
    std::cout << "  D50 = " << fittedLQ.D50() << " Gy" << std::endl;
    
    // Example 6: Query by cell type
    std::cout << "\n--- Example 6: Query tumor cell experiments ---" << std::endl;
    auto tumorExps = reader.getExperimentsByCellClass('t');
    std::cout << "Found " << tumorExps.size() << " tumor cell experiments" << std::endl;
    
    auto normalExps = reader.getExperimentsByCellClass('n');
    std::cout << "Found " << normalExps.size() << " normal cell experiments" << std::endl;
    
    // Example 7: Query by cell origin
    std::cout << "\n--- Example 7: Query by cell origin ---" << std::endl;
    auto humanExps = reader.getExperimentsByCellOrigin('h');
    std::cout << "Found " << humanExps.size() << " human cell experiments" << std::endl;
    
    auto rodentExps = reader.getExperimentsByCellOrigin('r');
    std::cout << "Found " << rodentExps.size() << " rodent cell experiments" << std::endl;
    
    // Example 8: List unique cell lines and ions
    std::cout << "\n--- Example 8: Available cell lines and ions ---" << std::endl;
    auto cellLines = reader.getUniqueCellLines();
    std::cout << "Cell lines (" << cellLines.size() << " total): ";
    for (size_t i = 0; i < 10 && i < cellLines.size(); ++i) {
        std::cout << cellLines[i];
        if (i < 9 && i < cellLines.size() - 1) std::cout << ", ";
    }
    if (cellLines.size() > 10) std::cout << ", ...";
    std::cout << std::endl;
    
    auto ions = reader.getUniqueIons();
    std::cout << "Ions (" << ions.size() << " total): ";
    for (size_t i = 0; i < ions.size(); ++i) {
        std::cout << ions[i];
        if (i < ions.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "PIDE Database Reader Test Complete!" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
