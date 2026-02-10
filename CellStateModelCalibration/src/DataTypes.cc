/**
 * @file DataTypes.cc
 * @brief Implementation of data type methods
 */

#include "DataTypes.hh"

namespace CellStateCalibration {

void CellStateModelParams::print() const {
    std::cout << "=== Cell State Model Parameters ===\n";
    std::cout << std::fixed;
    std::cout << "  E1 (healthy state):     " << std::setprecision(4) << E1 << "\n";
    std::cout << "  E3 (death threshold):   " << std::setprecision(4) << E3 << "\n";
    std::cout << "  sigma (transition width): " << std::setprecision(4) << sigma << "\n";
    std::cout << "  alpha (damage per DSB): " << std::setprecision(4) << alpha << "\n";
    std::cout << "  kappa (DSB/Gy):         " << std::setprecision(4) << kappa << "\n";
    std::cout << "\n";
    std::cout << "Derived quantities:\n";
    std::cout << "  Critical DSB count (E3/alpha): " << std::setprecision(2) << criticalDSBCount() << " DSBs\n";
    std::cout << "  Critical dose (E3/alpha/kappa): " << std::setprecision(3) << criticalDose() << " Gy\n";
    std::cout << "===================================\n";
}

} // namespace CellStateCalibration
