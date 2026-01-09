/**
 * @file PIDEDataReader.cc
 * @brief Implementation of PIDE database reader
 */

#include "PIDEDataReader.hh"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>
#include <iomanip>

namespace PIDE {

PIDEDataReader::PIDEDataReader() {}

bool PIDEDataReader::loadData(const std::string& pideDirectory) {
    std::string csvFile = pideDirectory + "/PIDE3.4.csv";
    std::string photonDatFile = pideDirectory + "/PIDE3.4_PhotonRawData.dat";
    std::string ionDatFile = pideDirectory + "/PIDE3.4_IonRawData.dat";
    
    bool success = true;
    
    // Try to load main database
    if (!loadMainDatabase(csvFile)) {
        std::cerr << "Warning: Could not load main database from " << csvFile << std::endl;
        std::cerr << "Please run convert_pide_xlsx.py to generate the CSV file." << std::endl;
        success = false;
    }
    
    // Load raw data files
    if (!loadPhotonRawData(photonDatFile)) {
        std::cerr << "Warning: Could not load photon raw data from " << photonDatFile << std::endl;
    }
    
    if (!loadIonRawData(ionDatFile)) {
        std::cerr << "Warning: Could not load ion raw data from " << ionDatFile << std::endl;
    }
    
    dataLoaded_ = success && !experiments_.empty();
    return dataLoaded_;
}

bool PIDEDataReader::loadMainDatabase(const std::string& csvFile) {
    std::ifstream file(csvFile);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    // Skip header line
    if (!std::getline(file, line)) {
        return false;
    }
    
    experiments_.clear();
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::vector<std::string> fields = splitLine(line, ',');
        if (fields.size() < 24) continue;
        
        ExperimentEntry entry;
        
        // Parse identifiers
        entry.expID = parseInt(fields[0]);
        entry.publicationID = parseInt(fields[1]);
        entry.publicationName = trim(fields[2]);
        entry.ionExpNum = parseInt(fields[3]);
        
        // Parse cell information
        entry.cellLine = trim(fields[4]);
        entry.cellClass = fields[5].empty() ? 'n' : fields[5][0];
        entry.cellOrigin = fields[6].empty() ? 'h' : fields[6][0];
        entry.cellCycle = trim(fields[7]);
        entry.dnaContent = parseDouble(fields[8]);
        
        // Parse radiation information
        entry.photonRadiation = trim(fields[9]);
        entry.photonExpNum = parseInt(fields[10]);
        entry.ion = trim(fields[11]);
        entry.charge = parseInt(fields[12]);
        entry.irradiationConditions = fields[13].empty() ? 'm' : fields[13][0];
        entry.LET_keV_um = parseDouble(fields[14]);
        entry.energy_MeV_u = parseDouble(fields[15]);
        
        // Parse LQ parameters from publication
        entry.photonLQ_paper.alpha = parseDouble(fields[16]);
        entry.photonLQ_paper.beta = parseDouble(fields[17]);
        entry.photonLQ_paper.isValid = !std::isnan(entry.photonLQ_paper.alpha);
        
        entry.ionLQ_paper.alpha = parseDouble(fields[18]);
        entry.ionLQ_paper.beta = parseDouble(fields[19]);
        entry.ionLQ_paper.isValid = !std::isnan(entry.ionLQ_paper.alpha);
        
        // Parse LQ parameters from PIDE fit
        entry.photonLQ_fit.alpha = parseDouble(fields[20]);
        entry.photonLQ_fit.beta = parseDouble(fields[21]);
        entry.photonLQ_fit.isValid = !std::isnan(entry.photonLQ_fit.alpha);
        
        entry.ionLQ_fit.alpha = parseDouble(fields[22]);
        entry.ionLQ_fit.beta = parseDouble(fields[23]);
        entry.ionLQ_fit.isValid = !std::isnan(entry.ionLQ_fit.alpha);
        
        experiments_.push_back(entry);
    }
    
    return !experiments_.empty();
}

bool PIDEDataReader::loadPhotonRawData(const std::string& datFile) {
    std::ifstream file(datFile);
    if (!file.is_open()) {
        return false;
    }
    
    photonRawData_.clear();
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        RawSurvivalCurve curve;
        
        // Read header: expID, publicationID, publicationName, photonExpNum
        iss >> curve.expID >> curve.publicationID >> curve.publicationName >> curve.experimentNum;
        
        // Read dose-survival pairs
        double dose, sf;
        std::string token;
        
        // Check for N/A
        iss >> token;
        if (token == "N/A") {
            curve.hasData = false;
        } else {
            curve.hasData = true;
            dose = std::stod(token);
            iss >> sf;
            curve.dataPoints.push_back({dose, sf});
            
            while (iss >> dose >> sf) {
                curve.dataPoints.push_back({dose, sf});
            }
        }
        
        photonRawData_[{curve.publicationID, curve.experimentNum}] = curve;
    }
    
    return !photonRawData_.empty();
}

bool PIDEDataReader::loadIonRawData(const std::string& datFile) {
    std::ifstream file(datFile);
    if (!file.is_open()) {
        return false;
    }
    
    ionRawData_.clear();
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        RawSurvivalCurve curve;
        
        // Read header: expID, publicationID, publicationName, ionExpNum, photonExpNum
        iss >> curve.expID >> curve.publicationID >> curve.publicationName 
            >> curve.experimentNum >> curve.photonExpNum;
        
        // Read dose-survival pairs
        double dose, sf;
        std::string token;
        
        // Check for N/A
        iss >> token;
        if (token == "N/A") {
            curve.hasData = false;
        } else {
            curve.hasData = true;
            dose = std::stod(token);
            iss >> sf;
            curve.dataPoints.push_back({dose, sf});
            
            while (iss >> dose >> sf) {
                curve.dataPoints.push_back({dose, sf});
            }
        }
        
        ionRawData_[curve.expID] = curve;
    }
    
    return !ionRawData_.empty();
}

std::vector<ExperimentEntry> PIDEDataReader::getExperimentsByCellLine(const std::string& cellLine) const {
    std::vector<ExperimentEntry> result;
    for (const auto& exp : experiments_) {
        if (exp.cellLine == cellLine) {
            result.push_back(exp);
        }
    }
    return result;
}

std::vector<ExperimentEntry> PIDEDataReader::getExperimentsByIon(const std::string& ion) const {
    std::vector<ExperimentEntry> result;
    for (const auto& exp : experiments_) {
        if (exp.ion == ion) {
            result.push_back(exp);
        }
    }
    return result;
}

std::vector<ExperimentEntry> PIDEDataReader::getExperimentsByLETRange(double minLET, double maxLET) const {
    std::vector<ExperimentEntry> result;
    for (const auto& exp : experiments_) {
        if (exp.LET_keV_um >= minLET && exp.LET_keV_um <= maxLET) {
            result.push_back(exp);
        }
    }
    return result;
}

std::vector<ExperimentEntry> PIDEDataReader::getExperimentsByCellClass(char cellClass) const {
    std::vector<ExperimentEntry> result;
    for (const auto& exp : experiments_) {
        if (exp.cellClass == cellClass) {
            result.push_back(exp);
        }
    }
    return result;
}

std::vector<ExperimentEntry> PIDEDataReader::getExperimentsByCellOrigin(char cellOrigin) const {
    std::vector<ExperimentEntry> result;
    for (const auto& exp : experiments_) {
        if (exp.cellOrigin == cellOrigin) {
            result.push_back(exp);
        }
    }
    return result;
}

std::vector<ExperimentEntry> PIDEDataReader::getPhotonExperiments() const {
    // Return unique photon experiments based on (publicationID, photonExpNum)
    std::map<std::pair<int, int>, ExperimentEntry> uniquePhoton;
    for (const auto& exp : experiments_) {
        auto key = std::make_pair(exp.publicationID, exp.photonExpNum);
        if (uniquePhoton.find(key) == uniquePhoton.end()) {
            uniquePhoton[key] = exp;
        }
    }
    
    std::vector<ExperimentEntry> result;
    for (const auto& pair : uniquePhoton) {
        result.push_back(pair.second);
    }
    return result;
}

RawSurvivalCurve PIDEDataReader::getPhotonRawData(int publicationID, int photonExpNum) const {
    auto key = std::make_pair(publicationID, photonExpNum);
    auto it = photonRawData_.find(key);
    if (it != photonRawData_.end()) {
        return it->second;
    }
    return RawSurvivalCurve();
}

RawSurvivalCurve PIDEDataReader::getIonRawData(int expID) const {
    auto it = ionRawData_.find(expID);
    if (it != ionRawData_.end()) {
        return it->second;
    }
    return RawSurvivalCurve();
}

std::vector<std::string> PIDEDataReader::getUniqueCellLines() const {
    std::set<std::string> cellLines;
    for (const auto& exp : experiments_) {
        cellLines.insert(exp.cellLine);
    }
    return std::vector<std::string>(cellLines.begin(), cellLines.end());
}

std::vector<std::string> PIDEDataReader::getUniqueIons() const {
    std::set<std::string> ions;
    for (const auto& exp : experiments_) {
        ions.insert(exp.ion);
    }
    return std::vector<std::string>(ions.begin(), ions.end());
}

const ExperimentEntry* PIDEDataReader::getExperimentByID(int expID) const {
    for (const auto& exp : experiments_) {
        if (exp.expID == expID) {
            return &exp;
        }
    }
    return nullptr;
}

double PIDEDataReader::compareSurvivalCurve(const std::vector<double>& simDoses,
                                             const std::vector<double>& simSF,
                                             int expID) const {
    if (simDoses.size() != simSF.size() || simDoses.empty()) {
        return -1.0;
    }
    
    // Get experimental data
    auto rawData = getIonRawData(expID);
    if (!rawData.hasData || rawData.dataPoints.empty()) {
        // Try to get LQ parameters instead
        const auto* exp = getExperimentByID(expID);
        if (!exp || !exp->ionLQ_fit.isValid) {
            return -1.0;
        }
        
        // Compare with LQ curve
        double mse = 0.0;
        for (size_t i = 0; i < simDoses.size(); ++i) {
            double expSF = exp->ionLQ_fit.survivalFraction(simDoses[i]);
            double diff = simSF[i] - expSF;
            mse += diff * diff;
        }
        return mse / simDoses.size();
    }
    
    // Interpolate experimental data to simulation dose points
    double mse = 0.0;
    int validPoints = 0;
    
    for (size_t i = 0; i < simDoses.size(); ++i) {
        double dose = simDoses[i];
        
        // Find bracketing experimental points
        double expSF = -1.0;
        for (size_t j = 0; j < rawData.dataPoints.size() - 1; ++j) {
            if (rawData.dataPoints[j].dose_Gy <= dose && 
                rawData.dataPoints[j + 1].dose_Gy >= dose) {
                // Linear interpolation
                double d1 = rawData.dataPoints[j].dose_Gy;
                double d2 = rawData.dataPoints[j + 1].dose_Gy;
                double s1 = rawData.dataPoints[j].survivalFraction;
                double s2 = rawData.dataPoints[j + 1].survivalFraction;
                
                if (d2 - d1 > 0) {
                    double t = (dose - d1) / (d2 - d1);
                    expSF = s1 + t * (s2 - s1);
                }
                break;
            }
        }
        
        if (expSF >= 0) {
            double diff = simSF[i] - expSF;
            mse += diff * diff;
            validPoints++;
        }
    }
    
    if (validPoints == 0) return -1.0;
    return mse / validPoints;
}

LQParameters PIDEDataReader::fitLQParameters(const std::vector<double>& doses,
                                              const std::vector<double>& survivalFractions) {
    LQParameters params;
    if (doses.size() != survivalFractions.size() || doses.size() < 3) {
        return params;
    }
    
    // Fit: -ln(SF) = alpha*D + beta*D^2
    // Using least squares: minimize sum of (y_i - alpha*x_i - beta*x_i^2)^2
    // where y_i = -ln(SF_i), x_i = D_i
    
    int n = doses.size();
    double sumX = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0;
    double sumY = 0, sumXY = 0, sumX2Y = 0;
    
    for (int i = 0; i < n; ++i) {
        if (survivalFractions[i] <= 0) continue;
        
        double x = doses[i];
        double x2 = x * x;
        double y = -std::log(survivalFractions[i]);
        
        sumX += x;
        sumX2 += x2;
        sumX3 += x * x2;
        sumX4 += x2 * x2;
        sumY += y;
        sumXY += x * y;
        sumX2Y += x2 * y;
    }
    
    // Solve 2x2 system:
    // [sumX2  sumX3 ] [alpha]   [sumXY ]
    // [sumX3  sumX4 ] [beta ]   [sumX2Y]
    
    double det = sumX2 * sumX4 - sumX3 * sumX3;
    if (std::abs(det) < 1e-10) {
        // Fallback to linear fit (beta = 0)
        params.alpha = sumXY / sumX2;
        params.beta = 0;
    } else {
        params.alpha = (sumX4 * sumXY - sumX3 * sumX2Y) / det;
        params.beta = (sumX2 * sumX2Y - sumX3 * sumXY) / det;
    }
    
    // Ensure non-negative parameters
    if (params.alpha < 0) params.alpha = 0;
    if (params.beta < 0) params.beta = 0;
    
    params.isValid = true;
    return params;
}

void PIDEDataReader::printSummary() const {
    std::cout << "\n========================================" << std::endl;
    std::cout << "PIDE Database Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::cout << "Total experiments: " << experiments_.size() << std::endl;
    std::cout << "Photon raw data entries: " << photonRawData_.size() << std::endl;
    std::cout << "Ion raw data entries: " << ionRawData_.size() << std::endl;
    
    auto cellLines = getUniqueCellLines();
    std::cout << "\nUnique cell lines: " << cellLines.size() << std::endl;
    
    auto ions = getUniqueIons();
    std::cout << "Unique ions: " << ions.size() << std::endl;
    std::cout << "  Ions: ";
    for (size_t i = 0; i < ions.size() && i < 10; ++i) {
        std::cout << ions[i];
        if (i < ions.size() - 1 && i < 9) std::cout << ", ";
    }
    if (ions.size() > 10) std::cout << ", ...";
    std::cout << std::endl;
    
    // Count by cell class
    int tumorCount = 0, normalCount = 0;
    for (const auto& exp : experiments_) {
        if (exp.cellClass == 't') tumorCount++;
        else normalCount++;
    }
    std::cout << "\nTumor cell experiments: " << tumorCount << std::endl;
    std::cout << "Normal cell experiments: " << normalCount << std::endl;
    
    // Count by cell origin
    int humanCount = 0, rodentCount = 0;
    for (const auto& exp : experiments_) {
        if (exp.cellOrigin == 'h') humanCount++;
        else rodentCount++;
    }
    std::cout << "\nHuman cell experiments: " << humanCount << std::endl;
    std::cout << "Rodent cell experiments: " << rodentCount << std::endl;
    
    // LET range
    double minLET = 1e10, maxLET = 0;
    for (const auto& exp : experiments_) {
        if (exp.LET_keV_um > 0) {
            if (exp.LET_keV_um < minLET) minLET = exp.LET_keV_um;
            if (exp.LET_keV_um > maxLET) maxLET = exp.LET_keV_um;
        }
    }
    std::cout << "\nLET range: " << std::fixed << std::setprecision(2) 
              << minLET << " - " << maxLET << " keV/μm" << std::endl;
    
    std::cout << "========================================\n" << std::endl;
}

std::vector<std::string> PIDEDataReader::splitLine(const std::string& line, char delimiter) const {
    std::vector<std::string> result;
    std::stringstream ss(line);
    std::string field;
    
    bool inQuotes = false;
    std::string currentField;
    
    for (char c : line) {
        if (c == '"') {
            inQuotes = !inQuotes;
        } else if (c == delimiter && !inQuotes) {
            result.push_back(currentField);
            currentField.clear();
        } else {
            currentField += c;
        }
    }
    result.push_back(currentField);
    
    return result;
}

double PIDEDataReader::parseDouble(const std::string& str) const {
    std::string trimmed = trim(str);
    if (trimmed.empty() || trimmed == "N/A" || trimmed == "NA" || trimmed == "nan") {
        return std::nan("");
    }
    try {
        return std::stod(trimmed);
    } catch (...) {
        return std::nan("");
    }
}

int PIDEDataReader::parseInt(const std::string& str) const {
    std::string trimmed = trim(str);
    if (trimmed.empty() || trimmed == "N/A" || trimmed == "NA") {
        return 0;
    }
    try {
        return std::stoi(trimmed);
    } catch (...) {
        return 0;
    }
}

std::string PIDEDataReader::trim(const std::string& str) const {
    size_t start = str.find_first_not_of(" \t\r\n\"");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\r\n\"");
    return str.substr(start, end - start + 1);
}

} // namespace PIDE
