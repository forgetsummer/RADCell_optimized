#ifndef CellStateModel_h
#define CellStateModel_h 1
#include <map>
#include <vector>
#include <string>
#include "Cell.hh"

using namespace std;

struct CellStateModelParameter
{
  double E1;
  double E2;
  double E3;
  double sigma;
  double alpha;
  double beta;
};

class CellStateModel
{
public:
    // Set the RNG seed for the current thread (reproducible Monte Carlo)
    void SetRandomSeed(unsigned long long seed);

    void CellStateModelParameterSetup(Cell theCell);// function for setting up the simulation parameters for cell state model
    void SetMisrepairRate(const std::string& cellType, double k_error); // Set misrepair rate for a cell type
    void TissueGeometryInitialization(double xDim, double yDim, double zDim,double gridSize);
    void CellTypeInitialiation(int cellId, Cell theCell);
    void CellPositionInitialization(int cellID, double cx, double cy, double cz);
    void CellPhaseInitialization(int cellID, string cphase);// call this function for each cell to initialize cell phase
    void CellPhaseInitializationRandom(int cellID); // To initialize cell phase randomly
    void CellStateInitialization(int cellID, string cState);// function initialize cell state, S1,S2,S3
    void CellPhaseUpdate(int cellID, bool proliferationState, double deltaT, int frequency);// function updating the cell phase with time going on
    void CellStateUpdate(int cellID, int DSBNum, double integralConcentration, double deltaT, int frequency); // function updating the cell state, deltaT in unit of seconds

    std::map<int, std::string> GetCellPhase();
    std::map<int, double> GetCellPhaseDuration();
    std::map<int, double> GetCellAge();
    std::map<int, std::string> GetCellState();
    std::map<int, double> GetCellStateDuration();
    std::map<int, double> GetCellStateAge();
    std::map<int, int> GetCellAncestryID();
    std::map<int, double> GetCellPositionX();
    std::map<int, double> GetCellPositionY();
    std::map<int, double> GetCellPositionZ();
    void SetUpContactInhibition(bool considerOrNot); // determine whether considering contact inhibition or not, set up true when consider it
    
    // Priority 2: Checkpoint hold diagnostics
    std::map<int, double> GetCheckpointHoldAge(); // Get checkpoint hold age for all cells
    int GetCheckpointHoldCount(); // Get number of cells currently in checkpoint hold
    
    // Checkpoint enable/disable (Probabilistic Checkpoint Design only)
    void SetCheckpointEnabled(bool enabled); // Enable or disable checkpoint mechanism
    bool IsCheckpointEnabled() const; // Check if checkpoint is enabled
    
    // Priority 1: Soft saturation enable/disable
    void SetSoftSaturationEnabled(bool enabled); // Enable or disable soft saturation (E3 + 3.0*sigma vs E3)
    bool IsSoftSaturationEnabled() const; // Check if soft saturation is enabled
    
    CellStateModelParameter GetCellStateModelParameters(double p_sp,double p_mis_sp, double avgN);
    
private:
    double d;// grid size of cell home
    int NX;
    int NY;
    int NZ;
    void CellPhaseTransition(int cellID, double deltaT, int frequency);//function updating the cell phase with time going on
    void CellStateTransition(int cellID, double E1,double E2, double E3,double sigma, double increaseTime,bool transitionType);
    void CellStateTransitionMisrepair(int cellID, double E1,double E2, double E3,double sigma, double increaseTime,bool transitionType);
    double GaussianSampling(double mu,double sigma);  
    double GaussianSampling2Pi(double mu, double sigma);
    double GaussianCDF(double x, double mu, double sigma);
    double InstantaneousStateJumpProb(double p_sp,double increaseTime,double ObservationTime);
    double DelayedStateJumpProb(double p_sp,double increaseTime,double ObservationTime);
    
    // Helper functions for probabilistic checkpoint model
    double CheckpointEngagementProbability(double E, double E_hold, double w);
    double SampleLogNormal(double mu, double sigma);
    double GatingFunction(double t, double T_hold, double q);
    double SaturatingLogMean(double Ecur, double E_hold, double w_E);
    bool considerContactInhibition = true;
    
    struct PhaseInfo // a struct containing cell age information
    {
        std::string phase; // age phase, {G1, S, G2, M}
        double age; // time cell has been staying at this phase
        double duration; // expected time the cell will stay at this phase, determined by GaussianSampling
        double telomereFraction;
        int ancestryID;// the ancestry id for the daughter cell, daughter cell will inherite mother ancestry ID
    };
    struct StateInfo
    {
        std::string state;// state , {S1,S2,S3}
        double age;
        double duration;
        double E;// cell state energy
    };

    struct CycleInfo
    {
        double mTG1; // mean time for staying at G1
        double sigmaTG1; // sigma for G1
        double mTS;
        double sigmaTS;
        double mTG2;
        double sigmaTG2;
        double mTM;
        double sigmaTM;
    };

    struct CellStateParaInfo
    {
        double alpha;// consant for translating DSB numbers to state energy, E = alpha*DSBs
        double beta;// constant for translating integral of concentration to state energy, E = beta*C
        double sigma;
        double E1;
        double E2;
        double E3;
        double fG1;
        double fS;
        double fG2;
        double fM;
        double f1;
        double f2;
        double lambda1;
        double lambda2;
        double To12;
        double To13;
        double To21;
        double To23;
        double k_error;  // misrepair error rate constant for stochastic misrepair channel
    };

    struct PositionInfo
    {
         int i;
         int j;
         int k;
    };
    std::map<std::string,PositionInfo>CheckCellContactInhibitionCondition(int ci, int cj, int ck);
    std::map<std::string,PositionInfo>CheckCellContactInhibitionCondition_new(int i, int j, int k);
    std::vector<std::vector<std::vector<int> > > cellHomeResidence;
    std::map<int,PositionInfo > cellPositionMap;//key is cell id
    std::map<int, PhaseInfo> cellPhaseMap; // key is cell id
    std::map<int, StateInfo> cellStateMap; // map for storing the cell state, {S1, S2, S3} 
    std::map<int, Cell> cellTypeMap;
    std::map<std::string, CycleInfo> cellCycleInfoMap; // key is cell type
    std::map<std::string, CellStateParaInfo> cellStateParaInfoMap; // key is cell type
    std::map<int, bool> cellProliferativeMap;// key is cell id, value is proliferative state, true or not 
    std::map<int, double> cellExternalPerturbationEnergyMap;// map for storing the external perturbation energy, key is cell id
    std::map<int, double> cellDamageEnergyMap;// map for storing damage energy for misrepair calculation, key is cell id
    std::map<int, int> cellInitialDSBMap;// map for storing initial DSB count for misrepair calculation (matches analytical model), key is cell id
    std::map<int, double> cellS2HoldAge_h;// map for tracking checkpoint hold age (hours since cell entered S2 checkpoint), key is cell id (Priority 2)
    
    // Probabilistic checkpoint data structures
    std::map<int, double> cellT_hold_h;  // Store sampled hold duration per cell (key: cellID) for probabilistic checkpoint
    std::map<int, bool> cellCheckpointEngaged;  // Track if cell engaged checkpoint (key: cellID) for probabilistic checkpoint

    // Checkpoint parameters (Probabilistic Checkpoint Design)
    bool checkpointEnabled = true;  // Enable/disable checkpoint mechanism
    
    // Probabilistic checkpoint parameters
    double k_hold = 1.5;        // Threshold parameter: E_hold = E3 - k_hold*sigma (default: 1.5)
    double w_hold = 8.0;        // Width parameter for logistic function (default: 0.8*sigma, typically 8.0) - increased to reduce over-engagement
    double eta = 1.0;           // Diversion strength parameter (default: 1.0, full diversion)
    // Hold duration is modeled as a lognormal random variable whose *median* scales with the
    // mean cell-cycle length (T_cycle = mTG1 + mTS + mTG2 + mTM). This keeps checkpoint arrest
    // times consistent across cell types with different cycle lengths.
    //
    // Median(T_hold) = f_hold * T_cycle, where f_hold interpolates between
    // holdMedianMinCycleFraction (low damage) and holdMedianMaxCycleFraction (high damage).
    double holdMedianMinCycleFraction = 0.4;  // e.g., 0.4 * T_cycle (≈4h if T_cycle≈10h)
    double holdMedianMaxCycleFraction = 1.0;  // e.g., 1.0 * T_cycle (≈10h if T_cycle≈10h)
    // Clamp sampled holds to avoid lognormal long tails producing multi-cycle arrests.
    double holdClampMinCycleFraction  = 0.1;  // minimum hold = 0.1 * T_cycle
    double holdClampMaxCycleFraction  = 1.5;  // maximum hold = 1.5 * T_cycle
    double w_E_hold = 10.0;      // Width parameter for saturation function (≈ sigma)
    double sigmaT_hold = 0.5;    // Standard deviation for lognormal hold duration (reduced from 1.0)
    double q_gate = 1.5;        // Power parameter for gating function (default: 1.5) - reduced from 3.0 to be less aggressive
    
    // Soft saturation parameters
    bool softSaturationEnabled = true;  // Enable/disable soft saturation (E3 + margin*sigma vs E3)
    double softSaturationMargin = 0.0;  // Margin in sigma units: Ecur >= (E3 + margin*sigma) [Set to 0 = hard saturation at E3]

};

#endif
