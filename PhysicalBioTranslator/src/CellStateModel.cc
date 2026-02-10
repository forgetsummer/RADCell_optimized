#include "CellStateModel.hh"
#include <stdlib.h>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "StandardNormalDistribution.hh"

// Thread-local RNG for uniform distribution [0,1)
// - Each thread has its own engine to avoid data races under OpenMP.
// - Default seeding uses std::random_device; for reproducibility use CellStateModel::SetRandomSeed().
static std::mt19937& GetRandomEngine() {
    thread_local std::mt19937 gen;
    thread_local bool seeded = false;
    if (!seeded) {
        std::random_device rd;
        gen.seed(rd());
        seeded = true;
    }
    return gen;
}

static double UniformRand() {
    thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(GetRandomEngine());
}

void CellStateModel::SetRandomSeed(unsigned long long seed)
{
    GetRandomEngine().seed(static_cast<std::mt19937::result_type>(seed));
    // Mark this thread as explicitly seeded
    // (GetRandomEngine uses a thread_local 'seeded' flag; reseeding here is enough for reproducibility.)
}


void CellStateModel::CellStateModelParameterSetup(Cell theCell)
{
    std::string cType = theCell.GetCellType();
    cellStateParaInfoMap[cType].alpha =theCell.GetCellStateParameters().alpha;// set up the alpha,for deltaE = alpha*DSB_num
    cellStateParaInfoMap[cType].beta = theCell.GetCellStateParameters().beta;// set up for the beta, for deltaE = beta* integralConcentration
    cellStateParaInfoMap[cType].E1 = theCell.GetCellStateParameters().E1;// mean state energy of S1 state
    cellStateParaInfoMap[cType].E2 = theCell.GetCellStateParameters().E2;// mean state energy of S2 state
    cellStateParaInfoMap[cType].E3 = theCell.GetCellStateParameters().E3;// mean state energy of S3 state
    cellStateParaInfoMap[cType].sigma = theCell.GetCellStateParameters().sigma;
    cellStateParaInfoMap[cType].fG1 = theCell.GetCellCycleParameters().fG1;
    cellStateParaInfoMap[cType].fS = theCell.GetCellCycleParameters().fS;
    cellStateParaInfoMap[cType].fG2 = theCell.GetCellCycleParameters().fG2;
    cellStateParaInfoMap[cType].fM = theCell.GetCellCycleParameters().fM;
    cellStateParaInfoMap[cType].f1=theCell.GetCellDNADamageRepairParameters().f1;
    cellStateParaInfoMap[cType].f2= theCell.GetCellDNADamageRepairParameters().f2;
    cellStateParaInfoMap[cType].lambda1= theCell.GetCellDNADamageRepairParameters().lambda1;
    cellStateParaInfoMap[cType].lambda2 = theCell.GetCellDNADamageRepairParameters().lambda2;
    cellStateParaInfoMap[cType].To12 = theCell.GetCellStateParameters().To12;
    cellStateParaInfoMap[cType].To13 = theCell.GetCellStateParameters().To13;
    cellStateParaInfoMap[cType].To21 = theCell.GetCellStateParameters().To21;
    cellStateParaInfoMap[cType].To23 = theCell.GetCellStateParameters().To23;
    // k_error is set separately via SetMisrepairRate (default is 0.0)

    cellCycleInfoMap[cType].mTG1 =theCell.GetCellCycleParameters().mTG1;
    cellCycleInfoMap[cType].sigmaTG1 = theCell.GetCellCycleParameters().sigmaTG1;
    cellCycleInfoMap[cType].mTS = theCell.GetCellCycleParameters().mTS;
    cellCycleInfoMap[cType].sigmaTS = theCell.GetCellCycleParameters().sigmaTS;
    cellCycleInfoMap[cType].mTG2 = theCell.GetCellCycleParameters().mTG2;
    cellCycleInfoMap[cType].sigmaTG2 = theCell.GetCellCycleParameters().sigmaTG2;
    cellCycleInfoMap[cType].mTM = theCell.GetCellCycleParameters().mTM;
    cellCycleInfoMap[cType].sigmaTM = theCell.GetCellCycleParameters().sigmaTM;
}

void CellStateModel::SetMisrepairRate(const std::string& cellType, double k_error)
{
    // Set the misrepair rate for a specific cell type
    // This enables stochastic misrepair channel in S2->S3 transitions
    cellStateParaInfoMap[cellType].k_error = k_error;
}

void CellStateModel::TissueGeometryInitialization(double xDim, double yDim, double zDim, double gridSize)
{
    d = gridSize;
    double cellHomeSizeX = gridSize;
    double cellHomeSizeY = gridSize;
    double cellHomeSizeZ = gridSize;
    int N_X=ceil(xDim/cellHomeSizeX);
    int N_Y=ceil(yDim/cellHomeSizeY);
    int N_Z=ceil(zDim/cellHomeSizeZ);
    if (N_Z==0)// in case of  ceil(0)=0
    {
        N_Z=1;
    }
    if (N_Y==0)
    {
        N_Y=1;
    }
    if (N_X==0)
    {
        N_X=1;
    }

    NX = N_X;
    NY = N_Y;
    NZ = N_Z;

    for (int i=0;i<N_X;i++)
    {
        cellHomeResidence.push_back(std::vector<std::vector<int> >());
        for (int j=0;j<N_Y;j++)
        {
            (cellHomeResidence[i]).push_back(std::vector<int>());
            for (int k=0;k<N_Z;k++)
            {
                (cellHomeResidence[i][j]).push_back(0);// intialize cellHomeResidence
            }
        }
    }

}

void CellStateModel::CellTypeInitialiation(int cellId, Cell theCell)
{
    cellTypeMap[cellId] = theCell;
}


double CellStateModel::GaussianSampling(double mu, double sigma)
{
    double Xf;
    double pi = 3.1415926535897;
    Xf =mu+ sigma*sqrt(-2*log(UniformRand()))*cos(pi/2.0*UniformRand());
    return Xf;
    
}
double CellStateModel::GaussianSampling2Pi(double mu, double sigma)
{
    double Xf;
    double pi = 3.1415926535897;
    Xf =mu+ sigma*sqrt(-2*log(UniformRand()))*cos(2.0*pi*UniformRand());
    return Xf;

}


double CellStateModel::GaussianCDF(double x, double mu, double sigma)
{
    double x_stdNORM=(x-mu)/sigma;
    double p_CDF;
    p_CDF = 0.5*(1+erf(x_stdNORM/sqrt(2)));
    return p_CDF;// return the probability of N(x,mu, sigma)

}

// Helper functions for probabilistic checkpoint model
double CellStateModel::CheckpointEngagementProbability(double E, double E_hold, double w)
{
    // Logistic function: π(E) = 1 / (1 + exp(-(E - E_hold)/w))
    if (w <= 0.0) return (E >= E_hold) ? 1.0 : 0.0;
    return 1.0 / (1.0 + exp(-(E - E_hold) / w));
}

double CellStateModel::SampleLogNormal(double mu, double sigma)
{
    // Sample lognormal distribution: if Z ~ N(0,1), then exp(mu + sigma*Z) ~ LogNormal(mu, sigma^2)
    double z = GaussianSampling(0.0, 1.0);
    return exp(mu + sigma * z);
}

double CellStateModel::GatingFunction(double t, double T_hold, double q)
{
    // Gating function: g(t; T_hold) = (t/T_hold)^q for t < T_hold, else 1
    if (T_hold <= 0.0) return 1.0;
    if (t >= T_hold) return 1.0;
    if (t <= 0.0) return 0.0;
    return pow(t / T_hold, q);
}

double CellStateModel::SaturatingLogMean(double Ecur, double E_hold, double w_E)
{
    // Saturating function: s(E) = 1 / (1 + exp(-(E - E_hold)/w_E))
    // Returns value in [0, 1] for use in μ_E calculation
    // This bounds T_hold to realistic timescales by saturating the log-mean parameter
    if (w_E <= 0.0) return (Ecur >= E_hold) ? 1.0 : 0.0;
    return 1.0 / (1.0 + exp(-(Ecur - E_hold) / w_E));
}




void CellStateModel::CellPositionInitialization(int cellID, double cx, double cy, double cz)
{
//     cout<<"cx ="<<cx<<endl;
//     cout<<"cy ="<<cy<<endl;
//     cout<<"cz ="<<cz<<endl;
    double c_i = floor(cx/d)-1;
    double c_j = floor(cy/d)-1;
    double c_k = floor(cz/d)-1;            
    if (c_i==-1)
    {
        c_i=0;
    }
    if (c_j==-1)
    {
        c_j=0;
    }
    if (c_k==-1)
    {
        c_k=0;
    }
    cellHomeResidence[c_i][c_j][c_k] = 1; // this cell home was occupied by the cell
    PositionInfo cellPosition;
    cellPosition.i = c_i;
    cellPosition.j = c_j;
    cellPosition.k = c_k;
    cellPositionMap[cellID] = cellPosition;

}


void CellStateModel::CellPhaseInitialization(int cellID, string cphase)
{
    PhaseInfo cellPhase;
    string cellType =cellTypeMap[cellID].GetCellType();
    if (cellCycleInfoMap.find(cellType)==cellCycleInfoMap.end())
    {
        cout<<endl<<"ERROR: "<<"THE CELL TYPE IS NOT EXISTING!"<<endl;
        exit(1);
    }

    
    if ((cphase!="G1") && ( cphase!="S") && (cphase!="G2")&& (cphase!="M") && (cphase!="G0"))
    {
        cout<<cphase<<endl;
        cout<<endl<<"ERROR: "<<"THIS PHASE DOES NOT EXIST!"<<endl;
        exit(1);
    }
    else
    {
        cellPhase.phase = cphase;
        cellPhase.age = 0;
        cellPhase.telomereFraction=1;// initially telomereFraction is 1
        cellPhase.ancestryID = cellID;
        if (cphase == "G1")
        {
            cellPhase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTG1,cellCycleInfoMap.at(cellType).sigmaTG1);
        }
        if (cphase =="S")
        {
            cellPhase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTS,cellCycleInfoMap.at(cellType).sigmaTS);
        }
        if (cphase == "G2")
        {
            cellPhase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTG2,cellCycleInfoMap.at(cellType).sigmaTG2);
        }
        if (cphase == "M")
        {
            cellPhase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTM,cellCycleInfoMap.at(cellType).sigmaTM);
        }
        if (cphase == "G0")
        {
            cellPhase.duration = 1E+18;// just initially give a very large number which means that cell will not proliferate
        }
    }
    
    cellPhaseMap[cellID] = cellPhase;
    
}


void CellStateModel::CellPhaseInitializationRandom(int cellID)
{
    string cellType = cellTypeMap[cellID].GetCellType();
    double eta = UniformRand();
    string phase;
    if (eta<=0.2)
    {
        phase = "G1";
    }
    if (0.2<eta && eta<=0.4)
    {
        phase = "G2";
    }
    if (0.4< eta && eta<=0.6)
    {
        phase = "S";
    }
    if (0.6<eta && eta <=0.8)
    {
        phase = "M";
    }
    if (0.8<eta && eta<=1.0)
    {
        phase = "G0";
    }
    PhaseInfo cellPhase;
    cellPhase.age = 0;
    cellPhase.phase = phase;
    cellPhase.telomereFraction = 1;
    cellPhase.ancestryID = cellID;
    
    if (phase == "G1")
    {
//         cout<<"meam time G1 "<<cellCycleInfoMap.at(cellType).mTG1<<" sigma time "<<cellCycleInfoMap.at(cellType).sigmaTG1<<endl;
        cellPhase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTG1,cellCycleInfoMap.at(cellType).sigmaTG1);
    }
    if (phase =="S")
    {
//         cout<<"meam time S "<<cellCycleInfoMap.at(cellType).mTS<<" sigma time "<<cellCycleInfoMap.at(cellType).sigmaTS<<endl;
        cellPhase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTS,cellCycleInfoMap.at(cellType).sigmaTS);
    }
    if (phase == "G2")
    {
//         cout<<"meam time G2 "<<cellCycleInfoMap.at(cellType).mTG2<<" sigma time "<<cellCycleInfoMap.at(cellType).sigmaTG2<<endl;
        cellPhase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTG2,cellCycleInfoMap.at(cellType).sigmaTG2);
    }
    if (phase == "M")
    {
//         cout<<"meam time M "<<cellCycleInfoMap.at(cellType).mTM<<" sigma time "<<cellCycleInfoMap.at(cellType).sigmaTM<<endl;
        cellPhase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTM,cellCycleInfoMap.at(cellType).sigmaTM);
    }
    if (phase == "G0")
    {
        cellPhase.duration = 1E+18;
    }
    
    cellPhaseMap[cellID] = cellPhase;

}

void CellStateModel::CellStateInitialization(int cellID,string cState)
{
    StateInfo cellState;
    cellState.state = cState;
    cellState.age =0;
    // IMPORTANT: Initialize state energy deterministically.
    // If left uninitialized, st.E can become NaN/garbage and silently prevent transitions,
    // severely biasing survival (especially at 0 Gy and after division).
    cellState.E = 0.0;
    if ((cState!="S1") && ( cState!="S2")&& (cState!="S3"))
    {
        cout<<cState<<endl;
        cout<<endl<<"ERROR: "<<"THIS STATE DOES NOT EXIST!"<<endl;
        cout<<"THE POSSIBLE STATES ARE: S1, S2, S3"<<endl;
        exit(1);
    }

    if(cState=="S1")
    {
        cellState.duration=1E+18;
    }
    if(cState=="S2")
    {
        double E2 = cellTypeMap[cellID].GetCellStateParameters().E2;
        double sigma = cellTypeMap[cellID].GetCellStateParameters().sigma;
        double f1 = cellTypeMap[cellID].GetCellDNADamageRepairParameters().f1;
        double f2 = cellTypeMap[cellID].GetCellDNADamageRepairParameters().f2;
        double lambda1 = cellTypeMap[cellID].GetCellDNADamageRepairParameters().lambda1;
        double lambda2 = cellTypeMap[cellID].GetCellDNADamageRepairParameters().lambda2;

        double E20= GaussianSampling2Pi(E2,sigma);// Sample a E20 of S2 state
        while (E20<0)
        {
            E20= GaussianSampling2Pi(E2,sigma);// Sample a E20 of S2 state
        }
        double pi = 3.1415926535897;
        double arrestedDuration = E20/sqrt(2*pi)/sigma*(f1/lambda1+f2/lambda2);
        cellState.duration = arrestedDuration;
        cellState.E = E20;
    }
    if(cState=="S22")
    {
        cellState.duration = 1E+18;
    }
    if(cState=="S3")
    {
        cellState.duration = 1E+18;
    }

    cellStateMap[cellID] = cellState;
    cellExternalPerturbationEnergyMap[cellID] = 0;// at the beginning, there is no external perturbation energy

}

std::map<std::string, CellStateModel::PositionInfo > CellStateModel::CheckCellContactInhibitionCondition_new(int i, int j, int k)
{
    std::map<std::string, PositionInfo> possible;

    // ---- 0) Defensive bounds check for the input coordinate ----
    if (i < 0 || i >= NX || j < 0 || j >= NY || k < 0 || k >= NZ)
        return possible;

    auto try_add = [&](const char* key, int ni, int nj, int nk)
    {
        // bounds check neighbor
        if (ni < 0 || ni >= NX || nj < 0 || nj >= NY || nk < 0 || nk >= NZ)
            return;

        // occupancy check (0 means empty in your logic)
        if (cellHomeResidence[ni][nj][nk] == 0)
        {
            PositionInfo p;
            p.i = ni; p.j = nj; p.k = nk;
            possible.emplace(key, p);
        }
    };

    // ---- 1) Check each axis only if that axis has >1 cell ----
    if (NX > 1) { try_add("right", i + 1, j,     k); try_add("left",  i - 1, j,     k); }
    if (NY > 1) { try_add("front", i,     j + 1, k); try_add("back",  i,     j - 1, k); }
    if (NZ > 1) { try_add("up",    i,     j,     k + 1); try_add("down", i,   j,     k - 1); }

    return possible;
}



map< string, CellStateModel::PositionInfo > CellStateModel::CheckCellContactInhibitionCondition(int i, int j, int k)
{
    std::map<std::string,PositionInfo> possibleCellHomeForMitosis;
    if (NX>1&&NY>1&&NZ==1)// if only in X-Y plane
    {
        if (i==0)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
        }
        if (i==NX-1)
        {
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (0<i && i<NX-1)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (j==0)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            
        }

        if (j==NY-1)
        {
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
        }
                                                    
        if (0<j && j<NY-1)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
            
        }
        
    }

    if (NY>1&&NZ>1&&NX==1)//if in Y-Z plane
    {
        if (j==0)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            
        }

        if (j==NY-1)
        {
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
        }
                                                    
        if (0<j && j<NY-1)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
            
        }
        if (k==0)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }
        }
        if (k==NZ-1)
        {
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        }
        if (0<k && k<NZ-1)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }  
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        } 
        
    }
    if (NX>1&&NZ>1&&NY==1)// if in X-Z plane
    {
            if (i==0)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
        }
        if (i==NX-1)
        {
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (0<i && i<NX-1)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (k==0)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }
        }
        if (k==NZ-1)
        {
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        }
        if (0<k && k<NZ-1)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }  
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        } 
        
    }
    if (NX>1&&NY>1&&NZ>1)// 3D case
    {
        if (i==0)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
        }
        if (i==NX-1)
        {
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (0<i && i<NX-1)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (j==0)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            
        }

        if (j==NY-1)
        {
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
        }
                                                    
        if (0<j && j<NY-1)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
            
        }
        if (k==0)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }
        }
        if (k==NZ-1)
        {
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        }
        if (0<k && k<NZ-1)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }  
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        } 
        
    }
    return possibleCellHomeForMitosis;
    
}

void CellStateModel::CellPhaseUpdate(int cellID, bool proliferationState, double deltaT, int frequency)
{
    //correct fix: move the guard ABOVE any [] access

        // 1) Guard first (BEFORE any operator[] access)
    auto itPhase = cellPhaseMap.find(cellID);
    auto itState = cellStateMap.find(cellID);
    auto itPos   = cellPositionMap.find(cellID);
    auto itT = cellTypeMap.find(cellID);

    if (itPhase == cellPhaseMap.end()) return;
    if (itState == cellStateMap.end()) return;
    if (itPos   == cellPositionMap.end()) return;
    if (itT == cellTypeMap.end()) return;


	    // 2) Safe read (no insertion)
	    const std::string& state = itState->second.state;

	    // Cell-cycle progression policy:
	    // - S1: proliferative/clonogenic
	    // - S2: repair/arrest (NON-proliferative)
	    // - S3: dead (NON-proliferative)
	    //
	    // Allowing S2 cells to keep cycling (even after a finite hold) lets them reach mitosis
	    // and "reset" damage on division, which can strongly overpredict survival.
	    const bool proliferativeFromCellState = (state == "S1");

	    // 3) Write is OK; here operator[] is acceptable if you want it to exist for all cells.
	    // If you want to avoid insertion too, use find/at pattern.
	    cellProliferativeMap[cellID] = (proliferationState && proliferativeFromCellState);//     cellProliferativeMap[cellID] =(proliferationState);//when the nutrient condtion and healthy state all true, then it goes to proliferation

    CellPhaseTransition(cellID, deltaT, frequency);
//     cout<<"CELL: "<<cellID<<" : "<<"AGE IS: "<<cellPhaseMap[cellID].age<<" PHASE IS: "<<cellPhaseMap[cellID].phase<<endl;
}

void CellStateModel::CellPhaseTransition(int cellID, double deltaT, int frequency)
{
    // ---- 0) Guard: do nothing if cell was deleted earlier this tick ----
    auto itType = cellTypeMap.find(cellID);
    auto itPro  = cellProliferativeMap.find(cellID);
    auto itPh   = cellPhaseMap.find(cellID);
    auto itPos  = cellPositionMap.find(cellID);
    auto itSt   = cellStateMap.find(cellID);

    if (itType == cellTypeMap.end()) return;
    if (itPro  == cellProliferativeMap.end()) return;
    if (itPh   == cellPhaseMap.end()) return;
    if (itPos  == cellPositionMap.end()) return;
    if (itSt   == cellStateMap.end()) return;

    const std::string cellType = itType->second.GetCellType();
    const double increaseTime  = deltaT * frequency / 3600.0; // seconds -> hours

    // If not proliferative, do nothing (your design)
    if (!itPro->second) return;

    // Use references for convenience (safe because we already checked find())
    PhaseInfo& phase = itPh->second;
    PositionInfo& pos = itPos->second;

    // ---- 1) Check telomere condition (your threshold) ----
    const bool allowableInCellTelomereLength = (phase.telomereFraction > 8.88E-16);

    // If telomere too short, go to G0
    if (!allowableInCellTelomereLength)
    {
        phase.phase = "G0";
        phase.age = 0;
        phase.duration = 1E+18;
        return;
    }

    // ---- Helper lambda to check if space is available for division ----
    auto hasSpaceForDivision = [&]() -> bool {
        if (!considerContactInhibition) return true;
        std::map<std::string, PositionInfo> possibleHomes = 
            CheckCellContactInhibitionCondition(pos.i, pos.j, pos.k);
        return !possibleHomes.empty();
    };

    // ---- 2) Phase transitions ----
    
    // G0 -> G1: Cell wakes up and enters cell cycle
    if (phase.phase == "G0")
    {
        phase.phase = "G1";
        phase.age = 0;
        phase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTG1,
                                          cellCycleInfoMap.at(cellType).sigmaTG1);
        return;
    }

    // G1 -> S: Check contact inhibition before committing to DNA synthesis
    // This is the first checkpoint - cells sense crowding before S phase
    if (phase.phase == "G1")
    {
        phase.age += increaseTime;
        if (phase.age > phase.duration)
        {
            // Check contact inhibition at G1->S transition
            if (hasSpaceForDivision())
            {
                // Space available - proceed to S phase
                phase.phase = "S";
                phase.age = 0;
                phase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTS,
                                                  cellCycleInfoMap.at(cellType).sigmaTS);
            }
            else
            {
                // No space - enter quiescence (G0)
                phase.phase = "G0";
                phase.age = 0;
                phase.duration = 1E+18;  // Long duration in G0
            }
        }
        return;
    }

    // S -> G2: DNA synthesis complete, proceed to G2
    // No contact inhibition check here - cell is already committed
    if (phase.phase == "S")
    {
        phase.age += increaseTime;
        if (phase.age > phase.duration)
        {
            phase.phase = "G2";
            phase.age = 0;
            phase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTG2,
                                              cellCycleInfoMap.at(cellType).sigmaTG2);
        }
        return;
    }

    // G2 -> M: Check contact inhibition before committing to mitosis
    // This is the second checkpoint - cells verify space before mitosis
    if (phase.phase == "G2")
    {
        phase.age += increaseTime;
        if (phase.age > phase.duration)
        {
            // Check contact inhibition at G2->M transition
            if (hasSpaceForDivision())
            {
                // Space available - proceed to M phase
                phase.phase = "M";
                phase.age = 0;
                phase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTM,
                                                  cellCycleInfoMap.at(cellType).sigmaTM);
            }
            else
            {
                // No space - enter quiescence (G0)
                phase.phase = "G0";
                phase.age = 0;
                phase.duration = 1E+18;
            }
        }
        return;
    }

    // M phase: Cell division (mitosis)
    // At this point, space was already verified at G2->M transition
    if (phase.phase == "M")
    {
        phase.age += increaseTime;

        if (phase.age <= phase.duration)
            return;

        // ---- Priority 3: Mitotic Catastrophe Check ----
        // Check residual damage before allowing division (only for S1/S2 cells)
        if (itSt != cellStateMap.end()) {
            const std::string& state = itSt->second.state;
            if (state == "S1" || state == "S2") {  // Only check viable cells
                double Ecur = itSt->second.E;
                
                // Use phase-specific E3 for mitotic threshold
                const double E3 = cellStateParaInfoMap.at(cellType).E3;
                const double sigma = cellStateParaInfoMap.at(cellType).sigma;
                double E3_mitotic = E3 * 0.7;  // Lower threshold for mitotic catastrophe
                
                // Use soft saturation (consistent with Priority 1)
                // Hard saturation: Ecur >= E3, Soft saturation: Ecur >= E3 + margin*sigma
                double p_death = 0.0;
                double mitotic_threshold = softSaturationEnabled ? 
                    (E3_mitotic + softSaturationMargin * sigma) : E3_mitotic;
                if (Ecur >= mitotic_threshold) {
                    p_death = 1.0;
                } else {
                    double E_mitotic = -std::fabs(Ecur - E3_mitotic) / (2.0 * sigma);
                    p_death = 2.0 * GaussianCDF(E_mitotic, 0, 1);
                    p_death = std::max(0.0, std::min(1.0, p_death));  // Clamp to [0,1]
                }
                
                if (UniformRand() < p_death) {
                    // Mitotic catastrophe - cell dies instead of dividing
                    itSt->second.state = "S3";
                    // Clean up checkpoint tracking (avoid stale entries in diagnostics)
                    cellS2HoldAge_h.erase(cellID);
                    cellT_hold_h.erase(cellID);
                    cellCheckpointEngaged.erase(cellID);
                    // Don't divide - cell is dead
                    return;
                }
            }
        }

        // ---- Division happens now ----
        // Recompute possible positions (space may have changed during M phase)
        std::map<std::string, PositionInfo> possibleCellHomeForMitosis;
        if (considerContactInhibition)
        {
            possibleCellHomeForMitosis = CheckCellContactInhibitionCondition(pos.i, pos.j, pos.k);
        }
        else
        {
            // No contact inhibition - allow division at original position
            possibleCellHomeForMitosis.emplace("original", pos);
        }

        // If no space available now (rare - space was checked at G2->M), go to G0
        if (possibleCellHomeForMitosis.empty())
        {
            phase.phase = "G0";
            phase.age = 0;
            phase.duration = 1E+18;
            return;
        }

        // Find max existing cellID
        int cellID_max = 0;
        for (auto it = cellPhaseMap.begin(); it != cellPhaseMap.end(); ++it)
            cellID_max = std::max(cellID_max, it->first);

        const int id1 = cellID_max + 1;
        const int id2 = cellID_max + 2;

        // Cache mother info before erase
        const double motherTel = phase.telomereFraction;
        const int motherAncestry = phase.ancestryID;
        const PositionInfo motherPos = pos;
        const Cell motherCellType = itType->second;

        // Daughter 1: takes mother's position
        cellPhaseMap[id1].phase = "G0";
        cellPhaseMap[id1].age = 0;
        cellPhaseMap[id1].duration = 1E+18;
        cellPhaseMap[id1].telomereFraction = 0.5 * motherTel;
        cellPhaseMap[id1].ancestryID = motherAncestry;

        cellStateMap[id1].state = "S1";
        cellStateMap[id1].age = 0;
        cellStateMap[id1].duration = 1E+18;
        cellStateMap[id1].E = 0.0;

        cellPositionMap[id1] = motherPos;
        cellHomeResidence[motherPos.i][motherPos.j][motherPos.k] = 1;
        cellTypeMap[id1] = motherCellType;

        // Daughter 2: random position from available spaces
        const int n  = (int)possibleCellHomeForMitosis.size();
        const int ii = (int)std::floor(UniformRand() * n);
        auto pick = possibleCellHomeForMitosis.begin();
        std::advance(pick, ii);

        cellPhaseMap[id2].phase = "G0";
        cellPhaseMap[id2].age = 0;
        cellPhaseMap[id2].duration = 1E+18;
        cellPhaseMap[id2].telomereFraction = 0.5 * motherTel;
        cellPhaseMap[id2].ancestryID = motherAncestry;

        cellStateMap[id2].state = "S1";
        cellStateMap[id2].age = 0;
        cellStateMap[id2].duration = 1E+18;
        cellStateMap[id2].E = 0.0;

        cellPositionMap[id2] = pick->second;
        cellHomeResidence[pick->second.i][pick->second.j][pick->second.k] = 1;
        cellTypeMap[id2] = motherCellType;

        // Delete mother cell
        cellPhaseMap.erase(cellID);
        cellStateMap.erase(cellID);
        cellPositionMap.erase(cellID);
        cellTypeMap.erase(cellID);
        cellProliferativeMap.erase(cellID);
        cellExternalPerturbationEnergyMap.erase(cellID);
        cellDamageEnergyMap.erase(cellID);
        cellInitialDSBMap.erase(cellID);
        cellS2HoldAge_h.erase(cellID);
        cellT_hold_h.erase(cellID);
        cellCheckpointEngaged.erase(cellID);

        return;
    }

    // If phase string is unexpected, do nothing
    return;
}

void CellStateModel::CellStateUpdate(int cellID,
                                     int DSBNum,
                                     double integralConcentration,
                                     double deltaT,
                                     int frequency)
{
    // ---- 0) Guard first: cell might have been deleted earlier this tick ----
    auto itType  = cellTypeMap.find(cellID);
    auto itPhase = cellPhaseMap.find(cellID);
    auto itState = cellStateMap.find(cellID);

    if (itType  == cellTypeMap.end())  return;
    if (itPhase == cellPhaseMap.end()) return;
    if (itState == cellStateMap.end()) return;

    const std::string cellType = itType->second.GetCellType();

    // seconds -> hours
    const double increaseTime = deltaT * frequency / 3600.0;

    const double alpha = cellStateParaInfoMap.at(cellType).alpha;
    const double beta  = cellStateParaInfoMap.at(cellType).beta;

    // external perturbation energy (store/update per cell)
    cellExternalPerturbationEnergyMap[cellID] = alpha * DSBNum + beta * integralConcentration;

    auto& st = itState->second;

    const double dt_h = increaseTime; // hours
    const double dEinj = cellExternalPerturbationEnergyMap.at(cellID);

    // bi-exponential repair parameters (must be available in your param struct)
    double f1 = cellStateParaInfoMap.at(cellType).f1;
    double f2 = cellStateParaInfoMap.at(cellType).f2;
    const double lambda1 = cellStateParaInfoMap.at(cellType).lambda1; // 1/hour
    const double lambda2 = cellStateParaInfoMap.at(cellType).lambda2; // 1/hour

    // normalize weights defensively (recommended)
    const double s = f1 + f2;
    if (s > 0.0) { f1 /= s; f2 /= s; } else { f1 = 1.0; f2 = 0.0; }

    // combined decay factor over this step
    const double decay =
        f1 * std::exp(-lambda1 * dt_h) +
        f2 * std::exp(-lambda2 * dt_h);

    // update state energy
    st.E = st.E * decay + dEinj;

    // Store initial DSB count (legacy/diagnostic; misrepair now uses current E)
    if (DSBNum > 0) {
        // New damage: store initial DSB count (for misrepair calculation)
        // Only set if not already set (to preserve initial value)
        if (cellInitialDSBMap.find(cellID) == cellInitialDSBMap.end()) {
            cellInitialDSBMap[cellID] = DSBNum;
        }
    }

    // instantaneous vs delayed transition type
    const bool transitionType = (DSBNum > 0);

    // baseline state energies
    const double E1    = cellStateParaInfoMap.at(cellType).E1;
    const double E2    = cellStateParaInfoMap.at(cellType).E2;
    const double E3    = cellStateParaInfoMap.at(cellType).E3;
    const double sigma = cellStateParaInfoMap.at(cellType).sigma;
    const double k_error = cellStateParaInfoMap.at(cellType).k_error;  // Check if misrepair is enabled

    // current phase string (safe via iterator)
    const std::string& ph = itPhase->second.phase;

    // Choose transition function: use misrepair if k_error > 0
    bool useMisrepair = (k_error > 1e-10);

    // ---- Checkpoint hold logic (Probabilistic Checkpoint) ----
    // Increment checkpoint hold age for cells in S2 checkpoint hold
    if (checkpointEnabled) {
        if (itState != cellStateMap.end() && itState->second.state == "S2") {
            auto itHoldAge = cellS2HoldAge_h.find(cellID);
            if (itHoldAge != cellS2HoldAge_h.end()) {
                // Cell is in checkpoint hold - increment hold age
                itHoldAge->second += increaseTime;
            }
        }
    }

    // ---- 1) Dispatch by phase (keep your original formulas) ----
    if (ph == "G1")
    {
        if (useMisrepair) {
            CellStateTransitionMisrepair(cellID, E1, E2, E3, sigma, increaseTime, transitionType);
        } else {
            CellStateTransition(cellID, E1, E2, E3, sigma, increaseTime, transitionType);
        }
    }
    else if (ph == "S")
    {
        const double fS = cellStateParaInfoMap.at(cellType).fS;

        const double E1_S = 0.0;
        const double E3_S = std::sqrt(E3 * E3 - 8.0 * sigma * sigma * std::log(fS));
        const double E2_S = E3_S - std::sqrt((E3 - E2) * (E3 - E2) - 8.0 * sigma * sigma * std::log(fS));

        if (useMisrepair) {
            CellStateTransitionMisrepair(cellID, E1_S, E2_S, E3_S, sigma, increaseTime, transitionType);
        } else {
            CellStateTransition(cellID, E1_S, E2_S, E3_S, sigma, increaseTime, transitionType);
        }
    }
    else if (ph == "G2")
    {
        const double fG2 = cellStateParaInfoMap.at(cellType).fG2;

        const double E1_G2 = 0.0;
        const double E3_G2 = std::sqrt(E3 * E3 - 8.0 * sigma * sigma * std::log(fG2));
        const double E2_G2 = E3_G2 - std::sqrt((E3 - E2) * (E3 - E2) - 8.0 * sigma * sigma * std::log(fG2));

        if (useMisrepair) {
            CellStateTransitionMisrepair(cellID, E1_G2, E2_G2, E3_G2, sigma, increaseTime, transitionType);
        } else {
            CellStateTransition(cellID, E1_G2, E2_G2, E3_G2, sigma, increaseTime, transitionType);
        }
    }
    else if (ph == "M")
    {
        const double fM = cellStateParaInfoMap.at(cellType).fM;

        const double E1_M = 0.0;
        const double E3_M = std::sqrt(E3 * E3 - 8.0 * sigma * sigma * std::log(fM));
        const double E2_M = E3_M - std::sqrt((E3 - E2) * (E3 - E2) - 8.0 * sigma * sigma * std::log(fM));

        if (useMisrepair) {
            CellStateTransitionMisrepair(cellID, E1_M, E2_M, E3_M, sigma, increaseTime, transitionType);
        } else {
            CellStateTransition(cellID, E1_M, E2_M, E3_M, sigma, increaseTime, transitionType);
        }
    }
    else if (ph == "G0")
    {
        // keep your original "treat G0 like G1" idea
        if (useMisrepair) {
            CellStateTransitionMisrepair(cellID, 0.0, E2, E3, sigma, increaseTime, transitionType);
        } else {
            CellStateTransition(cellID, 0.0, E2, E3, sigma, increaseTime, transitionType);
        }
    }
    else
    {
        // unknown phase string -> do nothing
        return;
    }
}


void CellStateModel::CellStateTransition(int cellID,
                                        double E1, double E2, double E3,
                                        double sigma,
                                        double increaseTime,
                                        bool transitionType) // true=instantaneous, false=delayed
{
    // ---- 0) Guard: cell might have been deleted earlier this tick ----
    auto itState = cellStateMap.find(cellID);
    if (itState == cellStateMap.end())
        return;

    // Reference (alias) to the map-stored StateInfo (no copy)
    auto& st = itState->second;

    // Helper to keep probabilities sane
    auto clamp01 = [](double x) {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    };


    auto itType  = cellTypeMap.find(cellID);
    const std::string cellType = itType->second.GetCellType();
    const double To12 =cellStateParaInfoMap.at(cellType).To12;
    const double To13 = cellStateParaInfoMap.at(cellType).To13;
    const double To21 = cellStateParaInfoMap.at(cellType).To21;
    const double To23 = cellStateParaInfoMap.at(cellType).To23;

    // =========================
    // S1 transitions
    // =========================
    if (st.state == "S1")
    {
        const double Ecur = st.E;

        // Overlap-based cumulative probabilities over observation window To
        double ES1ToS2 = -std::fabs(Ecur - E2) / (2.0 * sigma);
        double PS1ToS2 = 2.0 * GaussianCDF(ES1ToS2, 0, 1);

        double ES1ToS3 = -std::fabs(Ecur - E3) / (2.0 * sigma);
        double PS1ToS3 = 2.0 * GaussianCDF(ES1ToS3, 0, 1);

        // Soft saturation: only trigger at extreme margin (Priority 1)
        // Hard saturation: Ecur >= E3, Soft saturation: Ecur >= E3 + margin*sigma
        if (softSaturationEnabled) {
            if (Ecur >= E3 + softSaturationMargin * sigma) PS1ToS3 = 1.0;
        } else {
            if (Ecur >= E3) PS1ToS3 = 1.0;  // Hard saturation
        }

        // Map cumulative probability -> per-step probability
        double p12_step = 0.0;
        double p13_step = 0.0;
        if (transitionType) {
            p12_step = InstantaneousStateJumpProb(PS1ToS2, increaseTime, To12);
            p13_step = InstantaneousStateJumpProb(PS1ToS3, increaseTime, To13);
        } else {
            p12_step = DelayedStateJumpProb(PS1ToS2, increaseTime, To12);
            p13_step = DelayedStateJumpProb(PS1ToS3, increaseTime, To13);
        }

        p12_step = clamp01(p12_step);
        p13_step = clamp01(p13_step);

        // Probabilistic checkpoint: modify probabilities if enabled
        if (checkpointEnabled) {
            // Calculate checkpoint engagement probability
            const double E_hold = E3 - k_hold * sigma;
            const double pi_E = CheckpointEngagementProbability(Ecur, E_hold, w_hold);
            
            // Store baseline probabilities
            const double p12_0 = p12_step;
            const double p13_0 = p13_step;
            
            // Apply marginalization formulas
            // Use effective_eta = eta directly (default eta = 1.0)
            double effective_eta = eta;
            
            p13_step = (1.0 - effective_eta * pi_E) * p13_0;
            p12_step = p12_0 + effective_eta * pi_E * p13_0;
            
            // Re-clamp after modification
            p12_step = clamp01(p12_step);
            p13_step = clamp01(p13_step);
        }

        // Ensure total <= 1
        const double pSum = p12_step + p13_step;
        if (pSum > 1.0) {
            p12_step /= pSum;
            p13_step /= pSum;
        }

        // Competing transitions: S3 priority then S2 (same as earlier)
        const double u = UniformRand();
        if (u < p13_step) {
            st.state = "S3";
            // Clean up checkpoint tracking if exists
            if (checkpointEnabled) {
                cellS2HoldAge_h.erase(cellID);
                cellT_hold_h.erase(cellID);
                cellCheckpointEngaged.erase(cellID);
            }
            // optionally reset age/duration here if your model expects it
        }
        else if (u < p13_step + p12_step) {
            st.state = "S2";
            // Probabilistic checkpoint: determine if checkpoint was engaged
            if (checkpointEnabled) {
                const double E_hold = E3 - k_hold * sigma;
                const double pi_E = CheckpointEngagementProbability(Ecur, E_hold, w_hold);
                
                // Sample checkpoint engagement
                const double u_checkpoint = UniformRand();
	                if (u_checkpoint < pi_E) {
	                    // Checkpoint engaged: sample hold duration and store
	                    // Use a cell-cycle-scaled lognormal hold distribution to keep arrest times
	                    // consistent across cell types with different cycle lengths.
	                    const double s_E = SaturatingLogMean(Ecur, E_hold, w_E_hold);
	                    const CycleInfo& cyc = cellCycleInfoMap.at(cellType);
	                    const double T_cycle = cyc.mTG1 + cyc.mTS + cyc.mTG2 + cyc.mTM; // hours
	                    const double T_cycle_safe =
	                        (std::isfinite(T_cycle) && T_cycle > 0.0) ? T_cycle : 1.0;

	                    const double median_min_h =
	                        std::max(1e-6, holdMedianMinCycleFraction * T_cycle_safe);
	                    const double median_max_h =
	                        std::max(median_min_h, holdMedianMaxCycleFraction * T_cycle_safe);

	                    const double mu_min = std::log(median_min_h);
	                    const double mu_max = std::log(median_max_h);
	                    const double mu_E = mu_min + (mu_max - mu_min) * s_E;

	                    double T_hold = SampleLogNormal(mu_E, sigmaT_hold);
	                    const double clamp_min_h =
	                        std::max(1e-6, holdClampMinCycleFraction * T_cycle_safe);
	                    const double clamp_max_h =
	                        std::max(clamp_min_h, holdClampMaxCycleFraction * T_cycle_safe);
	                    if (T_hold < clamp_min_h) T_hold = clamp_min_h;
	                    else if (T_hold > clamp_max_h) T_hold = clamp_max_h;
	                    
	                    cellT_hold_h[cellID] = T_hold;
	                    cellCheckpointEngaged[cellID] = true;
	                    cellS2HoldAge_h[cellID] = 0.0;  // Initialize hold age
                } else {
                    // No checkpoint: cell enters S2 normally
                    cellCheckpointEngaged[cellID] = false;
                }
            }
            // optionally reset age/duration here if your model expects it
        }
        else {
            st.age += increaseTime;
        }
        return;
    }

    // =========================
    // S2 transitions
    // =========================
    if (st.state == "S2")
    {
        const double Ecur = st.E;

        // S2 -> S1 (repair/recovery)
        double ES2ToS1 = -std::fabs(Ecur - E1) / (2.0 * sigma);
        double PS2ToS1 = 2.0 * GaussianCDF(ES2ToS1, 0, 1);

        // S2 -> S3 (progression to lethal)
        double ES2ToS3 = -std::fabs(Ecur - E3) / (2.0 * sigma);
        double PS2ToS3 = 2.0 * GaussianCDF(ES2ToS3, 0, 1);

        // Soft saturation: only trigger at extreme margin (Priority 1)
        // Hard saturation: Ecur >= E3, Soft saturation: Ecur >= E3 + margin*sigma
        if (softSaturationEnabled) {
            if (Ecur >= E3 + softSaturationMargin * sigma) PS2ToS3 = 1.0;
        } else {
            if (Ecur >= E3) PS2ToS3 = 1.0;  // Hard saturation
        }

        // Calculate baseline transition probabilities
        double p21_step = 0.0;
        double p23_0_step = 0.0;
        if (transitionType) {
            p21_step = InstantaneousStateJumpProb(PS2ToS1, increaseTime, To21);
            p23_0_step = InstantaneousStateJumpProb(PS2ToS3, increaseTime, To23);
        } else {
            p21_step = DelayedStateJumpProb(PS2ToS1, increaseTime, To21);
            p23_0_step = DelayedStateJumpProb(PS2ToS3, increaseTime, To23);
        }
        
        // Probabilistic checkpoint: apply gating function if in checkpoint
        double p23_step = p23_0_step;
        if (checkpointEnabled) {
            auto itEngaged = cellCheckpointEngaged.find(cellID);
            if (itEngaged != cellCheckpointEngaged.end() && itEngaged->second) {
                // Cell is in checkpoint hold
                auto itHoldAge = cellS2HoldAge_h.find(cellID);
                auto itT_hold = cellT_hold_h.find(cellID);
                if (itHoldAge != cellS2HoldAge_h.end() && itT_hold != cellT_hold_h.end()) {
                    const double t = itHoldAge->second;  // Current hold age
                    const double T_hold = itT_hold->second;  // Sampled hold duration
                    const double g = GatingFunction(t, T_hold, q_gate);
                    // Apply gating to catastrophe probability
                    p23_step = g * p23_0_step;
                }
            }
        }

        p21_step = clamp01(p21_step);
        p23_step = clamp01(p23_step);

        const double pSum = p21_step + p23_step;
        if (pSum > 1.0) {
            p21_step /= pSum;
            p23_step /= pSum;
        }

        const double u = UniformRand();
        if (u < p23_step) {
            st.state = "S3";
            // Clean up checkpoint tracking when cell dies (Priority 2)
            if (checkpointEnabled) {
                cellS2HoldAge_h.erase(cellID);
                cellT_hold_h.erase(cellID);
                cellCheckpointEngaged.erase(cellID);
            }
        }
        else if (u < p23_step + p21_step) {
            st.state = "S1";
            // Clean up checkpoint tracking when cell repairs (Priority 2)
            if (checkpointEnabled) {
                cellS2HoldAge_h.erase(cellID);
                cellT_hold_h.erase(cellID);
                cellCheckpointEngaged.erase(cellID);
            }
        }
        else {
            st.age += increaseTime;
        }
        return;
    }

    // =========================
    // S3 transitions (absorbing)
    // =========================
    if (st.state == "S3")
    {
        st.age += increaseTime;
        return;
    }

    // Unknown state string -> do nothing
    return;
}



void CellStateModel::CellStateTransitionMisrepair(int cellID,
                                        double E1, double E2, double E3,
                                        double sigma,
                                        double increaseTime,
                                        bool transitionType) // true=instantaneous, false=delayed
{
    // ---- 0) Guard: cell might have been deleted earlier this tick ----
    auto itState = cellStateMap.find(cellID);
    if (itState == cellStateMap.end())
        return;

    // Reference (alias) to the map-stored StateInfo (no copy)
    auto& st = itState->second;

    // Helper to keep probabilities sane
    auto clamp01 = [](double x) {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    };


    auto itType  = cellTypeMap.find(cellID);
    const std::string cellType = itType->second.GetCellType();
    const double To12 =cellStateParaInfoMap.at(cellType).To12;
    const double To13 = cellStateParaInfoMap.at(cellType).To13;
    const double To21 = cellStateParaInfoMap.at(cellType).To21;
    const double To23 = cellStateParaInfoMap.at(cellType).To23;
    const CellStateParaInfo& par = cellStateParaInfoMap.at(cellType);

    // =========================
    // S1 transitions
    // =========================
    if (st.state == "S1")
    {
        const double Ecur = st.E;

        // Overlap-based cumulative probabilities over observation window To
        double ES1ToS2 = -std::fabs(Ecur - E2) / (2.0 * sigma);
        double PS1ToS2 = 2.0 * GaussianCDF(ES1ToS2, 0, 1);

        double ES1ToS3 = -std::fabs(Ecur - E3) / (2.0 * sigma);
        double PS1ToS3 = 2.0 * GaussianCDF(ES1ToS3, 0, 1);

        // Soft saturation: only trigger at extreme margin (Priority 1)
        // Hard saturation: Ecur >= E3, Soft saturation: Ecur >= E3 + margin*sigma
        if (softSaturationEnabled) {
            if (Ecur >= E3 + softSaturationMargin * sigma) PS1ToS3 = 1.0;
        } else {
            if (Ecur >= E3) PS1ToS3 = 1.0;  // Hard saturation
        }

        // Map cumulative probability -> per-step probability
        double p12_step = 0.0;
        double p13_step = 0.0;
        if (transitionType) {
            p12_step = InstantaneousStateJumpProb(PS1ToS2, increaseTime, To12);
            p13_step = InstantaneousStateJumpProb(PS1ToS3, increaseTime, To13);
        } else {
            p12_step = DelayedStateJumpProb(PS1ToS2, increaseTime, To12);
            p13_step = DelayedStateJumpProb(PS1ToS3, increaseTime, To13);
        }

        p12_step = clamp01(p12_step);
        p13_step = clamp01(p13_step);

        // Probabilistic checkpoint: modify probabilities if enabled
        if (checkpointEnabled) {
            // Calculate checkpoint engagement probability
            const double E_hold = E3 - k_hold * sigma;
            const double pi_E = CheckpointEngagementProbability(Ecur, E_hold, w_hold);
            
            // Store baseline probabilities
            const double p12_0 = p12_step;
            const double p13_0 = p13_step;
            
            // Apply marginalization formulas
            // Use effective_eta = eta directly (default eta = 1.0)
            double effective_eta = eta;
            
            p13_step = (1.0 - effective_eta * pi_E) * p13_0;
            p12_step = p12_0 + effective_eta * pi_E * p13_0;
            
            // Re-clamp after modification
            p12_step = clamp01(p12_step);
            p13_step = clamp01(p13_step);
        }

        // Ensure total <= 1
        const double pSum = p12_step + p13_step;
        if (pSum > 1.0) {
            p12_step /= pSum;
            p13_step /= pSum;
        }

        // Competing transitions: S3 priority then S2 (same as earlier)
        const double u = UniformRand();
        if (u < p13_step) {
            st.state = "S3";
            // Clean up checkpoint tracking if exists
            cellS2HoldAge_h.erase(cellID);
            cellT_hold_h.erase(cellID);
            cellCheckpointEngaged.erase(cellID);
            // optionally reset age/duration here if your model expects it
        }
        else if (u < p13_step + p12_step) {
            st.state = "S2";
            // Probabilistic checkpoint: determine if checkpoint was engaged
            if (checkpointEnabled) {
                const double E_hold = E3 - k_hold * sigma;
                const double pi_E = CheckpointEngagementProbability(Ecur, E_hold, w_hold);
                
                // Sample checkpoint engagement
                const double u_checkpoint = UniformRand();
	                if (u_checkpoint < pi_E) {
	                    // Checkpoint engaged: sample hold duration and store
	                    // Use a cell-cycle-scaled lognormal hold distribution to keep arrest times
	                    // consistent across cell types with different cycle lengths.
	                    const double s_E = SaturatingLogMean(Ecur, E_hold, w_E_hold);
	                    const CycleInfo& cyc = cellCycleInfoMap.at(cellType);
	                    const double T_cycle = cyc.mTG1 + cyc.mTS + cyc.mTG2 + cyc.mTM; // hours
	                    const double T_cycle_safe =
	                        (std::isfinite(T_cycle) && T_cycle > 0.0) ? T_cycle : 1.0;

	                    const double median_min_h =
	                        std::max(1e-6, holdMedianMinCycleFraction * T_cycle_safe);
	                    const double median_max_h =
	                        std::max(median_min_h, holdMedianMaxCycleFraction * T_cycle_safe);

	                    const double mu_min = std::log(median_min_h);
	                    const double mu_max = std::log(median_max_h);
	                    const double mu_E = mu_min + (mu_max - mu_min) * s_E;

	                    double T_hold = SampleLogNormal(mu_E, sigmaT_hold);
	                    const double clamp_min_h =
	                        std::max(1e-6, holdClampMinCycleFraction * T_cycle_safe);
	                    const double clamp_max_h =
	                        std::max(clamp_min_h, holdClampMaxCycleFraction * T_cycle_safe);
	                    if (T_hold < clamp_min_h) T_hold = clamp_min_h;
	                    else if (T_hold > clamp_max_h) T_hold = clamp_max_h;
	                    
	                    cellT_hold_h[cellID] = T_hold;
	                    cellCheckpointEngaged[cellID] = true;
	                    cellS2HoldAge_h[cellID] = 0.0;  // Initialize hold age
                } else {
                    // No checkpoint: cell enters S2 normally
                    cellCheckpointEngaged[cellID] = false;
                }
            }
            // optionally reset age/duration here if your model expects it
        }
        else {
            st.age += increaseTime;
        }
        return;
    }

    // =========================
    // S2 transitions (with misrepair channel)
    // =========================
    if (st.state == "S2")
    {
        const double Ecur = st.E;

        double ES2ToS1 = -std::fabs(Ecur - E1) / (2.0 * sigma);
        double PS2ToS1 = 2.0 * GaussianCDF(ES2ToS1, 0, 1);

        double ES2ToS3 = -std::fabs(Ecur - E3) / (2.0 * sigma);
        double PS2ToS3 = 2.0 * GaussianCDF(ES2ToS3, 0, 1);
        // Soft saturation: only trigger at extreme margin (Priority 1)
        // Hard saturation: Ecur >= E3, Soft saturation: Ecur >= E3 + margin*sigma
        if (softSaturationEnabled) {
            if (Ecur >= E3 + softSaturationMargin * sigma) PS2ToS3 = 1.0;
        } else {
            if (Ecur >= E3) PS2ToS3 = 1.0;  // Hard saturation
        }

        double p21_step = 0.0;
        double p23_energy_0_step = 0.0;
        if (transitionType) {
            p21_step        = InstantaneousStateJumpProb(PS2ToS1, increaseTime, To21);
            p23_energy_0_step = InstantaneousStateJumpProb(PS2ToS3, increaseTime, To23);
        } else {
            p21_step        = DelayedStateJumpProb(PS2ToS1, increaseTime, To21);
            p23_energy_0_step = DelayedStateJumpProb(PS2ToS3, increaseTime, To23);
        }

        p21_step        = clamp01(p21_step);
        p23_energy_0_step = clamp01(p23_energy_0_step);

        // Probabilistic checkpoint: apply gating function to catastrophe probability
        double p23_cat_step = p23_energy_0_step;
        if (checkpointEnabled) {
            auto itEngagedMis = cellCheckpointEngaged.find(cellID);
            if (itEngagedMis != cellCheckpointEngaged.end() && itEngagedMis->second) {
                // Cell is in checkpoint hold
                auto itHoldAgeMis = cellS2HoldAge_h.find(cellID);
                auto itT_holdMis = cellT_hold_h.find(cellID);
                if (itHoldAgeMis != cellS2HoldAge_h.end() && itT_holdMis != cellT_hold_h.end()) {
                    const double t = itHoldAgeMis->second;  // Current hold age
                    const double T_hold = itT_holdMis->second;  // Sampled hold duration
                    const double g = GatingFunction(t, T_hold, q_gate);
                    // Apply gating to catastrophe probability
                    p23_cat_step = g * p23_energy_0_step;
                }
            }
        }

        // =========================
        // MISREPAIR ON RECOVERY ATTEMPT (Version 2)
        // 
        // Radiobiological interpretation:
        //   - Misrepair represents mitotic catastrophe: cells leave arrest,
        //     try to divide with misrepaired damage, then fail.
        //   - Misrepair death is triggered when cell ATTEMPTS recovery (S2→S1),
        //     not as a constant background hazard.
        //
        // Logic:
        //   p_mis = probability of lethal misrepair given a recovery attempt
        //   p23_mis = p21 × p_mis      (death via failed recovery)
        //   p21_eff = p21 × (1-p_mis)  (successful recovery)
        //   p23 = p23_cat + p23_mis    (total death: catastrophe + misrepair)
        //
        // This produces LQ-like behavior:
        //   - Misrepair on recovery → linear αD component (low dose)
        //   - Catastrophe → quadratic βD² component (high dose)
        // =========================
        const double k_error_dsb = par.k_error; // per (DSB*hour)
        const double alpha_dsb   = par.alpha;   // energy per DSB

        const double k_error_energy =
            (alpha_dsb > 0.0) ? (std::max(0.0, k_error_dsb) / alpha_dsb) : 0.0; // 1/(energy*hour)

        // Misrepair probability (conditional on recovery attempt)
        const double lambda_mis = k_error_energy * std::max(0.0, Ecur); // 1/h
        double p_mis = 1.0 - std::exp(-lambda_mis * std::max(0.0, increaseTime));
        p_mis = clamp01(p_mis);

        // Misrepair death only occurs when cell attempts recovery
        double p23_mis_step = p21_step * p_mis;           // Death via failed recovery
        double p21_eff_step = p21_step * (1.0 - p_mis);   // Successful recovery

        // Total death = catastrophe + misrepair (mutually exclusive paths)
        double p23_step = p23_cat_step + p23_mis_step;
        p23_step = clamp01(p23_step);

        // Update p21 to effective recovery probability
        p21_step = p21_eff_step;

        // Renormalize if needed (should not exceed 1 with this logic, but safety check)
        const double pSum = p21_step + p23_step;
        if (pSum > 1.0) {
            p21_step /= pSum;
            p23_step /= pSum;
        }

        const double u = UniformRand();
        if (u < p23_step) {
            st.state = "S3";
            // Clean up checkpoint tracking when cell dies (Priority 2)
            if (checkpointEnabled) {
                cellS2HoldAge_h.erase(cellID);
                cellT_hold_h.erase(cellID);
                cellCheckpointEngaged.erase(cellID);
            }
        }
        else if (u < p23_step + p21_step) {
            st.state = "S1";
            // Clean up checkpoint tracking when cell repairs (Priority 2)
            if (checkpointEnabled) {
                cellS2HoldAge_h.erase(cellID);
                cellT_hold_h.erase(cellID);
                cellCheckpointEngaged.erase(cellID);
            }
        }
        else {
            st.age += increaseTime;
        }
        return;
    }

    // =========================
    // S3 transitions (absorbing)
    // =========================
    if (st.state == "S3")
    {
        st.age += increaseTime;
        return;
    }

    // Unknown state string -> do nothing
    return;
}


double CellStateModel::InstantaneousStateJumpProb(double p_sp, double increaseTime, double ObservationTime)
{
    //     double p = (1-pow((1-p_sp),increaseTime/ObservationTime));// probability in one time step
    double p = p_sp;
    return p;

}
double CellStateModel::DelayedStateJumpProb(double p_sp, double increaseTime, double ObservationTime)
{
        double p = (1-pow((1-p_sp),increaseTime/ObservationTime));// probability in one time step
//     double p = p_sp;
    return p;

}


map< int, string > CellStateModel::GetCellPhase()
{
    std::map<int, std::string> phaseMap;
    for (std::map<int,PhaseInfo>::iterator mitr_cell=cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
    {
        phaseMap[mitr_cell->first] = mitr_cell->second.phase;
    }
    return phaseMap;

}

map< int, double > CellStateModel::GetCellPhaseDuration()
{
    std::map<int, double> durationMap;
    for (std::map<int,PhaseInfo>::iterator mitr_cell=cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
    {
        durationMap[mitr_cell->first] = mitr_cell->second.duration;
    }
    return durationMap;

}
map< int, double > CellStateModel::GetCellAge()
{
    std::map<int, double> ageMap;
    for (std::map<int,PhaseInfo>::iterator mitr_cell=cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
    {
        ageMap[mitr_cell->first] = mitr_cell->second.age;
    }
    return ageMap;

}

map< int, string > CellStateModel::GetCellState()
{
    std::map<int, std::string> stateMap;
    for(std::map<int, StateInfo>::iterator mitr_cell=cellStateMap.begin();mitr_cell!=cellStateMap.end();mitr_cell++)
    {
        stateMap[mitr_cell->first] = mitr_cell->second.state;
    }

    return stateMap;
}
map< int, double > CellStateModel::GetCellStateDuration()
{
    std::map<int, double> durationMap;
    for(std::map<int, StateInfo>::iterator mitr_cell=cellStateMap.begin();mitr_cell!=cellStateMap.end();mitr_cell++)
    {
        durationMap[mitr_cell->first] = mitr_cell->second.duration;
    }    
    return durationMap;
}

map< int, double > CellStateModel::GetCellStateAge()
{
    std::map<int, double> ageMap;
    for(std::map<int, StateInfo>::iterator mitr_cell=cellStateMap.begin();mitr_cell!=cellStateMap.end();mitr_cell++)
    {
        ageMap[mitr_cell->first] = mitr_cell->second.age;
    }     
    return ageMap;
}

map< int, double > CellStateModel::GetCellPositionX()
{
    std::map<int, double> XPositionMap;
    for (std::map<int, PositionInfo>::iterator mitr_cell = cellPositionMap.begin();mitr_cell!=cellPositionMap.end();mitr_cell++)
    {
//         XPositionMap[mitr_cell->first] = d/2+mitr_cell->second.i*d;
        XPositionMap[mitr_cell->first] = mitr_cell->second.i*d;// when take the grid points as center of cell
    }
    return XPositionMap;

}
map< int, double > CellStateModel::GetCellPositionY()
{
    std::map<int, double> YPositionMap;
    for (std::map<int, PositionInfo>::iterator mitr_cell = cellPositionMap.begin();mitr_cell!=cellPositionMap.end();mitr_cell++)
    {
//         YPositionMap[mitr_cell->first] = d/2+mitr_cell->second.j*d;
        YPositionMap[mitr_cell->first] = mitr_cell->second.j*d;//when take the grid points as center of cell
    }
    return YPositionMap;
}

map< int, double > CellStateModel::GetCellPositionZ()
{
    std::map<int, double> ZPositionMap;
    for (std::map<int, PositionInfo>::iterator mitr_cell = cellPositionMap.begin();mitr_cell!=cellPositionMap.end();mitr_cell++)
    {
//         ZPositionMap[mitr_cell->first] = d/2 + mitr_cell->second.k*d;
        ZPositionMap[mitr_cell->first] =  mitr_cell->second.k*d;//when take the grid points as center of cell
    }
    return ZPositionMap;
}

// Priority 2: Checkpoint hold diagnostics
std::map<int, double> CellStateModel::GetCheckpointHoldAge()
{
    return cellS2HoldAge_h;  // Return copy of checkpoint hold age map
}

int CellStateModel::GetCheckpointHoldCount()
{
    return static_cast<int>(cellS2HoldAge_h.size());  // Number of cells in checkpoint hold
}

void CellStateModel::SetCheckpointEnabled(bool enabled)
{
    checkpointEnabled = enabled;
    if (!enabled) {
        // Clear checkpoint tracking when disabling
        cellS2HoldAge_h.clear();
        cellT_hold_h.clear();
        cellCheckpointEngaged.clear();
    }
}

bool CellStateModel::IsCheckpointEnabled() const
{
    return checkpointEnabled;
}

void CellStateModel::SetSoftSaturationEnabled(bool enabled)
{
    softSaturationEnabled = enabled;
}

bool CellStateModel::IsSoftSaturationEnabled() const
{
    return softSaturationEnabled;
}

map< int, int > CellStateModel::GetCellAncestryID()
{
    std::map<int, int> cellAncestryIDMap;
    for (std::map<int, PhaseInfo>::iterator mitr_cell = cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
    {
        cellAncestryIDMap[mitr_cell->first] = mitr_cell->second.ancestryID;
    }
    return cellAncestryIDMap;

}

void CellStateModel::SetUpContactInhibition(bool considerOrNot)
{
    considerContactInhibition = considerOrNot;

}

CellStateModelParameter CellStateModel::GetCellStateModelParameters(double p_sp, double p_mis_sp, double avgN)
{
    // this function is used to generate the parameters for cell state model
  
  // the first argument is the probability of spontaneous death
  // the second argument is the probability of misrepair 
  
  
  StandardNormalDistribution snd;
  
  double E1 = 0;
  
  double x = snd.NormalCDFInverse(p_sp/2);
  double E3 = 4*x*x;
  double sigma = sqrt(E3);
  
//   cout<<"E3 "<<E3<<endl;
//   cout<<"sigma "<<sigma<<endl;
  double p_arr_sp = p_sp/p_mis_sp;
  cout<<"p_arr_sp "<<p_arr_sp<<endl;
  double x1 = snd.NormalCDFInverse(p_arr_sp/2);
  cout<<"x1 = "<<x1<<endl;
  double E2 = -2*sigma*x1;
  
//   cout<<"E2 "<<E2<<endl;

  
  double alpha = E3/avgN;
  
//   cout<<"alpha "<<alpha<<endl;
  
  
  // for calculating beta, we use an approximate method
  // we use the fitting parameters of the equations governing the bystander signal diffusion
  
  std::vector<double> avC;
  double min_avC = 0.1;
  double max_avC = 62.155;
  double deltaC = (max_avC - min_avC)/200;
  
  for (int i = 0; i<200; i++)
  {
    avC.push_back(min_avC + (i)*deltaC);
  }
  
  
//   for (int i = 0; i<avC.size(); i++)
//   {
//     cout<<"avC "<<i<<" is "<<avC[i]<<endl;
//   }
//   
  
  std::vector<double> bys_f;
  
  for (int i = 0; i<avC.size(); i++)
  {
    double temp = pow(log(62.155/avC[i]),0.1228);
    bys_f.push_back(1- 0.6844*exp(0.282*temp));
  }
  
//   for (int i = 0; i<bys_f.size(); i++)
//   {
//     cout<<"bys_f "<<i <<" is "<<bys_f[i]<<endl;
//   }

  
  std::vector<double>betaVec;
  
  for (int i = 0; i<avC.size(); i++)
  {
    double x = snd.NormalCDFInverse(bys_f[i]/2);
    betaVec.push_back((E3 + 2*sigma*x)/avC[i]);
  }
  
//   ofstream file("checkbeta.csv");
  
//   for (int i = 0; i<betaVec.size(); i++)
//   {
//     cout<<"betaVec i "<<i<<" is "<<betaVec[i]<<endl;
//     file<<avC[i]<<","<< betaVec[i]<<endl;
//   }
//   
  
  
  
//     exit(0);
  double sum = 0;
  
  for (int i = 0; i<betaVec.size(); i++)
  {
    sum = sum + betaVec[i];
  }
  
  double beta = sum/betaVec.size();
  
  CellStateModelParameter cpr;
  cpr.alpha = alpha;
  cpr.beta = beta;
  cpr.E1 = E1;
  cpr.E2 = E2;
  cpr.E3 = E3;
  cpr.sigma = sigma;
  
  return cpr;
  

}
