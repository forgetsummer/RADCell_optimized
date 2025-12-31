#include "CellStateModel.hh"
#include <stdlib.h>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "StandardNormalDistribution.hh"

// Static random number generator for uniform distribution [0,1)
static std::mt19937& GetRandomEngine() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

static double UniformRand() {
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(GetRandomEngine());
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

    cellCycleInfoMap[cType].mTG1 =theCell.GetCellCycleParameters().mTG1;
    cellCycleInfoMap[cType].sigmaTG1 = theCell.GetCellCycleParameters().sigmaTG1;
    cellCycleInfoMap[cType].mTS = theCell.GetCellCycleParameters().mTS;
    cellCycleInfoMap[cType].sigmaTS = theCell.GetCellCycleParameters().sigmaTS;
    cellCycleInfoMap[cType].mTG2 = theCell.GetCellCycleParameters().mTG2;
    cellCycleInfoMap[cType].sigmaTG2 = theCell.GetCellCycleParameters().sigmaTG2;
    cellCycleInfoMap[cType].mTM = theCell.GetCellCycleParameters().mTM;
    cellCycleInfoMap[cType].sigmaTM = theCell.GetCellCycleParameters().sigmaTM;
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
    const bool proliferativeFromCellState = (itState->second.state == "S1");

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

    // ---- 2) "Space allowable" check at the CURRENT location ----
    // Note: Here we only decide whether cell should remain in cycling vs go to G0.
    // We'll recompute candidate positions again at actual division time.
    bool allowableInSpace = true;
    if (considerContactInhibition)
    {
        auto possibleNow = CheckCellContactInhibitionCondition(pos.i, pos.j, pos.k);
        allowableInSpace = !possibleNow.empty();
    }
    else
    {
        allowableInSpace = true; // by your design: ignore space constraint
    }

    // If environment not tolerable, go to G0
    if (!allowableInCellTelomereLength || !allowableInSpace)
    {
        phase.phase = "G0";
        phase.age = 0;
        phase.duration = 1E+18;
        return;
    }

    // ---- 3) Phase transitions ----
    if (phase.phase == "G0")
    {
        phase.phase = "G1";
        phase.age = 0;
        phase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTG1,
                                          cellCycleInfoMap.at(cellType).sigmaTG1);
        return;
    }

    if (phase.phase == "G1")
    {
        phase.age += increaseTime;
        if (phase.age > phase.duration)
        {
            phase.phase = "S";
            phase.age = 0;
            phase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTS,
                                              cellCycleInfoMap.at(cellType).sigmaTS);
        }
        return;
    }

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

    if (phase.phase == "G2")
    {
        phase.age += increaseTime;
        if (phase.age > phase.duration)
        {
            phase.phase = "M";
            phase.age = 0;
            phase.duration = GaussianSampling(cellCycleInfoMap.at(cellType).mTM,
                                              cellCycleInfoMap.at(cellType).sigmaTM);
        }
        return;
    }

    if (phase.phase == "M")
    {
        phase.age += increaseTime;

        if (phase.age <= phase.duration)
            return;

        // ---- 4) Division happens now ----
        // Recompute possible positions NOW (not earlier) to avoid stale neighbor list.
        std::map<std::string, PositionInfo> possibleCellHomeForMitosis;
        if (considerContactInhibition)
        {
            possibleCellHomeForMitosis = CheckCellContactInhibitionCondition(pos.i, pos.j, pos.k);
        }
        else
        {
            // Your current design: ignore space => allow "original".
            // NOTE: This allows both daughters in same voxel unless you later resolve overlaps.
            possibleCellHomeForMitosis.emplace("original", pos);
        }

        // Defensive: if CI enabled but no space now, push to G0 instead of crashing / undefined.
        if (possibleCellHomeForMitosis.empty())
        {
            phase.phase = "G0";
            phase.age = 0;
            phase.duration = 1E+18;
            return;
        }

        // Find max existing cellID (you can later optimize with nextCellID)
        int cellID_max = 0;
        for (auto it = cellPhaseMap.begin(); it != cellPhaseMap.end(); ++it)
            cellID_max = std::max(cellID_max, it->first);

        const int id1 = cellID_max + 1;
        const int id2 = cellID_max + 2;

        // Cache mother info before erase
        const double motherTel = phase.telomereFraction;
        const int motherAncestry = phase.ancestryID;
        const PositionInfo motherPos = pos;

        // Daughter 1
        cellPhaseMap[id1].phase = "G0";
        cellPhaseMap[id1].age = 0;
        cellPhaseMap[id1].duration = 1E+18;
        cellPhaseMap[id1].telomereFraction = 0.5 * motherTel;
        cellPhaseMap[id1].ancestryID = motherAncestry;

        cellStateMap[id1].state = "S1";
        cellStateMap[id1].age = 0;
        cellStateMap[id1].duration = 1E+18;

        cellPositionMap[id1] = motherPos; // first daughter takes mother's place
        cellHomeResidence[motherPos.i][motherPos.j][motherPos.k] = 1; // keep occupied

        // Daughter 2 (random pick)
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

        cellPositionMap[id2] = pick->second;
        cellHomeResidence[pick->second.i][pick->second.j][pick->second.k] = 1;

        // Delete mother
        cellPhaseMap.erase(cellID);
        cellStateMap.erase(cellID);
        cellPositionMap.erase(cellID);
        // (Optionally erase proliferative flag too, if you maintain it per-cell)
        cellProliferativeMap.erase(cellID);

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

    // instantaneous vs delayed transition type
    const bool transitionType = (DSBNum > 0);

    // baseline state energies
    const double E1    = cellStateParaInfoMap.at(cellType).E1;
    const double E2    = cellStateParaInfoMap.at(cellType).E2;
    const double E3    = cellStateParaInfoMap.at(cellType).E3;
    const double sigma = cellStateParaInfoMap.at(cellType).sigma;

    // current phase string (safe via iterator)
    const std::string& ph = itPhase->second.phase;

    // ---- 1) Dispatch by phase (keep your original formulas) ----
    if (ph == "G1")
    {
        CellStateTransition(cellID, E1, E2, E3, sigma, increaseTime, transitionType);
    }
    else if (ph == "S")
    {
        const double fS = cellStateParaInfoMap.at(cellType).fS;

        const double E1_S = 0.0;
        const double E3_S = std::sqrt(E3 * E3 - 8.0 * sigma * sigma * std::log(fS));
        const double E2_S = E3_S - std::sqrt((E3 - E2) * (E3 - E2) - 8.0 * sigma * sigma * std::log(fS));

        CellStateTransition(cellID, E1_S, E2_S, E3_S, sigma, increaseTime, transitionType);
    }
    else if (ph == "G2")
    {
        const double fG2 = cellStateParaInfoMap.at(cellType).fG2;

        const double E1_G2 = 0.0;
        const double E3_G2 = std::sqrt(E3 * E3 - 8.0 * sigma * sigma * std::log(fG2));
        const double E2_G2 = E3_G2 - std::sqrt((E3 - E2) * (E3 - E2) - 8.0 * sigma * sigma * std::log(fG2));

        CellStateTransition(cellID, E1_G2, E2_G2, E3_G2, sigma, increaseTime, transitionType);
    }
    else if (ph == "M")
    {
        const double fM = cellStateParaInfoMap.at(cellType).fM;

        const double E1_M = 0.0;
        const double E3_M = std::sqrt(E3 * E3 - 8.0 * sigma * sigma * std::log(fM));
        const double E2_M = E3_M - std::sqrt((E3 - E2) * (E3 - E2) - 8.0 * sigma * sigma * std::log(fM));

        CellStateTransition(cellID, E1_M, E2_M, E3_M, sigma, increaseTime, transitionType);
    }
    else if (ph == "G0")
    {
        // keep your original "treat G0 like G1" idea
        CellStateTransition(cellID, 0.0, E2, E3, sigma, increaseTime, transitionType);
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

        // Saturation (your design)
        if (Ecur >= E3) PS1ToS3 = 1.0;

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
            // optionally reset age/duration here if your model expects it
        }
        else if (u < p13_step + p12_step) {
            st.state = "S2";
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

        if (Ecur >= E3) PS2ToS3 = 1.0;

        double p21_step = 0.0;
        double p23_step = 0.0;
        if (transitionType) {
            p21_step = InstantaneousStateJumpProb(PS2ToS1, increaseTime, To21);
            p23_step = InstantaneousStateJumpProb(PS2ToS3, increaseTime, To23);
        } else {
            p21_step = DelayedStateJumpProb(PS2ToS1, increaseTime, To21);
            p23_step = DelayedStateJumpProb(PS2ToS3, increaseTime, To23);
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
        }
        else if (u < p23_step + p21_step) {
            st.state = "S1";
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


