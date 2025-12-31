// 
// This code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015
#ifndef Cell_h
#define Cell_h 1

#include <string>
#include <vector>

struct CellCycleParameter
  {
      double mTG1; // mean time for staying at G1
      double sigmaTG1; // sigma for G1
      double mTS; // mean time for staying at S
      double sigmaTS; // sigma for S
      double mTG2; // mean time for staying at G2
      double sigmaTG2; // sigma for G2
      double mTM; // mean time for staying at M
      double sigmaTM; // sigma for M
      double fG1 ; // relative radiation sensitivity for phase G1
      double fG2 ; // relative radiation sensitivity for phase G2
      double fM ; // relative radiation sensitivity for phase M
      double fS; // relative radiation sensitivity for phase S
  };
  struct CellStateParameter
  {
    double E1;// mean state energy for S1
    double E2;// mean state enery for S2
    double E3;// mean state energy for S3
    double sigma; // signa for gaussian distribution of state energy
    double alpha;//constant alpha for calcualting Delta E
    double beta; //  constant beta for calcualting Delta E
    double To12; // observation window width for state transition: S1 to S2
    double To13; // observation window width for state transition: S1 to S3
    double To21; // observation window width for state transition: S2 to S1
    double To23;// observation window width for state transtiion: S2 to S3
  };

 struct CellDNADamageRepairParameter
 {
    double f1;// faction for fast repair using two exponential distribution for DNA damage repair
    double f2; // faction for slow repair
    double lambda1; // unit in 1/hour
    double lambda2; // unit in 1/hour
 };
struct CellBystanderSignalParameter
{
    int Rt;// Receptor threshold, unit in #
    double mu; //Secretion rate unit in #/s
};


class Cell // creating a class named cell
{
public:
//   Cell(); // constructor
  void CellConstruct(std::string cType, std::string cOrganelle, std::string cShape, std::string cColor);
  std::string GetCellType(){return cellType;}
  std::string GetCellShape(){return cellShape;}
  
  std::vector<std::string> & GetCellOrganelle () {return cellOrganelle;}
  std::vector<std::string> & GetCellOrganelleColor() {return colorOfOrganelle;}
  void SetCellParametersForCellStateModeling(CellCycleParameter cyclePara, CellStateParameter statePara, CellDNADamageRepairParameter dnaRepairPara, CellBystanderSignalParameter bystanderSignalPara);
  CellCycleParameter GetCellCycleParameters() {return theCellCyleParameter;}
  CellStateParameter GetCellStateParameters() {return theCellStateParameter;}
  CellDNADamageRepairParameter GetCellDNADamageRepairParameters(){return theCellDNADamageRepairParameter;}
  CellBystanderSignalParameter GetCellBystanderSignalParameters(){return theCellBystanderSignalParameter;}

private:
  std::string cellType;// cell type
  std::vector<std::string> cellOrganelle; // a vector storing cell organelle
  std::vector<std::string> colorOfOrganelle;// a vector storing cell
  std::string cellShape; // cell shape

  CellCycleParameter theCellCyleParameter;
  CellStateParameter theCellStateParameter;
  CellDNADamageRepairParameter theCellDNADamageRepairParameter;
  CellBystanderSignalParameter theCellBystanderSignalParameter;

};
#endif
