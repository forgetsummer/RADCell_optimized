// 
// This is code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 

#include "Cell.hh"
#include <sstream>
#include <string>

using namespace std;


void Cell::CellConstruct(string cType, string cOrganelle, string cShape, string cColor)
{
  cellType=cType;
  cellShape=cShape;
  stringstream organList(cOrganelle);
  stringstream colorList(cColor);
  string organ;
  string color;
  while(organList>>organ)
  {
    cellOrganelle.push_back(organ);
  }
  while (colorList>>color)
  {
    colorOfOrganelle.push_back(color);
  }
  
}

void Cell::SetCellParametersForCellStateModeling(CellCycleParameter cyclePara, CellStateParameter statePara, CellDNADamageRepairParameter dnaRepairPara, CellBystanderSignalParameter bystanderSignalPara)
{
    theCellCyleParameter = cyclePara;
    theCellStateParameter = statePara;
    theCellDNADamageRepairParameter = dnaRepairPara;
    theCellBystanderSignalParameter = bystanderSignalPara;
}
