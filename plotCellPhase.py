import matplotlib.pyplot as plt
import numpy as np
import  math
from pylab import rcParams
import os
import sys

fileName = 'Mono_Electron_Plane_Eng_1PNum_1000000_cellPhase_dose_'

path = os.getcwd()
fileNameList = []
for i in os.listdir(path):
  if os.path.isfile(os.path.join(path,i)) and fileName in i:
    fileNameList.append(i)

print(fileNameList)

#sys.exit(0)

def func_cellPhase(data):
     
  print(type(data))
  print(data.shape)
  time = data[:,0]
  G0Num = data[:,1]
  G1Num=data[:,2]
  SNum=data[:,3]
  G2Num=data[:,4]
  num_mitosis = data[:,5]
  num_total = data[:,6]
  G0Ratio = G0Num/num_total
  G1Ratio=G1Num/num_total
  SRatio=SNum/num_total
  G2Ratio=G2Num/num_total
  mitosisRatio= num_mitosis/num_total
  print(mitosisRatio)
     
  return (time, G0Ratio,G1Ratio,SRatio,G2Ratio,mitosisRatio) 

def PlotCellPhaseInfo(fileName):
  data = np.genfromtxt(fileName, delimiter=',',skip_header=0)

  plt.figure(figsize=(9.33, 7))

  
  plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[1],'m.',label='G0 Phase Cell')
  plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[2],'k.',label='G1 Phase Cell')
  plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[3],'b.',label='S Phase Cell')
  plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[4],'g.',label='G2 Phase Cell')   
  plt.plot(func_cellPhase(data)[0],func_cellPhase(data)[5],'r.',label='M Phase Cell')
 
 
  plt.tick_params(axis='x', labelsize=15)
  plt.tick_params(axis='y', labelsize=15)
  plt.xlabel('Time, minutes',fontsize=15)
  plt.ylabel('Cell Phase Ratio',fontsize=15)
  plt.legend(bbox_to_anchor=(0.2, 1.15), loc='upper left', ncol=1)
  #plt.legend(loc='upper right')
  plt.tick_params(axis='x', labelsize=15)
  plt.tick_params(axis='y', labelsize=15)
  plt.xlabel('Time, minutes',fontsize=15)
  plt.ylabel('Cell Phase Ratio',fontsize=15)
  plt.legend(loc='best')
 

  plt.savefig(fileName[20:-4] + '.svg')

 
 
for i in range(len(fileNameList)):
  PlotCellPhaseInfo(fileNameList[i])
  
  
plt.show()
