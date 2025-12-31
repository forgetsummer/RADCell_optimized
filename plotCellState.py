import matplotlib.pyplot as plt
import numpy as np
import  math

from pylab import rcParams

import os
import sys


fileName = 'Mono_Electron_Plane_Eng_1PNum_1000000_cellState_dose_'

path = os.getcwd()
fileNameList = []
for i in os.listdir(path):
    if os.path.isfile(os.path.join(path,i)) and fileName in i:
        fileNameList.append(i)
        
print(fileNameList)

#sys.exit(0)
        
def PlotCellState(fileName):
  print(fileName)
  data = np.genfromtxt(fileName, delimiter=',',skip_header=0)
  
  plt.figure(figsize=(9.33, 7))
  
  print(type(data))
  print(data.shape)
  time = data[:,0]
  S1Num=data[:,1]
  S21Num=data[:,2]
  S22Num=data[:,3]
  S3Num = data[:,4]
  num_total =S1Num + S21Num + S22Num + S3Num
  S1Ratio=S1Num/num_total
  S21Ratio=S21Num/num_total
  S22Ratio=S22Num/num_total
  S3Ratio= S3Num/num_total

  S2Ratio = S21Ratio+ S22Ratio
  print('the ratio is ', S1Ratio, S21Ratio,S22Ratio,S3Ratio)

  plt.ylim(-0.1,1.1)
  plt.plot(time,S1Ratio,'b.',label='S1 State Cell')
  plt.plot(time,S2Ratio,'g.',label='S2 State Cell')
  plt.plot(time,S3Ratio,'r.',label='S3 State Cell')
  plt.tick_params(axis='x', labelsize=15)
  plt.tick_params(axis='y', labelsize=15)
  plt.xlabel('Time,minutes',fontsize=15)
  plt.ylabel('Cell State Ratio',fontsize=15)
  plt.tick_params(axis='x', labelsize=15)
  plt.tick_params(axis='y', labelsize=15)
  plt.legend(loc='best')
  plt.savefig(fileName[20:-4] + '.svg')

for i in range(len(fileNameList)):
  PlotCellState(fileNameList[i])
  
  
plt.show()
