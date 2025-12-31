import matplotlib.pyplot as plt
import numpy as np
import  math

import os
import sys



fileName = 'Mono_Electron_Plane_Eng_1PNum_1000000_cellSurvival_'

path = os.getcwd()
fileNameList = []
for i in os.listdir(path):
    if os.path.isfile(os.path.join(path,i)) and fileName in i:
        fileNameList.append(i)
        
print fileNameList


def getCellSF(fileName):
    data = np.genfromtxt(fileName, delimiter=',',skip_header=0)
    print type(data)
    print data.shape
    dose = data[:,0]
    totalCellNum=data[:,1]
    SF=data[:,2]
    return dose, SF # return the cell Survival fraction of cells

sfFileList = fileNameList

def sfErrorAnalysis(sfFileList):
  dose, sf1 = getCellSF(sfFileList[0])
  sf = np.zeros(len(sf1)) 
  sigma = np.zeros(len(sf1))
  
  for i in range(len(sfFileList)):
    dose, sf_temp = getCellSF(sfFileList[i])
    for j in range(len(sf)):
      sf[j] = sf[j] + sf_temp[j]
      
 
  for j in range(len(sf)):
    sf[j] = sf[j]/len(sfFileList)
    
  for i in range(len(sfFileList)):
    dose, sf_temp = getCellSF(sfFileList[i])
    for j in range(len(sf)):
      sigma[j] = sigma[j] + (sf_temp[j] - sf[j])**2
  
  for j in range(len(sf)):
    sigma[j] = np.sqrt(sigma[j]/len(sfFileList))
  
  
  return sf, sigma
    
#def sfErrorAnalysis(sf1,sf2,sf3):
    #print  'size',len(sf1)
    #sf =np.zeros(len(sf1))
    #sigma =np.zeros(len(sf1))
    #for i in range(0,len(sf1)):
        #sf[i] = 1.0/3*(sf1[i]+sf2[i]+sf3[i])
        #sigma[i] = np.sqrt(1.0/3*((sf1[i]-sf[i])**2+(sf2[i]-sf[i])**2+(sf3[i]-sf[i])**2))
    
        #print sf[i], sigma[i]
    #return  sf, sigma




def plotSFCurve(dose,sf,sfError,color,label):  
    print 'the sf is ',sf,'the size of sf', len(sf)
    print 'the error is',sfError
    print 'the dose is ',dose,'the size of dose',len(dose)
    plt.errorbar(dose,sf,sfError,fmt = color+'.-',markersize=10, linewidth=2.0,label = label)
    plt.yscale('log')   
    
    
dose1, sf1 = getCellSF(sfFileList[0])

sf, sigma =  sfErrorAnalysis(sfFileList)

print sf
print sigma

#sys.exit(0)

plotSFCurve(dose1,sf,sigma,'b','sf considering RIBE')
#plt.plot(dose3,sf3)
#plt.yticks(np.arange(0,1.1,0.05))  
plt.rcParams['figure.figsize'] = (9.33, 7) # set the default plot size
plt.hold(True)

plt.xlim([0,8.5])
#plt.ylim([1E-6,1])
plt.tick_params(axis='x', labelsize=15)
plt.tick_params(axis='y', labelsize=15)
plt.xlabel('Dose, Gy',fontsize=15)
plt.ylabel('Cell Survival Fraction',fontsize=15)
plt.grid(True,which="both",ls="-")
plt.legend(loc='best')
plt.savefig('SF_RIBE_Erotron_Monolayer.svg')

plt.show()