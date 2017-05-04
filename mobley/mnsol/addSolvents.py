import os
from os import listdir
from os.path import join, isfile
import pandas as pd
nameThisFile = 'FullData.csv'
home = os.path.expanduser('~')
path = home+'/repos/testasymmetry/mobley/mnsol/'
os.chdir(path)

files = [f for f in listdir(path) if isfile(join(path,f))]

driverFile = pd.read_csv(os.path.join(path,'New_data.csv'))

dontWant = []
for f in files:
    if f.split('.')[-1] != 'csv':
        dontWant.append(f)
dontWant.append('New_data.csv')
dontWant.append('mobley_sa.csv')
dontWant.append(nameThisFile)

solvents = list(driverFile)
row1 = solvents[0]
solutes = list(driverFile[row1])
del solvents[0]

addInfo = {}
for f in files:
    solventInfo = {}
    if f.split('.')[0] in solvents:
        # No Ions in New_Data
        File = open(f,'r')
        for line in File:
            solventInfo[line.split(',')[0]] = line.split(',')[1].strip()
        File.close()
        addInfo[f.split('.')[0]] = solventInfo
    elif f.split('_')[0] in solvents:
        # Ions in New_Data
        File = open(f,'r')
        for line in File:
            solventInfo[line.split(',')[0]] = line.split(',')[1].strip()
        File.close()
        addInfo[f.split('_')[0]] = solventInfo
    elif (f.split('.')[0] not in solvents and f.split('.')[0] not in dontWant and f not in dontWant):
        # No Ions not in New_Data
        File = open(f,'r')
        for line in File:
            try:
                solventInfo[line.split(',')[0]] = line.split(',')[1].strip()
            except:
                pass
        File.close()
        addInfo[f.split('.')[0]] = solventInfo
    elif (f.split('_')[0] not in solvents and f.split('_')[0] not in dontWant and f not in dontWant):
        # Ions not in New_Data
        File = open(f,'r')
        for line in File:
            try:
                solventInfo[line.split(',')[0]] = line.split(',')[1].strip()
            except:
                pass
        File.close()
        addInfo[f.split('_')[0]] = solventInfo

foo = pd.DataFrame()     
for key in addInfo.keys():
    solutes = addInfo[key].keys()
    dgList = addInfo[key].values()
    
    datTemp = pd.DataFrame(data=zip(*[solutes,dgList]), columns=[row1,key])
    foo = foo.combine_first(datTemp)

final = driverFile.combine_first(foo)
cols = final.columns.tolist()
colInd = len(cols)
cols = cols[-(colInd-2):] + cols[:-(colInd-2)]
final = final[cols]
final.to_csv(path+nameThisFile)
