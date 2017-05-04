import os
from os import listdir
from os.path import join, isfile
import pandas as pd
home = os.path.expanduser('~')
path = home+'/repos/testasymmetry/mobley/mnsol/'
os.chdir(path)

files = [f for f in listdir(path) if isfile(join(path,f))]

driverFile = pd.read_csv(os.path.join(path,'New_data.csv'))

dontWant = ['New_data', 'mobley_sa', '.DS_Store']

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
    elif (f.split('.')[0] not in solvents and f.split('.')[0] not in dontWant):
        # No Ions not in New_Data
        File = open(f,'r')
        for line in File:
            try:
                solventInfo[line.split(',')[0]] = line.split(',')[1].strip()
            except:
                pass
        File.close()
        addInfo[f.split('.')[0]] = solventInfo
    elif f.split('_')[0] not in solvents:
        # Ions not in New_Data
        File = open(f,'r')
        for line in File:
            try:
                solventInfo[line.split(',')[0]] = line.split(',')[1].strip()
            except:
                pass
        File.close()
        addInfo[f.split('_')[0]] = solventInfo
        
        
        
        

