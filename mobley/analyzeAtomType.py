import os
import pandas as pd
Home = os.path.expanduser('~')
prmPath = Home + '/Dropbox/lab/Projects/slic-jctc-mnsol/nlbc-mobley/prmcrd'
csvPath = Home + '/repos/testasymmetry/mobley/'

solventList = ['Water','Octanol','Chloroform','Hexane','Xylene','Hexadecane','Cyclohexane','Toluene','Carbontet']
for solvent in solventList:
    dataFile = open((solvent+'-Data.txt'),'r')
    soluteNames = open(('mnsol/'+solvent.lower()+'.csv'),'r')

    errors = [] # Get the errors from the Solvent-Data.txt file
    for line in dataFile:
        if line.split(',')[1].strip() != 'errfinal':
            errors.append(line.split(',')[1].strip())

    solutes = [] # Get the solute names from the solvent.csv file
    for line in soluteNames:
        solutes.append(line.split(',')[0])
    dataFile.close()
    soluteNames.close()
    
    typeLineIndex = [] # Get the line index for atom types from each solute.prtmtop file for all solutes from solvent.csv
    for solute in solutes:
        prmFile = open(os.path.join(prmPath,(solute+'.prmtop')))
        for i, line in enumerate(prmFile):
            if 'FLAG AMBER_ATOM_TYPE' in line:
                typeLineIndex.append(i + 2)            
        prmFile.close()

    index = 0  
    atomTypes = [] # gets the atom types in each prmtop file
    for solute in solutes:
        prmFile = open(os.path.join(prmPath,(solute+'.prmtop')))
        lines = prmFile.readlines()
        temp = lines[typeLineIndex[index]].strip()
        atomTypes.append(temp.split())
        prmFile.close()
        index += 1
        
    errDict = {} # Creates a dictionary with keys = abs(error) and vals = [solute_name, atom_types]
    index = 0
    for error in errors:
        data = [solutes[index],atomTypes[index]]
        errDict[abs(float(error))] = data
        index += 1
        
    colNames = ['Solute Name', 'Atom Type']
    dFrame = pd.DataFrame.from_dict(errDict,orient='index')
    dFrame.columns = colNames
    dFrame.to_csv(os.path.join(csvPath,('errorByAtomTypeFor'+solvent+'.csv')),header=True) # Stores the dictionary as a csv file    