import os
import pandas as pd
from collections import Counter
import linecache
Home = os.path.expanduser('~')

solventList = ['Hexadecane','Cyclohexane','Chloroform','Carbontet']
soluteDict = {}
for solvent in solventList:
    soluteNames = open(('mnsol/'+solvent.lower()+'.csv'),'r')

    solutes = []
    for line in soluteNames:
    	solutes.append(line.split(',')[0])
    soluteNames.close() 
    soluteDict[solvent] = solutes

unionSolutes = []
for vallists in soluteDict.values():
	for val in vallists:
		if val not in unionSolutes:
			unionSolutes.append(val)

check = Counter(unionSolutes)

for val in check.values():
	if val > 1:
		print('Duplicate Entry')

soluteData = {}
for solute in unionSolutes:
	errorForSolutePerSolvent = []
	for solvent in solventList:
		there = 0
		soluteNames = open(('mnsol/'+solvent.lower()+'.csv'),'r')
		for i,line in enumerate(soluteNames):
			if solute == line.split(',')[0]:
				errorLine = linecache.getline(solvent+'-Data.txt',i+2)
				errorForSolutePerSolvent.append(errorLine.split(',')[1])
				there = 1
				break
		if not there:
			errorForSolutePerSolvent.append(' ')
	soluteData[solute] = errorForSolutePerSolvent

    	soluteNames.close()

dFrame = pd.DataFrame.from_dict(soluteData,orient='index')
dFrame.columns = solventList
dFrame.to_csv(os.path.join(Home+'/repos/testasymmetry/mobley','ErrorOfSoluteBySolvent.csv'),header=True)