dataFile = open('Water-Data.txt','r')
soluteNames = open('mnsol/water.csv','r')

errors = []
for line in dataFile:
    if line.split(',')[1].strip() != 'errfinal':
        errors.append(line.split(',')[1].strip())

solutes = []   
for line in soluteNames:
    solutes.append(line.split(',')[0])

dataFile.close()
soluteNames.close()

