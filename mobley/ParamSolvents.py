import matlab.engine
import os

Home = os.environ['HOME']
os.chdir(Home+'/repos/testasymmetry/mobley')

Run_Type = 'run'

if Run_Type == 'param':
	Solvent_List = ['paramWater', 'paramOctanol', 'paramHexadecane', 'paramChloroform',
 		'paramCyclohexane','paramCarbontet', 'paramHexane', 'paramToluene', 'paramXylene']
elif Run_Type == 'run':
	Solvent_List = ['runWater', 'runOctanol', 'runHexadecane', 'runChloroform',
 		'runCyclohexane','runCarbontet', 'runHexane', 'runToluene', 'runXylene']

Solvent_List = Solvent_List

if Run_Type == 'param':
	numRuns = 10
	eng = matlab.engine.start_matlab()
	for i in range(numRuns):
		eng.eval('genRandomTestSet',nargout=0)
		if len(Solvent_List) <= 9:
			for solvent in Solvent_List:
				eng.eval(solvent,nargout=0)
		else:
			eng.eval(Solvent_List,nargout=0)

elif Run_Type == 'run':
	eng = matlab.engine.start_matlab()
	eng.eval('getOptValues',nargout=0)
	if len(Solvent_List) <= 9:
			for solvent in Solvent_List:
				eng.eval(solvent,nargout=0)
	else:
		eng.eval(Solvent_List,nargout=0)


