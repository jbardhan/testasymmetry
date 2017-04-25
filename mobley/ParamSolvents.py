import matlab.engine
import os

Home = os.environ['HOME']
os.chdir(Home+'/repos/testasymmetry/mobley')

Run_Type = 'param'

if Run_Type == 'param':
	Solvent_List = ['paramWater', 'paramOctanol', 'paramHexadecane', 'paramChloroform',
 		'paramCyclohexane','paramCarbontet', 'paramHexane', 'paramToluene', 'paramXylene']
elif Run_Type == 'run':
	Solvent_List = ['runWater', 'runOctanol', 'runHexadecane', 'runChloroform',
 		'runCyclohexane','runCarbontet', 'runHexane', 'runToluene', 'runXylene']

Solvent_List = Solvent_List

if len(Solvent_List) <= 9:
	eng = matlab.engine.start_matlab()
	eng.eval('genRandomTestSet',nargout=0)
	for solvent in Solvent_List:
		eng.eval(solvent,nargout=0)
else:
	eng = matlab.engine.start_matlab()
	eng.eval('genRandomTestSet',nargout=0)
	eng.eval(Solvent_List,nargout=0)

