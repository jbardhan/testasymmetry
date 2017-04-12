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

if len(Solvent_List) <= 9:
	for solvent in Solvent_List:
		eng = matlab.engine.start_matlab()
		eng.eval(solvent,nargout=0)
		eng.quit()
else:
	eng = matlab.engine.start_matlab()
	eng.eval(Solvent_List,nargout=0)
	eng.quit()

