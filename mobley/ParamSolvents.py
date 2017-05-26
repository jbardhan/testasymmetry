import matlab.engine
import os

Home = os.environ['HOME']
os.chdir(Home+'/repos/testasymmetry/mobley')

Run_Type = 'run'

if Run_Type == 'param':
	Solvent_List = ['paramPentane', 'paramHexane', 'paramHeptane', 'paramOctane', 'paramNonane', 'paramDecane']
elif Run_Type == 'run':
	Solvent_List = ['runPentane', 'runHexane', 'runHeptane', 'runOctane', 'runNonane', 'runDecane']

Solvent_List = Solvent_List[:]

if len(Solvent_List) <= 6:
	for solvent in Solvent_List:
		eng = matlab.engine.start_matlab()
		eng.eval(solvent,nargout=0)
		eng.quit()
else:
	eng = matlab.engine.start_matlab()
	eng.eval(Solvent_List,nargout=0)
	eng.quit()

