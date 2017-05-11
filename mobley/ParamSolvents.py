import matlab.engine
import os

Home = os.environ['HOME']
os.chdir(Home+'/repos/testasymmetry/mobley')

Run_Type = 'run'

if Run_Type == 'param':
	Solvent_List = ['paramButanol', 'paramWater', 'paramOctanol', 'paramMethanol', 'paramPropanol', 'paramEthanol']
elif Run_Type == 'run':
	Solvent_List = ['runButanol', 'runWater', 'runOctanol', 'runMethanol', 'runPropanol', 'runEthanol']

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

