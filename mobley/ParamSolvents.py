import matlab.engine
import os

Home = os.environ['HOME']
os.chdir(Home+'/repos/testasymmetry/mobley')

Run_Type = 'param'

if Run_Type == 'param':
	Solvent_List = ['paramWater', 'paramOctanol', 'paramMethanol', 'paramAcetonitrile', 'paramNitromethane', 'paramPropanol',
 		        'paramDichloroethane','paramDimethylsulfoxide', 'paramDimethylformamide', 'paramButanol', 'paramEthanol', 'paramNitrobenzene']
elif Run_Type == 'run':
	Solvent_List = ['runWater', 'runOctanol', 'runMethanol', 'runAcetonitrile', 'runNitromethane', 'runPropanol',
 		        'runDichloroethane','runDimethylsulfoxide', 'runDimethylformamide', 'runButanol', 'runEthanol', 'runNitrobenzene']

Solvent_List = Solvent_List[:]

if len(Solvent_List) <= 12:
	for solvent in Solvent_List:
		eng = matlab.engine.start_matlab()
		eng.eval(solvent,nargout=0)
		eng.quit()
else:
	eng = matlab.engine.start_matlab()
	eng.eval(Solvent_List,nargout=0)
	eng.quit()

