import matlab.engine
import os

Home = os.environ['HOME']
os.chdir(Home+'/repos/testasymmetry/mobley')

Run_Type = 'run'

if Run_Type == 'param':
	Solvent_List = ['paramWater','paramCarbontet','paramOctanol','paramToluene','paramXylene','paramHexane','paramHexadecane','paramCyclohexane','paramBenzene','paramChlorobenzene','paramDecane','paramBromobenzene','paramDiethylether','paramDibutylether','paramIsooctane','paramIodobenzene']
elif Run_Type == 'run':
	Solvent_List = ['runWater','runOctanol','runCarbontet','runToluene','runBenzene','runXylene','runHexane','runHexadecane','runCyclohexane','runChlorobenzene','runBromobenzene','runDiethylether','runDibutylether','runIsooctane','runDecane','runIodobenzene']

Solvent_List = Solvent_List[14:]

if len(Solvent_List) <= 16:
	for solvent in Solvent_List:
		eng = matlab.engine.start_matlab()
		eng.eval(solvent,nargout=0)
		eng.quit()
else:
	eng = matlab.engine.start_matlab()
	eng.eval(Solvent_List,nargout=0)
	eng.quit()

