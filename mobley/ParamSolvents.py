import matlab.engine
import os

Home = os.environ['HOME']
os.chdir(Home+'/repos/testasymmetry/mobley')

Run_Type = 'param'

if Run_Type == 'param':
	Solvent_List = ['paramWater','paramCarbontet','paramOctanol','paramToluene','paramBenzene','paramXylene','paramHexane','paramHexadecane','paramCyclohexane','paramChlorobenzene','paramDecane','paramBromobenzene','paramDiethylether','paramDibutylether','paramIsooctane','paramIodobenzene']
elif Run_Type == 'run':
	Solvent_List = ['runWater','runCarbontet','runOctanol','runToluene','runBenzene','runXylene','runHexane','runHexadecane','runCyclohexane','runChlorobenzene','runDecane','runBromobenzene','runDiethylether','runDibutylether','runIsooctane','runIodobenzene']

Solvent_List = Solvent_List[:]

if len(Solvent_List) <= 16:
	for solvent in Solvent_List:
		eng = matlab.engine.start_matlab()
		eng.eval(solvent,nargout=0)
		eng.quit()
else:
	eng = matlab.engine.start_matlab()
	eng.eval(Solvent_List,nargout=0)
	eng.quit()

