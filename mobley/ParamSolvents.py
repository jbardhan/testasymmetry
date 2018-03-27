import matlab.engine
import os
 
Home = os.environ['HOME']
os.chdir(Home+'/repos/testasymmetry/mobley')
 
Run_Type = 'param'
 
if Run_Type == 'param':
    Solvent_List = ['paramWater']
elif Run_Type == 'run':
    Solvent_List = ['runWater']
 
Solvent_List = Solvent_List
 
if len(Solvent_List) <= 1:
    for solvent in Solvent_List:
        eng = matlab.engine.start_matlab()
        eng.eval(solvent,nargout=0)
        eng.quit()
else:
    eng = matlab.engine.start_matlab()
    eng.eval(Solvent_List,nargout=0)
    eng.quit()