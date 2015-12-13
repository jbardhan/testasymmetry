function Error_MSA = ObjectiveFunction_MSA(x)
%ERROR          Returns the deviance of the experimental MD FEP energy results from the
%               calculated FEP free energy.  OBJECTIVEFUNCTION(Params)
%               takes in i situations and calculates the difference between the MD,
%               experimental energy of solvation, 
%               MD(i) and the calculated energy of solvation E(i).
%               
%               

addpath('..')
addpath('../born/')
addpath('../../pointbem/')

delS = [-199 -143 -101 -92 -78 -138 -89 -79 -66]';
%t = [25 30]; 
%for i = 1:length(t)-1 (t(i+1) - t(i)) 

MD = [ -483 - 5*delS(1)    % Li+
       -403 - 5*delS(2)    % Na+
       -333 - 5*delS(3)    % K+
       -316 - 5*delS(4)    % Rb+
       -288 - 5*delS(5)    % Cs+
       -396 - 5*delS(6)    % F-
       -310 - 5*delS(7)    % Cl-
       -291 - 5*delS(8)    % Br-
       -264 - 5*delS(9)    % I-
       ]';
%end  

Params = struct('alpha',x(1), 'beta', x(2),'EfieldOffset',x(3));
   
F = Level_2_MSA(Params); % Calls the function Level_2.m which calculates the energy using bornPicardNoStern.

Error_MSA = MD - F;