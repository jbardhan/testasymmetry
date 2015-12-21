function Error_MSA = ObjectiveFunction_MSA(x)
%ERROR_MSA      Returns the deviance of the experimental MD FEP energy
%               results from Fawcett Chp. 3 and
%               calculated FEP free energy.  OBJECTIVEFUNCTION(Params)
%               takes in i situations and calculates the difference between the MD,
%               experimental energy of solvation, 
%               MD(i) and the calculated energy of solvation E(i).
%               
%               

addpath('..')
addpath('../born/')
addpath('../../pointbem/')

delS = [-199 -143 -101 -92 -78 -138 -89 -79 -66].*(0.239/10^3)'; % Delta S at 25 C from Fawcett Ch. 3
delG = [-483 -403 -333 -316 -288 -396 -310 -291 -264].*0.239'; % Delta G at 25 C from Fawcett Ch. 3
t = 0; % Desired Temp in C
del_t = (t-25); 

MD = [ delG(1)-del_t*delS(1)    % Li+
       delG(2)-del_t*delS(2)    % Na+
       delG(3)-del_t*delS(3)    % K+
       delG(4)-del_t*delS(4)    % Rb+
       delG(5)-del_t*delS(5)    % Cs+
       delG(6)-del_t*delS(6)    % F-
       delG(7)-del_t*delS(7)    % Cl-
       delG(8)-del_t*delS(8)    % Br-
       delG(9)-del_t*delS(9)    % I-
       ]'; % The constant is Delta G at 25 C from Fawcett Ch. 3
 

Params = struct('alpha',x(1), 'beta', x(2),'EfieldOffset',x(3));
   
F = Level_2_MSA(Params); % Calls the function Level_2_MSA.m which calculates the energy using bornPicardNoStern.

Error_MSA = MD - F;