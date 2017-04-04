function Error = ObjectiveFunction(x)
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


MD = [ -93.4    % Na+
       -175.7   % Na-
       -57.0    % Cl+
       -95.3    % Cl-
       -73.4    % K+
       -128.3895 % K-
       -66.78   % Rb+
       -114.1   % Rb-
       -108.6   % Mg+
       -218.5   % Mg-
       -60.42   % Cs+
       -101.9   % Cs-
       -88.91   % Ca+
       -163.4   % Ca-
       -67.03   % Ba+
       -115.1   % Ba-
       -99.05   % Zn+
       -191.2   % Zn-
       -89.08   % Cd+
       -164.3   % Cd-
       ]';
   
Params = struct('alpha',x(1), 'beta', x(2),'EfieldOffset',x(3));
   
F = Level_2(Params); % Calls the function Level_2.m which calculates the energy using bornPicardNoStern.

Error = MD - F;
