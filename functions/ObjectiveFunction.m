function Error = ObjectiveFunction(i,Params,E)
%ERROR          Returns the error via lsm of the MD FEP energy minus the
%               calculated FEP free energy.  OBJECTIVEFUNCTION(i,Params,E)
%               takes in i situations and calculates the square of the
%               difference between the MD, experimental energy of solvation, 
%               MD(i) and the calculated energy of solvation E(i).
%               Eventually, we will want to minimize this function by
%               varying 'Params' to get a best fit to the data.
%               


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
       ];

for i = 1 : length(MD)
    Error(i) = [(MD(i) - E(i)).^2];
end