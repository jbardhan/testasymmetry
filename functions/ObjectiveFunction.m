function Error = ObjectiveFunction(i,Params,E)
%ERROR          Returns the error via lsm of the MD FEP energy minus the
%               calculated FEP free energy.  OBJECTIVEFUNCTION(i,Params,E)
%               takes in i situations and calculates the square of the
%               difference between the MD, experimental energy of solvation, 
%               MD(i) and the calculated energy of solvation E(i).
%               Eventually, we will want to minimize this function by
%               varying 'Params' to get a best fit to the data.
%               


MD = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];

for i = 1 : length(MD)
    Error(i) = [(MD(i) - E(i)).^2];
end