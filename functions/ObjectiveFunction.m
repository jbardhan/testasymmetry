function Error = ObjectiveFunction(i,Params,E)
%ERROR          Returns the error via lsm of the MD FEP energy minus the
%               calculated FEP free energy.


MD = [1 2 3 4 5 6 7 8 ];

for i = 1 : length(MD)
    Error(i) = [(MD(i) - E(i)).^2];
end