NaReference = [-93.4 -175.7]';
ClReference = [-57.0 -95.3]';
KReference  = [-73.4 -128.3895]';
RbReference = [-66.78 -114.1]';
MgReference = [-108.6 -218.5]';
CsReference = [-60.42 -101.9]';
CaReference = [-88.91 -163.4]';
BaReference = [-67.03 -115.1]';
ZnReference = [-99.05 -191.2]';
CdReference = [-89.08 -164.3]';

% for the ProblemSet structure, create a charge vector for each
% test.  Here, that means Natoms simulations.  Each one (the ith)
% has the q(i) charge and no others.
chargeDistribution = [1 -1];