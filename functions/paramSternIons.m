% Path information
addpath('..');
addpath('../functions/');
addpath('../born/');
addpath('../../pointbem');
addpath('../../panelbem');

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
global UsefulConstants ProblemSet
epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
staticpotential = 10.7;
kappa = 0.0;  % for now, this should be zero, meaning non-ionic solutions!
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);

% here we define the actual params for the NLBC test
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

% the PDB defines where the atoms are located
% the CRG defines what the atom charges are in the "wild type" 
pqrData = struct('xyz', [0 0 0], 'q', 1, 'R', 1);

% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.
loadBornReferenceAndChargeDistribution

addProblem('Na',pqrData,'../born/Na_2.srf',chargeDistribution, ...
	   NaReference);
addProblem('K',pqrData,'../born/K_2.srf',chargeDistribution, ...
	   KReference);
addProblem('Rb',pqrData,'../born/Rb_2.srf',chargeDistribution, ...
	   RbReference);
addProblem('Cs',pqrData,'../born/Cs_2.srf',chargeDistribution, ...
	   CsReference);
addProblem('Mg',pqrData,'../born/Mg_2.srf',chargeDistribution, ...
	   MgReference);
addProblem('Ca',pqrData,'../born/Ca_2.srf',chargeDistribution, ...
	   CaReference);
addProblem('Ba',pqrData,'../born/Ba_2.srf',chargeDistribution, ...
	   BaReference);
addProblem('Zn',pqrData,'../born/Zn_2.srf',chargeDistribution, ...
	   ZnReference);
addProblem('Cd',pqrData,'../born/Cd_2.srf',chargeDistribution, ...
	   CdReference);
addProblem('Cl',pqrData,'../born/Cl_2.srf',chargeDistribution, ...
	   ClReference);

% now ProblemSet is an "array" of length one. 

% calculateProblem runs all the tests associated with a given entry
% in the ProblemSet array.  for example:
%  1) for a Born ion we might have one test per example (say, +1
%  for sodium)
%  2) for a Born ion we might also have TWO tests (+1 and -1)
%  3) for an amino acid we might have ONE test per example
%  (charging all)
%  4) for an amino acid we might have Natoms tests per example
%  (charging each atom individually)
%  5) for an amino acid we might have 2*Natoms test per example
%  (charging each individually, and charging all BUT that one)

x0 = [0.5 -60 -0.5];
lb = [0 -Inf -Inf];
ub = [Inf 0 0];

options = optimoptions('lsqnonlin');
options = optimoptions(options,'Display', 'off');
y = @(x)ObjectiveFromBEM(x);
x = lsqnonlin(y,x0,lb,ub,options);


