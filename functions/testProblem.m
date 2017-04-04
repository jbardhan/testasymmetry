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
pdbFile = '../saltresidues/asp/asp.pdb';
crgFile = '../saltresidues/asp/jd2.crg';
pqrData = loadPdbAndCrg(pdbFile,crgFile);
% the PQR data structure is widely used and we'll make a lot of use
% of it.  

% this file holds all the references to the mesh data (triangles)
% -- have a look at the file format at some point, and see if you
% can read the initializeProblem.m function
srfSternFile   = '../saltresidues/asp/asp_scaledcharmm_stern_1.srf';

% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.
loadTestReferenceAndChargeDistribution

% this is the first of the main two steps.  it adds the given
% problem geometry to the list of tests.  
addProblem('asp_deprot_1',pqrData,srfSternFile,chargeDistribution, ...
	   referenceE);

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
[calculatedE,referenceE]=calculateProblem(ProblemSet(1), ...
					  asymParams);


