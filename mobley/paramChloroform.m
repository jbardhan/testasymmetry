% Path information
Home = getenv('HOME');
addpath(sprintf('%s/repos/pointbem',Home));
addpath(sprintf('%s/repos/panelbem',Home));
addpath(sprintf('%s/repos/testasymmetry',Home));
addpath(sprintf('%s/repos/testasymmetry/functions',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley',Home));
addpath(sprintf('%s/repos/testasymmetry/born',Home));

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet
epsIn  =  1;
Tbase = 293; mytemp=Tbase;
KelvinOffset = 273.15;
epsOut = 4.81;  % http://macro.lsu.edu/howto/solvents/Dielectric%20Constant%20.htm
conv_factor = 332.112;
staticpotential = 2.0; % this only affects charged molecules;

kappa = 0.0;  % should be zero, meaning non-ionic solutions!
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);

% here we define the actual params for the NLBC test
asymParams = struct('alpha',0.5, 'beta', -100,'EfieldOffset',1); 
     

analogs = {'3_methyl_1h_indole','toluene','methanol','ethanamide', ...
	   'n_butylamine','acetic_acid','propanoic_acid'};

%analogs = {'1_methyl_imidazole','2_methylpropane', ...
%	   '3_methyl_1h_indole','acetic_acid','ethanamide', ...
%	   'ethanol','methane','methanethiol','methanol', ...
%	   'methyl_ethyl_sulfide','n_butane','n_butylamine', ...
%	   'p_cresol','propane','propanoic_acid','toluene'};

for i=1:length(analogs)
  chdir(analogs{i});
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',analogs{i});
  loadVilla02chloroform
  chargeDist{i} = chargeDistribution;
  referenceData{i} = referenceE;
  chdir('..');
  addProblem(analogs{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i});
end

pqrData = struct('xyz', [0 0 0], 'q', 1, 'R', 1);

% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.
LoadIonChloroformReferenceAndChargeDataAtTemp

addProblem('Na',pqrData,'../born/Na_2.srf',CationChargePlusOne, ...
	   NaReference);
addProblem('K',pqrData,'../born/K_2.srf',CationChargePlusOne, ...
	   KReference);
addProblem('Cs',pqrData,'../born/Cs_2.srf',CationChargePlusOne, ...
	   CsReference);
%addProblem('F',pqrData,'../born/K_2.srf',AnionChargeMinusOne, ...
%	   FReference);
addProblem('Cl',pqrData,'../born/Cl_2.srf',AnionChargeMinusOne, ...
	   ClReference);

length(ProblemSet)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0.9786  -40.2471  -0.5590];
lb = [-4 -200 -100];
ub = [+4 +200 +100];

options = optimoptions('lsqnonlin','MaxIter',6);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEM(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
