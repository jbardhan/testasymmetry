% Path information
addpath('..');
addpath('../');
addpath('../functions/');
addpath('../born/');
addpath('../../pointbem');
addpath('../../panelbem');

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
global UsefulConstants ProblemSet
epsIn  =  1;
epsOut = 2;
conv_factor = 332.112;
staticpotential = 0;
kappa = 0.0;  % for now, this should be zero, meaning non-ionic solutions!
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);

% here we define the actual params for the NLBC test
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

% '1_methyl_imidazole'   -> his
% '2_methylpropane' ->  LEU (isobutane)
% '3_methyl_1h_indole'  -> trp
% 'acetic_acid' -> asp
% 'ethanamide' -> ASN (acetamide)
% 'ethanol' -> thr
% 'methane' ->  ALA
% 'methanethiol' -> cys
% 'methanol' -> ser
% 'methyl_ethyl_sulfide' -> met
% 'n_butane' -> ILE (butane)
% 'n_butylamine' ->lys
% 'p_cresol' -> tyr
% 'propane'  ->  VAL
% 'propanoic_acid' ->  GLU (propionic acid)
% 'toluene' ->   phe

analogs = {'1_methyl_imidazole','2_methylpropane', ...
	   '3_methyl_1h_indole','acetic_acid','ethanamide', ...
	   'ethanol','methane','methanethiol','methanol', ...
	   'methyl_ethyl_sulfide','n_butane','n_butylamine', ...
	   'p_cresol','propane','propanoic_acid','toluene'};

for i=1:length(analogs)
  chdir(analogs{i});
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_1.srf',analogs{i});
  loadTestReferenceAndChargeDistribution
  chargeDist{i} = chargeDistribution;
  loadVilla02cyclohexane
  referenceData{i} = referenceE;
  chdir('..');

  addProblem(analogs{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0.5 -60 -0.5];
lb = [0 -Inf -Inf];
ub = [Inf 0 0];

options = optimoptions('lsqnonlin');
options = optimoptions(options,'Display', 'off');
y = @(x)ObjectiveFromBEM(x);
x = lsqnonlin(y,x0,lb,ub,options);
