% Path information
addpath('/Users/jbardhan/repos/pointbem');
addpath('/Users/jbardhan/repos/panelbem');
addpath('/Users/jbardhan/repos/testasymmetry');
addpath('/Users/jbardhan/repos/testasymmetry/functions');
addpath('/Users/jbardhan/repos/testasymmetry/mobley');
addpath('/Users/jbardhan/repos/testasymmetry/born/');

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
global UsefulConstants ProblemSet
epsIn  =  1;
Tbase = 300; mytemp=Tbase;
KelvinOffset = 273.15;
epsOut = epsilon_t(Tbase);
fprintf('Using temperature %f K.\n',Tbase);
conv_factor = 332.112;
staticpotential = 10.7;
kappa = 0.0;  % for now, this should be zero, meaning non-ionic solutions!
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);

% here we define the actual params for the NLBC test
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

%asymParams = struct('alpha',0.326658, 'beta', -47.770653, ...
%		    'EfieldOffset',-1.198492);

% on 3/7/16 with 12 iters, no discharge, at density=2
%asymParams = struct('alpha',0.3233,'beta',-46.9515,...
% 'EfieldOffset',-1.1802); 

chdir('saltresidues93/phe/');
pqrData = loadPdbAndCrg('phe.pdb','phe.crg');
pqrAll{1} = pqrData;
srfFile{1} = 'saltresidues93/phe/phe_scaledcharmm_stern_2.srf';
loadAtomReferenceAndChargeDistribution
chargeDist{1} = chargeDistribution;
referenceData{1} = referenceE;
addProblem('phe',pqrAll{1},srfFile{1},chargeDist{1}, ...
	   referenceData{1});
chdir('../..');

chdir('saltresidues93/arg/');
pqrData = loadPdbAndCrg('arg.pdb','jr1.crg');
pqrAll{2} = pqrData;
srfFile{2} = 'saltresidues93/arg/arg_scaledcharmm_stern_2.srf';
loadAtomReferenceAndChargeDistribution
chargeDist{2} = chargeDistribution;
referenceData{2} = referenceE;
addProblem('arg',pqrAll{2},srfFile{2},chargeDist{2}, ...
	   referenceData{2});
chdir('../..');

chdir('saltresidues93/asp/');
pqrData = loadPdbAndCrg('asp.pdb','jd2.crg');
pqrAll{3} = pqrData;
srfFile{3} = 'saltresidues93/asp/asp_scaledcharmm_stern_2.srf';
loadAtomReferenceAndChargeDistribution
chargeDist{3} = chargeDistribution;
referenceData{3} = referenceE;
addProblem('asp',pqrAll{3},srfFile{3},chargeDist{3}, ...
	   referenceData{3});
chdir('../..');

pqrData = struct('xyz', [0 0 0], 'q', 1, 'R', 1);

% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.
LoadExperimentReferenceAndChargeDataAtTemp

addProblem('Na',pqrData,'born/Na_2.srf',CationChargePlusOne, ...
	   NaReference);
addProblem('K',pqrData,'born/K_2.srf',CationChargePlusOne, ...
	   KReference);
addProblem('Rb',pqrData,'born/Rb_2.srf',CationChargePlusOne, ...
	   RbReference);
addProblem('Cs',pqrData,'born/Cs_2.srf',CationChargePlusOne, ...
	   CsReference);
addProblem('Cl',pqrData,'born/Cl_2.srf',AnionChargeMinusOne, ...
	   ClReference);

length(ProblemSet)

[calculatedE_initial,referenceE_initial] = CalculateEnergiesFromBEM(asymParams); 

%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0.5 -60 -0.5];
lb = [0 -Inf -Inf];
ub = [Inf 0 0];

options = optimoptions('lsqnonlin','MaxIter',4);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEM(x);
[x,resnorm,residual,exitflag,output] = lsqnonlin(y,x0,lb,ub,options);

NEWasymParams = MakeParamsStruct(x);
[calculatedE_final,referenceE_final] = CalculateEnergiesFromBEM(asymParams); 
