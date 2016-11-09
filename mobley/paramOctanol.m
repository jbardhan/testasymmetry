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
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
saveMemory = 0;
writeLogfile = 0;
logfileName = 'junklogfile';

epsIn  =  1;
Tbase = 300; 
mytemp=Tbase;
KelvinOffset = 273.15;
epsOut = 10.3; % Zhao+Abraham J. Org. Chem 2005
conv_factor = 332.112;
staticpotential = 2.0; % this only affects charged molecules;
kappa = 0.0;  % should be zero, meaning non-ionic solutions!


% the staticpotential below should not be used any more, please check
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);
     
[mol_list,dG_list,surfArea_list]=textread('mnsol/octanol.csv',...
					  '%s %f %f','delimiter',',');

%testset = {'2_methylpropane','ethane','methane','n_butane','n_heptane','n_hexane',...
%	   'n_octane','n_pentane','propane'};

testset  = {'toluene','ethanol','butanone','nitromethane','n_octane','14_dioxane','octan_1_ol'};
	    


% all octanol available side chain analogues 
%testset = {'2_methylpropane', 'acetic_acid', 'ethanol', 'methane', 'methanol',...
% 'n_butane', 'n_butylamine', 'p_cresol', 'propane', 'propanoic_acid','toluene'};

% complete list of side chain analogues. not available for all solvents
%testset = {'1_methyl_imidazole','2_methylpropane', ...
%	   '3_methyl_1h_indole','acetic_acid','ethanamide', ...
%	   'ethanol','methane','methanethiol','methanol', ...
%	   'methyl_ethyl_sulfide','n_butane','n_butylamine', ...
%	   'p_cresol','propane','propanoic_acid','toluene'};
curdir=pwd;
for i=1:length(testset)
  dir=sprintf('%s/Dropbox/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',getenv('HOME'),testset{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;%chargeDistribution;
  foo = strcmp(mol_list,testset{i});
  index = find(foo);
  if length(index) ~= 1
    fprintf('error finding refdata!\n');
    keyboard
  end
  referenceData{i} = dG_list(index);
  surfArea{i} = surfArea_list(index);
  chdir(curdir);
  addProblemSA(testset{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end

pqrData = struct('xyz', [0 0 0], 'q', 1, 'R', 1);

% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.

NaReference = -97.3; NaR = 0.92*1.41075; NaSurfArea = 4*pi*NaR^2;
KReference  = -79.9; KR = 0.92*1.76375; KSurfArea = 4*pi*KR^2;
ClReference = -66.6; ClR = 0.92*2.27; ClSurfArea = 4*pi*ClR^2;

addProblemSA('Na',pqrData,'../born/Na_2.srf',1, ...
	   NaReference,NaSurfArea);
addProblemSA('K',pqrData,'../born/K_2.srf',1, ...
	   KReference,KSurfArea);
addProblemSA('Cl',pqrData,'../born/Cl_2.srf',-1, ...
	   ClReference,ClSurfArea);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0.5  -60  -0.5 0.0 -0.03 1.6];
lb = [-2 -200 -100 -20 -0.1 0];
ub = [+2 +200 +100 +20.1 +0.1 +4];

options = optimoptions('lsqnonlin','MaxIter',8);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEMSA(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
[err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);