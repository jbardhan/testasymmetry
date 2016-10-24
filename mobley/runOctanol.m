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
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet
epsIn  =  1;
Tbase = 300; mytemp=Tbase;
KelvinOffset = 273.15;
epsOut = 10.3; % Zhao+Abraham J. Org. Chem 2005
conv_factor = 332.112;
staticpotential = 2.0; % this only affects charged molecules;

kappa = 0.0;  % should be zero, meaning non-ionic solutions!
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);

% here we define the actual params for the NLBC test
asymParams = struct('alpha',0.5, 'beta', -100,'EfieldOffset',1); 
     
[mol_list,dG_list,surfArea_list]=textread('mnsol/octanol.csv',...
					  '%s %f %f','delimiter',',');

curdir = pwd;
for i=1:length(mol_list)
  dir=sprintf('/Volumes/Bardhan2TB/nlbc-mobley/nlbc_test/%s',mol_list{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;
  referenceData{i} = dG_list(i);
  surfArea{i} = surfArea_list(i);
  chdir(curdir);
  addProblemSA(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end

x = [1.98 -86.3423 -2.5825 18.8733 -0.0162 1.5785];
[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA(x);