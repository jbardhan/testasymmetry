% Path information
clear all
Home = getenv('HOME');
repo_path = sprintf('%s/repos',Home);
dropbox_path = sprintf('%s/Dropbox',Home);
addpath(sprintf('%s/repos/pointbem',Home));
addpath(sprintf('%s/repos/panelbem',Home));
addpath(sprintf('%s/repos/testasymmetry',Home));
addpath(sprintf('%s/repos/testasymmetry/functions',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley/slic-COSMO',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley/slic-COSMO/tests',Home));
addpath(sprintf('%s/repos/testasymmetry/born',Home));

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
logfileName = 'water.out';
epsOut = 78.36;
allData = readtable('../all_Bondii_data.csv'); 

ParamWatInfo = load('OptCosmoBondiiFinal_4_HB.mat');
x = ParamWatInfo.x;
training_set = ParamWatInfo.training_set;
mol_list = allData.solute;
dG_list = allData.dG_expt;

saveMemory = 1;
writeLogfile = 1;
epsIn  =  1;
Tbase = 300; mytemp=Tbase;

KelvinOffset = 273.15;
conv_factor = 332.112;
staticpotential = 0.0; % this only affects charged molecules;

kappa = 0.0;  % should be zero, meaning non-ionic solutions!
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);
%COSMO H-bond atom types
allHbondTypes = {'n_amine','n_amide','n_nitro',...
                 'n_other','o_carbonyl','o_ester',...
                 'o_nitro','o_hydroxyl','fluorine',...
                 'h_oh','h_nh','h_other'};

% here we define the actual params for the NLBC test
solventAreas = allData{495,9:79};
solventATypes = allData{495,80:144};
solventHbData = allData{495,145:end};
solventVdWA = allData{495,11};
solventVdWV = allData{495,12};
temperature = 24.85 + KelvinOffset;
curdir = pwd;
for i=1:length(mol_list)
  dir=sprintf('%s/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',dropbox_path,mol_list{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;%chargeDistribution;
  soluteAtomAreas{i} = allData{i,9:79};
  soluteAtomTypes{i} = {allData{i,80:144}};
  soluteHbondData{i} = allData{i,145:end};
  solute_VdWA{i} = allData{i,11};
  solute_VdWV{i} = allData{i,12};
  solventAtomAreas{i} = solventAreas;
  solventAtomTypes{i} = {solventATypes};
  solventHbondData{i} = solventHbData;
  solvent_VdWA{i} = solventVdWA;
  solvent_VdWV{i} = solventVdWV;
  atom_vols{i} = allData{i,14};
  temp{i} = temperature;
  referenceData{i} = dG_list(i);
  newHB{i}=1;
  chdir(curdir);
  addProblemCosmo(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},...
                  soluteAtomAreas{i},soluteAtomTypes{i},soluteHbondData{i},...
                  solute_VdWV{i},solute_VdWA{i},...
                  solventAtomAreas{i},solventAtomTypes{i},solventHbondData{i},...
                  solvent_VdWV{i},solvent_VdWA{i},...
                  atom_vols{i},temp{i},newHB{i});



end
disp_mob = allData.disp_mobley; 
cav_mob = allData.cav_mobley; 
np_mob = allData.np_mobley; 
es_mob = allData.es_mobley; 
np_SLIC = allData.np_SLIC; 
es_SLIC= allData.es_SLIC;
% xNoIonsInParam = [  1.8570  -75.6855   -2.9884    8.7302    0.0113    0.4098];
% xIonsInParam   = [  0.5479 -97.8176 -1.0219 -9.7033 0.0049 2.4240];
% xIonsAndMDStat = [    0.4004  -44.3595   -1.1232   -0.7356    0.0139    2.5034];
% xIonsAnd8Iter  = [   0.6260 -112.0291   -0.9296  -12.1084    0.0037    2.7546];
% 
% x = xIonsAnd8Iter;
[err,calc,ref,es,np,hb,disp,disp_slsl,disp_svsl,disp_svsv,cav,comb]=ObjectiveFromBEMCosmo(x);
save('RunCosmoBondiiFinal_4_HB.mat','mol_list','training_set','err','calc','ref','es','np','hb',...
    'disp','disp_slsl','disp_svsl','disp_svsv','cav','comb',...
     'disp_mob','cav_mob','np_mob','es_mob','np_SLIC','es_SLIC');
