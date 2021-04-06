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
addpath(sprintf('%s/repos/testasymmetry/born',Home));

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSetNp" for COSMO nonpolar calculation
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSetNp saveMemory writeLogfile logfileName
logfileName = 'water.out';
epsOut = 78.36;
allData = readtable('all_Bondii_data.csv'); 

ParamWatInfo = load('OptCosmoBondiinphb_1.mat');
x = ParamWatInfo.x;
training_set = ParamWatInfo.training_set;
mol_list = allData.solute;
dG_list = allData.dG_expt;
disp_mob = allData.disp_mobley; 
cav_mob = allData.cav_mobley; 
np_mob = allData.np_mobley; 
np_SLIC = allData.np_SLIC; 

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
         
% here we define the actual params for the NLBC test
solventAreas = allData{495,9:79};
solventATypes = allData{495,80:144};
solventHbData = allData{495,145:158};
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
  solute_VdWA{i} = allData{i,11};
  solute_VdWV{i} = allData{i,12};
  soluteHbondData{i} = allData{i,145:158};
  solventAtomAreas{i} = solventAreas;
  solventAtomTypes{i} = {solventATypes};
  solventHbondData{i} = solventHbData;
  solvent_VdWA{i} = solventVdWA;
  solvent_VdWV{i} = solventVdWV;
  atom_vols{i} = allData{i,14};
  temp{i} = temperature;
  referenceData{i} = np_mob(i);
  newHB{i}=0;
  chdir(curdir);
  addProblemCosmoNpHB(training_set{i},chargeDist{i},referenceData{i},...
                  soluteAtomAreas{i},soluteAtomTypes{i},soluteHbondData{i},...
                  solute_VdWV{i},solute_VdWA{i},...
                  solventAtomAreas{i},solventAtomTypes{i},solventHbondData{i},...
                  solvent_VdWV{i},solvent_VdWA{i},...
                  atom_vols{i},temp{i},newHB{i});



end
% xNoIonsInParam = [  1.8570  -75.6855   -2.9884    8.7302    0.0113    0.4098];
% xIonsInParam   = [  0.5479 -97.8176 -1.0219 -9.7033 0.0049 2.4240];
% xIonsAndMDStat = [    0.4004  -44.3595   -1.1232   -0.7356    0.0139    2.5034];
% xIonsAnd8Iter  = [   0.6260 -112.0291   -0.9296  -12.1084    0.0037    2.7546];
% 
% x = xIonsAnd8Iter;
[err,calc,ref,es,np,hb,disp,disp_slsl,disp_svsl,disp_svsv,cav,comb]=ObjectiveFromBEMCosmoNp(x);
rmse = rms(np - np_mob);
save('RunCosmoBondiinphb_1.mat','mol_list','training_set','err','calc','ref','np',...
    'disp','disp_slsl','disp_svsl','disp_svsv','cav','comb',...
    'disp_mob','cav_mob','np_mob','np_SLIC');
