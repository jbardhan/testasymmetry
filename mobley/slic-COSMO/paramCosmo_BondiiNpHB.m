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
% variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSetNp saveMemory writeLogfile logfileName
saveMemory = 0;
writeLogfile = 0;
logfileName = 'junklogfile';

allData = readtable('all_Bondii_data.csv'); 

%%
%%
%COSMO Dispersion atom types
%all_atom_types = {'br','c-sp','c-sp2','c-sp3','cl',...
%     'f','h','i','n-sp','n-sp3',...
%     'n-sp2','o-sp2','o-sp3-h',...
%     'o-sp2-n','o-sp3','p','s'};
    
%COSMO H-bond atom types
allHbondTypes = {'n_amine','n_amide','n_nitro',...
     'n_other','o_carbonyl','o_ester',...
     'o_nitro','o_hydroxyl','fluorine',...
     'h_oh','h_nh','h_other'};
    
epsOutWater = 78.34;% from mnsol Database
epsIn = 1;
epsOut = 12.85;%Wakai05 and Krossing06 ;Data_eps.Var2(14);
KelvinOffset = 273.15;
conv_factor = 332.112;
staticpotential = 0.0; % this only affects charged molecules;
kappa = 0.0; % should be zero, meaning non-ionic solutions!

% the staticpotential below should not be used any more, please check
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);
%    'all_atom_types',all_atom_types);
 

training_set = {'4_bromophenol', 'ethanamide', 'teflurane', '4_chloroaniline',...
   '2_methylpropane', '222_trifluoroethanol', '2_fluorophenol',...
   '2_iodopropane', 'iodobenzene', '1_nitropentane', '3_cyanophenol',...
   'pyridine','4_nitroaniline','14_dioxane','acetic_acid','butan_1_ol',...
   'methyl_acetate','propanone','triethyl_phosphate','trimethyl_phosphate',...
   'methanethiol','dimethyl_sulfate','piperidine','ethylamine','N_methylacetamide',...
   'nitromethane','nonanal','benzaldehyde','methanol',...
   '3_methyl_1h_indole','anthracene','124_trimethylbenzene','2_naphthylamine'...
   ,'4_formylpyridine','cyclohexylamine','dimethyl_sulfide','hex_1_ene',...
   'n_butanethiol','naphthalene','33_dimethylbutan_2_one','333_trimethoxypropionitrile',...
   'chloroethane','diethyl_sulfide','ethene','imidazole',...
   'methyl_octanoate','n_octane','n_propylbenzene',...
   'p_cresol','propanoic_acid','tetrahydropyran','trichloroethene',...
   '2_methoxyaniline','2_methylhexane','2_nitropropane','26_dimethylpyridine',...
   'benzene','but_1_ene','but_1_yne','m_xylene','methane',...
   'n_pentylamine','p_dibromobenzene'};
  
dG_list = allData.dG_expt; 
dG_disp_mob = allData.disp_mobley; 
dG_cav_mob = allData.cav_mobley; 
dG_np_mob = allData.np_mobley; 
dG_es_mob = allData.es_mobley; 
dG_es_best = allData.best_elec; 
dG_np_SLIC = allData.np_SLIC; 
dG_es_SLIC = allData.es_SLIC; 
all_solutes = allData.solute;
mol_list = all_solutes;
solventAreas = allData{495,9:79};
solventATypes = allData{495,80:144};
solventHbData = allData{495,145:158};
solventVdWA = allData{495,11};
solventVdWV = allData{495,12};
soluteAreas = allData{:,9:79};
soluteATypes = allData{:,80:144};
soluteHbData = allData{:,145:158};
soluteVdWA = allData{:,11};
soluteVdWV = allData{:,12};
temperature = 24.85 + KelvinOffset;
curdir=pwd;

for i=1:length(training_set)
  dir=sprintf('%s/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',dropbox_path,training_set{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;%chargeDistribution;
  foo = strcmp(mol_list,training_set{i});
  index = find(foo);
  if length(index) ~= 1
    fprintf('error finding refdata!\n');
    keyboard
  end
  soluteAtomAreas{i} = allData{index,9:79};
  soluteAtomTypes{i} = {allData{index,80:144}};
  soluteHbondData{i} = allData{index,145:158};
  referenceData{i} = dG_list(index) - dG_es_best(index);
  solute_VdWA{i} = allData{index,11};
  solute_VdWV{i} = allData{index,12};
  solventAtomAreas{i} = solventAreas;
  solventAtomTypes{i} = {solventATypes};
  solventHbondData{i} = solventHbData;
  solvent_VdWA{i} = solventVdWA;
  solvent_VdWV{i} = solventVdWV;
  atom_vols{i} = allData{index,14};
  newHB{i}=0; 
     
  temp{i} = temperature;
  chdir(curdir);
  addProblemCosmoNpHB(training_set{i},chargeDist{i},referenceData{i},...
      soluteAtomAreas{i},soluteAtomTypes{i},soluteHbondData{i},...
      solute_VdWV{i},solute_VdWA{i},...
      solventAtomAreas{i},solventAtomTypes{i},solventHbondData{i},...
      solvent_VdWV{i},solvent_VdWA{i},...
      atom_vols{i},temp{i},newHB{i});
end


% The following script is specialized to this example. We'll
% handle generating others. Not complicated, but it's not self-explanatory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,id]=ismember(training_set,mol_list);
disp_mob = allData.disp_mobley(id); 
cav_mob = allData.cav_mobley(id); 
np_mob = allData.np_mobley(id); 
es_mob = allData.es_mobley(id); 
np_SLIC = allData.np_SLIC(id); 
es_SLIC= allData.es_SLIC(id);

x0_1 = [0.81 1.02 1.05 1.74 0.55 0.13 0.00 0.49 0.78 ... % disp x(6:14)  
0.37 1.54 0.06 0.38 0.72 0.55 0.35 1.10 0.34 ... % disp x(15:23)  
0.59 0.09 0.22 0.08 0.05 0.07 0.35 0.83 0.11 1.19 3.74 0.04 ... % hb x(24:35)
1.44 ... % comb x(36)      
2.69 ]; % cav x(37)      

x0_2 =[ 0.23 0.47 0.57 0.38 0.22 0.74 0.74 0.61 0.30 ... % disp x(6:14)    
   0.36 0.63 0.40 0.60 0.97 0.27 0.04 0.78 0.83 ... % disp x(15:23)    
   0.60 0.43 0.74 0.28 0.30 0.87 0.15 0.59 0.52 0.33 0.26 0.35 ... % hb x(24:35) 
   5.96 ... % comb x(36)            
   0.65 ]; % cav x(37)            
x0_3 =[ 0.57 0.83 0.78 0.29 0.45 0.88 0.64 0.07 0.18 ... % disp x(6:14)    
   0.34 0.24 0.58 0.04 0.25 0.19 0.95 0.95 0.44 ... % disp x(15:23)    
   0.79 0.22 0.55 0.48 0.39 0.75 0.69 0.17 0.67 0.50 0.54 0.34 ... % hb x(24:35) 
   1.67 ... % comb x(36)            
   0.19 ]; % cav x(37)            
x0_4 =[ 0.27 0.53 0.96 0.21 0.10 0.11 0.45 0.84 0.01 ... % disp x(6:14)    
   0.33 0.65 0.42 0.97 0.96 0.64 0.29 0.69 0.84 ... % disp x(15:23)    
   0.81 0.34 0.09 0.47 0.66 0.67 0.36 0.03 0.82 0.32 0.85 0.24 ... % hb x(24:35) 
   6.64 ... % comb x(36)            
   0.25 ]; % cav x(37)            
x0_5 =[ 0.31 0.84 0.94 0.28 0.68 0.18 0.27 0.92 0.15 ... % disp x(6:14)    
   0.35 0.74 0.59 0.53 0.88 0.50 0.74 0.64 0.11 ... % disp x(15:23)    
   0.93 0.93 0.79 0.39 0.15 0.98 0.34 0.64 0.87 0.35 0.43 0.79 ... % hb x(24:35) 
   0.54 ... % comb x(36)            
   0.42 ]; % cav x(37)            
x0_6 =[ 0.09 0.57 0.01 0.33 0.08 0.39 0.63 0.29 0.54 ... % disp x(6:14)    
   0.11 0.66 0.00 0.92 0.22 0.78 0.90 0.82 0.82 ... % disp x(15:23)    
   0.64 0.43 0.36 0.06 0.24 0.55 0.74 0.39 0.24 0.32 0.94 0.76 ... % hb x(24:35) 
   2.34 ... % comb x(36)            
   0.44 ]; % cav x(37)            
x0_7 =[ 0.83 0.80 0.48 0.03 0.00 0.59 0.70 0.03 0.77 ... % disp x(6:14)    
   0.93 0.19 0.84 0.71 0.49 0.99 0.81 0.69 0.46 ... % disp x(15:23)    
   0.28 0.80 0.66 0.76 0.89 0.52 0.36 0.23 0.23 0.14 0.41 0.25 ... % hb x(24:35) 
   9.96 ... % comb x(36)            
   0.37 ]; % cav x(37)            
x0_8 =[ 0.55 0.97 0.83 0.85 0.24 0.59 0.73 0.54 0.59 ... % disp x(6:14)    
   0.60 0.64 0.75 0.56 0.66 0.95 0.15 0.15 0.34 ... % disp x(15:23)    
   0.55 0.59 0.29 0.04 0.35 0.21 0.84 0.14 0.19 0.76 0.07 0.16 ... % hb x(24:35) 
   6.72 ... % comb x(36)            
   0.48 ]; % cav x(37)            
x0_9 =[ 0.50 0.53 0.31 0.13 0.35 0.51 0.70 0.60 0.56 ... % disp x(6:14)    
   0.41 0.34 0.55 0.88 0.78 0.35 0.13 0.99 0.54 ... % disp x(15:23)    
   0.14 0.20 0.66 0.74 0.94 0.50 0.70 0.51 0.91 0.70 0.34 0.25 ... % hb x(24:35) 
   5.92 ... % comb x(36)            
   0.82 ]; % cav x(37)            
x0_10 =[ 0.70 0.27 0.28 0.85 0.43 0.84 0.42 0.39 0.33 ... % disp x(6:14)    
   0.50 1.00 0.31 0.97 0.73 0.58 0.40 0.47 0.12 ... % disp x(15:23)    
   0.66 0.47 0.45 0.30 0.19 0.14 0.62 0.32 0.38 0.55 0.36 0.64 ... % hb x(24:35) 
   8.24 ... % comb x(36)            
   0.20 ]; % cav x(37)            
x0_11 =[ 0.22 0.23 0.73 0.04 0.03 0.48 0.30 0.36 0.18 ... % disp x(6:14)    
   0.63 0.43 0.68 0.91 0.82 0.89 0.85 0.48 0.03 ... % disp x(15:23)    
   0.99 0.89 0.48 0.91 0.01 0.04 0.20 0.85 0.06 0.63 0.74 0.66 ... % hb x(24:35) 
   4.57 ... % comb x(36)            
   0.71 ]; % cav x(37)            
x0_12 =[ 0.80 0.89 0.04 0.11 0.79 0.27 0.65 0.51 0.92 ... % disp x(6:14)    
   0.12 0.99 0.17 0.22 0.04 0.48 0.10 0.68 0.74 ... % disp x(15:23)    
   0.84 0.13 0.09 0.66 1.00 0.58 0.34 0.33 0.73 0.83 0.09 0.66 ... % hb x(24:35) 
   3.36 ... % comb x(36)            
   0.17 ]; % cav x(37)            

x0=x0_1;

ub = [4 4 4 4 4 4 4 4 4 ...
  4 4 4 4 4 4 4 4 1 ...
  2 2 2 2 2 2 2 2 2 2 2 2 ...
  20 ...
  4];
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
  0 0 0 0 0 0 0 0 0 0 0 0 ...
  1 ...
  0.25];
 
options = optimoptions('lsqnonlin','MaxIter',30);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEMCosmoNpHB(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np,hb,disp,disp_slsl,disp_svsl,disp_svsv,cav,comb]=ObjectiveFromBEMCosmoNpHB(x);
[err0,calc0,ref0,es0,np0,hb0,disp0,disp_slsl0,disp_svsl0,disp_svsv0,cav0,comb0]=ObjectiveFromBEMCosmoNpHB(x0);
[~,id]=ismember(training_set,mol_list);
disp_mob = allData.disp_mobley(id); 
cav_mob = allData.cav_mobley(id); 
np_mob = allData.np_mobley(id); 
es_mob = allData.es_mobley(id); 
np_SLIC = allData.np_SLIC(id); 
es_SLIC= allData.es_SLIC(id);
rmse = rms(calc-ref);
save('OptCosmoBondiinphb_1.mat','x','training_set','mol_list','rmse','ref','calc','es','np','hb','disp',...
  'disp_slsl','disp_svsl','disp_svsv','comb','cav',...
  'disp_mob','cav_mob','np_mob','es_mob','np_SLIC',...
  'x0','calc0','es0','np0','hb0','disp0', 'disp_slsl0','disp_svsl0',...
  'disp_svsv0','comb0','cav0','epsOut','x');




