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
global UsefulConstants ProblemSetNp saveMemory writeLogfile logfileName
saveMemory = 0;
writeLogfile = 0;
logfileName = 'junklogfile';

allData = readtable('../all_Bondii_data.csv'); 

%%
%%
%COSMO Dispersion atom types
%all_atom_types = {'br','c-sp','c-sp2','c-sp3','cl',...
%                  'f','h','i','n-sp','n-sp3',...
%                  'n-sp2','o-sp2','o-sp3-h',...
%                  'o-sp2-n','o-sp3','p','s'};
              
%COSMO H-bond atom types
allHbondTypes = {'n_amine','n_amide','n_nitro',...
                 'n_other','o_carbonyl','o_ester',...
                 'o_nitro','o_hydroxyl','fluorine',...
                 'h_oh','h_nh','h_other'};
             
epsOutWater = 78.34;% from mnsol Database
epsIn  =  1;
epsOut = 12.85;%Wakai05 and Krossing06  ;Data_eps.Var2(14);
KelvinOffset = 273.15;
conv_factor = 332.112;
staticpotential = 0.0; % this only affects charged molecules;
kappa = 0.0;  % should be zero, meaning non-ionic solutions!

% the staticpotential below should not be used any more, please check
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);
%             'all_atom_types',all_atom_types);
 

training_set  = {'4_bromophenol', 'ethanamide', 'teflurane', '4_chloroaniline',...
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
dG_np_SLIC = allData.np_SLIC; 
dG_es_SLIC = allData.es_SLIC; 
all_solutes = allData.solute;
mol_list = all_solutes;
solventAreas = allData{495,9:79};
solventATypes = allData{495,80:144};
solventHbData = allData{495,145:end};
solventVdWA = allData{495,11};
solventVdWV = allData{495,12};
soluteAreas = allData{:,9:79};
soluteATypes = allData{:,80:144};
soluteHbData = allData{:,145:end};
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
  soluteHbondData{i} = allData{index,145:end};
  referenceData{i} = dG_np_mob(index);
  solute_VdWA{i} = allData{index,11};
  solute_VdWV{i} = allData{index,12};
  solventAtomAreas{i} = solventAreas;
  solventAtomTypes{i} = {solventATypes};
  solventHbondData{i} = solventHbData;
  solvent_VdWA{i} = solventVdWA;
  solvent_VdWV{i} = solventVdWV;
  atom_vols{i} = allData{index,14};
  temp{i} = temperature;
  chdir(curdir);
  addProblemCosmoNp(training_set{i},chargeDist{i},referenceData{i},...
                  soluteAtomAreas{i},soluteAtomTypes{i},...
                  solute_VdWV{i},solute_VdWA{i},...
                  solventAtomAreas{i},solventAtomTypes{i},...
                  solvent_VdWV{i},solvent_VdWA{i},...
                  atom_vols{i},temp{i});
end


% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,id]=ismember(training_set,mol_list);
disp_mob = allData.disp_mobley(id); 
cav_mob = allData.cav_mobley(id); 
np_mob = allData.np_mobley(id); 
es_mob = allData.es_mobley(id); 
np_SLIC = allData.np_SLIC(id); 
es_SLIC= allData.es_SLIC(id);

x0_1 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.2 ... % disp
      12 ... %combinatorial x(36)
      0.5]; %cavity correction coeff

x0_2 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.2 ... % disp
      12 ... %combinatorial x(36)
      1.0]; %cavity correction coeff

x0_3 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.2 ... % disp
      12 ... %combinatorial x(36)
      1.5]; %cavity correction coeff

x0_4 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.2 ... % disp
      12 ... %combinatorial x(36)
      2.0]; %cavity correction coeff

x0_5 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.3 ... % disp
      12 ... %combinatorial x(36)
      0.5]; %cavity correction coeff

x0_6 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.3 ... % disp
      12 ... %combinatorial x(36)
      1.0]; %cavity correction coeff

x0_7 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.3 ... % disp
      12 ... %combinatorial x(36)
      1.5]; %cavity correction coeff

x0_8 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.3 ... % disp
      12 ... %combinatorial x(36)
      2.0]; %cavity correction coeff

x0_9 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.4 ... % disp
      12 ... %combinatorial x(36)
      0.5]; %cavity correction coeff

x0_10 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.4 ... % disp
      12 ... %combinatorial x(36)
      1.0]; %cavity correction coeff

x0_11 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.4 ... % disp
      12 ... %combinatorial x(36)
      1.5]; %cavity correction coeff

x0_12 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.4 ... % disp
      12 ... %combinatorial x(36)
      2.0]; %cavity correction coeff

x0_13 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.5 ... % disp
      12 ... %combinatorial x(36)
      0.5]; %cavity correction coeff

x0_14 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.5 ... % disp
      12 ... %combinatorial x(36)
      1.0]; %cavity correction coeff

x0_15 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.5 ... % disp
      12 ... %combinatorial x(36)
      1.5]; %cavity correction coeff

x0_16 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.5 ... % disp
      12 ... %combinatorial x(36)
      2.0]; %cavity correction coeff

x0_17 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.6 ... % disp
      12 ... %combinatorial x(36)
      0.5]; %cavity correction coeff

x0_18 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.6 ... % disp
      12 ... %combinatorial x(36)
      1.0]; %cavity correction coeff

x0_19 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.6 ... % disp
      12 ... %combinatorial x(36)
      1.5]; %cavity correction coeff

x0_20 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.6 ... % disp
      12 ... %combinatorial x(36)
      2.0]; %cavity correction coeff

x0_21 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.2 ... % disp
      12 ... %combinatorial x(36)
      2.5]; %cavity correction coeff

x0_22 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.3 ... % disp
      12 ... %combinatorial x(36)
      2.5]; %cavity correction coeff

x0_23 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.4 ... % disp
      12 ... %combinatorial x(36)
      2.5]; %cavity correction coeff

x0_24 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.5 ... % disp
      12 ... %combinatorial x(36)
      2.5]; %cavity correction coeff

x0_25 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.6 ... % disp
      12 ... %combinatorial x(36)
      2.5]; %cavity correction coeff

x0_26 = [64170 23743 37793 42479 42223 8157 0 106000 5161 ...  % disp x(6:23)
      36000 27800 7700 13000 11800 8200 56000 67000 0.7 ... % disp
      12 ... %combinatorial x(36)
      2.5]; %cavity correction coeff



ub = [100000 100000 100000 100000 100000 100000 100000 200000 100000 ...
      100000 100000 100000 100000 100000 100000 100000 100000 1 ...
      40 ...
      4];
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
      1 ...
      0.25];
      
options = optimoptions('lsqnonlin','MaxIter',16);
options = optimoptions(options,'Display', 'iter');

x0=x0_17;

y = @(x)ObjectiveFromBEMCosmoNp(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,np,disp,disp_slsl,disp_svsl,disp_svsv,cav,comb]=ObjectiveFromBEMCosmoNp(x);
[err0,calc0,ref0,np0,disp0,disp_slsl0,disp_svsl0,disp_svsv0,cav0,comb0]=ObjectiveFromBEMCosmoNp(x0);
rmse = rms(calc-ref);
save('OptCosmoBondiiNp_17.mat','x','training_set','mol_list','rmse','ref','calc','np','disp',...
     'disp_slsl','disp_svsl','disp_svsv','comb','cav',...
     'disp_mob','cav_mob','np_mob','es_mob','np_SLIC',...
     'x0','calc0','np0','disp0', 'disp_slsl0','disp_svsl0',...
     'disp_svsv0','comb0','cav0','epsOut');
