close all
clear all
clc
% Path information

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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%                                Set these values before running the code                      %%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ionflag=0; % ionflag=0 : ions data are not included in the testset 
           % ionflag=1 : ions data are included in the testset

temp_min= 4.85; % lower bound of the temperature interval
temp_max=44.85; % upper bound in the temperature interval
tempdiv=5; % number of temperatures that we calculate solvation free energies 
temp_C=linspace(temp_min,temp_max,tempdiv); % temperatures in Celcius    
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Note: we calculate everything at 298 K which is equal to 24.85C, We
%%%% will calculate \Delta G of ions which is available at 25C at 24.85 C
%%%% to use them in our objective function. But in the postprocessing
%%%% script (whaterthermo.m) we will compare 
                   
% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.

for j=1:tempdiv
    clear global
    loadConstants
    convertKJtoKcal = 1/joulesPerCalorie;
    global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
    saveMemory = 0;
    writeLogfile = 0;
    logfileName = 'junklogfile';
    epsIn  =  1; % dielectric constant of the solutes
    KelvinOffset = 273.15;
    conv_factor = 332.112;
    staticpotential = 0.0; % this only affects charged molecules;
    kappa = 0.0;  % should be zero, meaning non-ionic solutions!
    epsOut = (-1.410e-6)*temp_C(j)^3+(9.398e-4)*temp_C(j)^2-0.40008*temp_C(j)+87.740;  
    % temperature dependence of the dielectric constant of Water T in C
    % from Mollerup15 (Modeling the permittivity of electrolyte solutions)
    allHbondTypes = {'n_amine','n_amide','n_nitro',...
                 'n_other','o_carbonyl','o_ester',...
                 'o_nitro','o_hydroxyl','fluorine',...
                 'h_oh','h_nh','h_other'};
    UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa',kappa,...
                             'conv_factor',conv_factor,'staticpotential',staticpotential);
                         
    training_set  = {'heptan_4_one','pentan_1_ol','morpholine','14_dioxane','tetrahydrofuran','pentan_1_ol','4_methylpentan_2_one',...
                    'heptan_2_one','piperidine','methane','cyclopentanol','diethyl_disulfide','buta_13_diene','butanone','n_propylbenzene','chloroethylene','morpholine',...
                    'pentan_1_ol','octan_1_ol','n_propylamine','4_methylpentan_2_one','2_methylpropan_1_ol','heptan_1_ol','n_butane','2_methylpropan_2_ol','pentan_3_one',...
                    'pyrene','cyclopropane','2_methylpyridine','morpholine','2_methylbutan_2_ol','diethyl_disulfide','cyclopentanol','hexan_3_ol','4_methylpyridine',...
                    '3_nitrophenol','p_cresol','methane','2_butoxyethanol','cyclopropane','octan_1_ol','methane','n_hexane','2_methylpropane','12_dimethoxyethane',...
                    'pentan_3_one','n_propylbenzene','butan_1_ol','pyridine','p_cresol','dimethyl_sulfide','ethylbenzene','prop_2_en_1_ol','12_dimethoxyethane','toluene',...
                    'butanone','26_dimethylpyridine','pentan_1_ol','ethanamide','cyclohexane','N_methylmorpholine','2_propoxyethanol','ethene','2_methylpyridine','buta_13_diene','propene'};
                
    thermo_data = readtable('thermo_full.csv'); 
    allData = readtable('all_Bondii_data.csv'); 
    
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
    t_ref_aca=24.85;
    expt_mol_ref = thermo_data.Var1;
    dG_ref = thermo_data.Var2;
    dS_ref = thermo_data.Var3;
    Cp_ref = thermo_data.Var4;
    
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
      solute_VdWA{i} = allData{index,11};
      solute_VdWV{i} = allData{index,12};
      solventAtomAreas{i} = solventAreas;
      solventAtomTypes{i} = {solventATypes};
      solventHbondData{i} = solventHbData;
      solvent_VdWA{i} = solventVdWA;
      solvent_VdWV{i} = solventVdWV;
      atom_vols{i} = allData{index,14};
      temp{i} = temp_C(j)+KelvinOffset;
      newHB{i}=5;
      expt_foo = strcmp(expt_mol_ref,training_set{i});
      expt_index = find(expt_foo);
      if length(expt_index) ~= 1
        fprintf('error finding refdata!\n');
        keyboard
      end
      referenceData{i} = dG_ref(expt_index) - dS_ref(expt_index)*(temp_C(j)-t_ref_aca) + ...
      Cp_ref(expt_index) * ((temp_C(j) - t_ref_aca) - ...
      (temp_C(j) + KelvinOffset) * log(((temp_C(j) + KelvinOffset))/((t_ref_aca + KelvinOffset))));
      
      chdir(curdir);
      addProblemCosmo(training_set{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},...
                      soluteAtomAreas{i},soluteAtomTypes{i},soluteHbondData{i},...
                      solute_VdWV{i},solute_VdWA{i},...
                      solventAtomAreas{i},solventAtomTypes{i},solventHbondData{i},...
                      solvent_VdWV{i},solvent_VdWA{i},...
                      atom_vols{i},temp{i},newHB{i});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    initialGuessInfo = load('RunCosmo_final_3.mat');
    
    x0 = initialGuessInfo.x;
       
    ub = [+2 +200 +100 +20  +0.1 ...
          100000 100000 100000 100000 100000 100000 100000 200000 100000 ...
          100000 100000 100000 100000 100000 100000 100000 100000 1 ...
          2 2 2 2 2 2 2 2 2 2 2 2 ...
          40 ...
          2];
    lb = [-2 -200 -100 -20  -0.1 ...
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
          0 0 0 0 0 0 0 0 0 0 0 0 ...
          0 ...
          0];

    options = optimoptions('lsqnonlin','MaxIter',8);
    options = optimoptions(options,'Display', 'iter');

    y = @(x)ObjectiveFromBEMCosmo(x);
    [x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
    [err,calc,ref,es,np,hb,disp,disp_slsl,disp_svsl,disp_svsv,cav,comb]=ObjectiveFromBEMCosmo(x);
    [err0,calc0,ref0,es0,np0,hb0,disp0,disp_slsl0,disp_svsl0,disp_svsv0,cav0,comb0]=ObjectiveFromBEMCosmo(x0);
    [~,id]=ismember(training_set,mol_list);
    
    xvec(j,:)=x;
    refvec(j,:)=ref;
    calcvec(j,:)=calc;
    esvec(j,:)=es;
    npvec(j,:)=np;
    dispvec(j,:)=disp;
    cavvec(j,:)=cav;
    hbvec(j,:)=hb;
    xvec0(j,:)=x0;
    refvec0(j,:)=ref0;
    calcvec0(j,:)=calc0;
    esvec0(j,:)=es0;
    npvec0(j,:)=np0;
    dispvec0(j,:)=disp0;
    cavvec0(j,:)=cav0;
    hbvec0(j,:)=hb0;
    rmse(j,:)=rms(ref-calc);
    

end
save('OptWater_Cosmo_thermo','xvec','refvec','calcvec','esvec','npvec','dispvec','cavvec','combvec','hbvec',...
    'xvec0','refvec0','calcvec0','esvec0','npvec0','dispvec0','cavvec0','combvec0','hbvec0',...
    'x0vec','calc0vec','es0vec','np0vec','testset','dS','Cp','tempvec','ionflag','aca_num','ion_num','t_ref_aca','t_ref_ion','rmse');

