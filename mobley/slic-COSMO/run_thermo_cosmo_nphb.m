close all
clear all
clc

% Path information

Home = getenv('HOME');
repo_path=sprintf('%s/repos',Home);
dropbox_path=sprintf('%s/Dropbox',Home);
addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));                    
addpath(sprintf('%s/testasymmetry/mobley/reference-data',repo_path));  
addpath(sprintf('%s/repos/testasymmetry/mobley/slic-COSMO',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley/slic-COSMO/tests',Home));                  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%                                Set these values before running the code                      %%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_min= 4.85; % lower bound of the temperature interval
temp_max=44.85; % upper bound in the temperature interval
tempdiv=5; % number of temperatures that we calculate solvation free energies 
temp_C=linspace(temp_min,temp_max,tempdiv); % temperatures in Celcius            
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the optimized parameters from file
optimizedParamFile = 'OptWater_Cosmo_thermo_nphb.mat';
xES =  [0.95	-22.1	-1.17	-1.19	-0.1
        0.95	-22.1	-1.17	-1.19	-0.1
        0.95	-22.1	-1.17	-1.19	-0.1
        0.95	-22.1	-1.17	-1.19	-0.1
        0.95	-22.1	-1.17	-1.19	-0.1];
ParamWatInfo=load(optimizedParamFile);
x = [xES ParamWatInfo.xvec];
allData = readtable('all_Bondii_data.csv'); 
%allData = readtable('dG_data.csv'); 
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
%temperature = 24.85 + KelvinOffset;
curdir=pwd;

for j=1:tempdiv
    clear global
    loadConstants
    convertKJtoKcal = 1/joulesPerCalorie;
    global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
    saveMemory = 1;
    writeLogfile = 1;
    logfileName = 'junklogfile';
    epsIn  =  1; % dielectric constant of the solutes
    KelvinOffset = 273.15;
    conv_factor = 332.112;
    staticpotential = 0.0; % this only affects charged molecules;
    kappa = 0.0;  % should be zero, meaning non-ionic solutions!
    epsOut = (-1.410e-6)*temp_C(j)^3+(9.398e-4)*temp_C(j)^2-0.40008*temp_C(j)+87.740;  
        % temperature dependence of the dielectric constant of Water T in C
    %dG_list = dGdata{:,j}; 
    UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa',kappa,...
                             'conv_factor',conv_factor,'staticpotential',staticpotential);
    
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
      temp{i} = temp_C(j)+KelvinOffset;
      referenceData{i} = dG_list(i);
      newHB{i}=5;
      chdir(curdir);
      addProblemCosmo(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},...
                      soluteAtomAreas{i},soluteAtomTypes{i},soluteHbondData{i},...
                      solute_VdWV{i},solute_VdWA{i},...
                      solventAtomAreas{i},solventAtomTypes{i},solventHbondData{i},...
                      solvent_VdWV{i},solvent_VdWA{i},...
                      atom_vols{i},temp{i},newHB{i});
    end

    [err(j,:),calcE(j,:),refE(j,:),esE(j,:),npE(j,:),...
     hbE(j,:),dispE(j,:),disp_slsl(j,:),disp_svsl(j,:),...
     disp_svsv(j,:),cavE(j,:),combE(j,:)]=ObjectiveFromBEMCosmo(x(j,:));
end
save('RunWater_Cosmo_nphb.mat','mol_list','x',...
     'err','calcE','refE','esE','npE','hbE',...
     'dispE','cavE','combE','temp_C');

RunWater=load('RunWater_Cosmo_nphb.mat');


dGfunc=struct(); % structure that has the information of the solutes and the linear function that 
                 % fits the calculated values of dG at different temperatures
calcE = RunWater.calcE; % calculated values for dG at different temperatures
refE = RunWater.refE; % calculated values for dG at different temperatures

    
esE=RunWater.esE;
npE=RunWater.npE;
cavE=RunWater.cavE;
combE=RunWater.combE;
dispE=RunWater.dispE;
hbE=RunWater.hbE;
x = RunWater.x;
TEMP=RunWater.temp_C';
temp_K=TEMP+273.15;
[m,index]=ismember(24.85,TEMP);
mol_list=RunWater.mol_list;

for i=1:length(mol_list)
    f = @(R) (R(1)-R(2)*(temp_K-298)+R(3)*((temp_K-298)-temp_K.*log(temp_K./298)))-calcE(:,i);
    R0=[refE(index,i),1,1];
    options=optimoptions('lsqnonlin','StepTolerance',1e-6);
    options=optimoptions(options,'OptimalityTolerance',1e-6);
    options=optimoptions(options,'FunctionTolerance',1e-6);
    [R,resnorm,residual,exitflag,output]=lsqnonlin(f,R0,[],[],options);
    dGfunc(i).name=mol_list(i); 
    dGfunc(i).dg=R(1);
    dGfunc(i).ds=R(2);
    dGfunc(i).cp=R(3); 
    dsvec(i)=dGfunc(i).ds*1000;
    cpvec(i)=dGfunc(i).cp*1000;
    resnorm(i)=resnorm;
    exitflag(i)=exitflag;
    output(i)=output;
end
disp_mob = allData.disp_mobley; 
cav_mob = allData.cav_mobley; 
np_mob = allData.np_mobley; 
es_mob = allData.es_mobley; 
np_SLIC = allData.np_SLIC; 
es_SLIC= allData.es_SLIC;
 
output_name='RunWater_thermo_Cosmo_nphb.mat';
save(output_name,'mol_list','x',...
     'calcE','refE','esE','npE','hbE',...
     'disp_mob','cav_mob','np_mob','es_mob','np_SLIC','es_SLIC',...
     'dispE','cavE','combE','dGfunc','dsvec','cpvec');



    