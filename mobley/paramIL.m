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
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
saveMemory = 0;
writeLogfile = 0;
logfileName = 'junklogfile';

Data_eps = readtable('reference-data/eps_data.csv'); 
% Dielectric constant of ILs at 298K
% check reference-data/RTIL.xls for more details on dG and epsOut of ILs
% Reference articles are available at Dropbox/lab/
% ref1 = Wakai05, ref2 = Rybinska14, ref3 = Hunger09, ref4 = Huang11 
% ref5 = Nakamura10
% Data_eps.Var2 = dielectric constant 
% Var2.(1) emimtf2n_ref1     (2) emimtf2n_ref2      (3) emimtf2n_ref4  
%      (4) bmimtf2n_ref2     (5) bmimtf2n_ref2
%      (5) bmimbf4_ref1      (6) bmimbf4_ref2       (7) bmimbf4_ref3 
%      (8) bmimpf6_ref1      (9) bmimpf6_ref2       (10)bmimpf6_ref3   
%      (11)bmimtfo_ref2      (12)bmimtfo_ref4    
%      (13)bmimcl_ref3    
%      (14)omimtf2n_ref2
%
% example: Data_eps.Var2(11) = dielectric constant of bmimtfo from ref2

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
     
Data_p = readtable('reference-data/paluch_all.csv');
Data_l = readtable('reference-data/latif_all.csv');

% When using Paluch data: Data=Data_p
% Paluch Data: L1:water L2:octanol L3:emimtf2n L4:bmimtf2n L5:omimtf2n 
%
%                  L1  L2  L3  L4  L5
%  dG_298   Data.Var2 | 3 | 4 | 5 | 6
%     328           7 | 8 | 9 |10 |11
%     358          12 |13 |14 |15 |16
%  dH_328          17 |18 |19 |20 |21
%-TdS_328          22 |23 |24 |25 |26
%
% example: Data.Var13 = dG @ 358K for octanol(L2)
%
% When using Latif data: Data=Data_l
% Latif Data: L1:bmimbf4 L2:bmimtf2n L3:bmimpf6 L4:bmimtfo L5:bmimcl
%
%                  L1  L2  L3  L4  L5
%  dG_298   Data.Var2 | 3 | 4 | 5 | 6
%
% example: Data.Var4 = dG @ 298K for bmimpf6(L3)
% Data.Var1 is the solute list in both cases

Data = Data_p;
mol_list = Data.Var1;
dG_list = Data.Var6; %omimtf2n_paluch
OptFileName = sprintf('Opt_%s.mat','omimtf2n_paluch_ref2'); %16.8 Nkamura10

fid = fopen('~/repos/testasymmetry/mobley/mnsol/mobley_sa.csv','r');
Data = textscan(fid,'%s %f','delimiter',',');
fclose(fid);
all_solutes = Data{1};
all_surfAreas = Data{2};
[m, index] = ismember(mol_list,all_solutes);
surfArea_list = all_surfAreas(index);

testset = {'methane','propane','n_butane','2_methylpropane','toluene',...
           'methanol','ethanol','p_cresol','methanethiol',...
           'methyl_ethyl_sulfide','ethanamide','3_methyl_1h_indole'};

curdir=pwd;
for i=1:length(testset)
  dir=sprintf('%s/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',dropbox_path,testset{i});
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


% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0.5 -60 -0.5   -0.5*tanh(- -0.5)     0 -0.03 1.6];
lb = [-2 -200 -100 -20  -0.1  -0.1  -4];
ub = [+2 +200 +100 +20  +0.1  +0.1  +4];

options = optimoptions('lsqnonlin','MaxIter',8);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEMSA(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
[err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);
rmse = rms(calc-ref);
save(OptFileName,'x','rmse','ref','calc','es','np','x0','calc0','es0','np0','epsOut');

