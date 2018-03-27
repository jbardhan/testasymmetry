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
epsOut = 78.36; % from MNSol
 
mytemp=Tbase;
KelvinOffset = 273.15;
conv_factor = 332.112;
staticpotential = 0.0; % this only affects charged molecules;
kappa = 0.0;  % should be zero, meaning non-ionic solutions!
 
 
% the staticpotential below should not be used any more, please check
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
             kappa,'conv_factor',conv_factor,...
             'staticpotential',staticpotential);
 
testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane','Li','Na','K','Rb','Cs','Cl','Br','I'};  % test set without florine
      
       
fid = fopen('~/repos/testasymmetry/mobley/mnsol/mobley_dG_AND_sa_and_vol_fixed.csv','r');
Data = textscan(fid,'%s %f  %f  %f  %f  %f  %f  %f','delimiter',',');
fclose(fid);
mol_list = Data{1};
dG_list = Data{2};
all_surfAreas = Data{3};
[m, index] = ismember(testset,mol_list);
surfArea_list = all_surfAreas(index);
         
 
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
  referenceData{i} = dG_list(i);
  surfArea{i} = surfArea_list(i);
  chdir(curdir);
  addProblemSA(testset{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end
 
 
% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
x0 = [0.0001 0.0001 0.0001   0.0001     0.0001 0.00001 0.00001];
lb = [-0.001 -0.001 -0.001 -0.001  -0.001  -0.1  -2];
ub = [+0.001 +0.001 +0.001 +0.001  +0.001  +0.1  +2];
 
options = optimoptions('lsqnonlin','MaxIter',8);
options = optimoptions(options,'Display', 'iter');
 
y = @(x)ObjectiveFromBEMSA(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
[err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);
 
save('OptWater_SLIC_Off','x','ref','calc','es','np','x0','calc0','es0','np0');

