% Uncomment 68, 86, 89, 90 and 95 and comment out 69, 79-81, 87, 91, 92 and 96 
% ... if linear SA model is being used and viceversa for ES only
clear all
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
epsOut = 78.36; % from MNSol
KelvinOffset = 273.15;
conv_factor = 332.112;
staticpotential = 0.0; % this only affects charged molecules;
kappa = 0.0;  % should be zero, meaning non-ionic solutions!

% the staticpotential below should not be used any more, please check
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,'staticpotential',staticpotential);
     
fid = fopen('mnsol/water_es.csv','r'); 
Data = textscan(fid,'%s %f','delimiter',',');
fclose(fid);
mol_list = Data{1};
dG_list = Data{2};
fid = fopen('mnsol/mobley_sa_Chris.csv','r');
Data = textscan(fid,'%s %f %f %f %f %f %f','delimiter',',');
fclose(fid);
all_solutes = Data{1};
all_surfAreas = Data{2};
es_mobley = Data{6};
np_mobley = Data{7};
[m, index] = ismember(mol_list,all_solutes);
surfArea_list = all_surfAreas(index);
testset  = {'anthracene','cyclohexane','ethane','ethanol'...
    'methanol','n_heptane','n_hexane','n_octane',...
    'n_pentane','nitromethane','propan_1_ol','pyrene'};

curdir=pwd;
for i=1:length(testset)
  dir=sprintf('%s/Dropbox/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',getenv('HOME'),testset{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;
  foo = strcmp(mol_list,testset{i});
  index = find(foo);
  if length(index) ~= 1
    fprintf('error finding refdata!\n');
    keyboard
  end
  referenceData{i} = dG_list(index);
  surfArea{i} = surfArea_list(index);
  chdir(curdir);
  %addProblemSA(testset{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
  addProblem(testset{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i});
end


% The following script is specialized to this example.  We'll
% handle generating others.  Not complicated, but it's not self-explanatory.

x0 = [0.5 -60 -0.5  -0.5*tanh(- -0.5) 0 0 0];
lb = [0 -200 -100 -1  -0.1  x0(6)-0.001 x0(7)-0.001];
ub = [2 +200 +100 +1  +0.1  x0(6)+0.001 x0(7)+0.001];
x0 = x0(1:5);
lb = lb(1:5);
ub = ub(1:5);

options = optimoptions('lsqnonlin','MaxIter',8);
options = optimoptions(options,'Display', 'iter');

%y = @(x)ObjectiveFromBEMSA(x);
y = @(x)ObjectiveFromBEM(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
%[err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);
%[err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
err0 = ObjectiveFromBEM(x0);
err  = ObjectiveFromBEM(x);
rmse_train = rms(err);

%save('OptWaterSA','x','ref','calc','es','np','x0','calc0','es0','np0','rmse_train');
save('OptWater_es','x','err','rmse_train');

