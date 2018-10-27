% Path information
Home = getenv('HOME');
addpath(sprintf('%s/repos/pointbem',Home));
addpath(sprintf('%s/repos/panelbem',Home));
addpath(sprintf('%s/repos/testasymmetry',Home));
addpath(sprintf('%s/repos/testasymmetry/functions',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley',Home));
addpath(sprintf('%s/repos/testasymmetry/born',Home));

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet2" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants2 ProblemSet2 saveMemory writeLogfile logfileName
logfileName = 'two-region.out';
epsOut = 78.36;

curdir = pwd;
chdir('mobley');
ParamWatInfo = load('OptWater_thermo');
x = ParamWatInfo.xvec(3,:);
x_PB = [0 0 0 0 0 ParamWatInfo.xvec(3,6:7)];

chdir(curdir);

fid = fopen('two-region/two-region.csv','r');
Data = textscan(fid,'%s %s %f %f %f','delimiter',',');
fclose(fid);
mol_list1 = Data{1};
mol_list2 = Data{2};
all_surfAreas1 = Data{3};
all_surfAreas2 = Data{4};
dG_list = Data{5};
surfArea_list1 = all_surfAreas1;
surfArea_list2 = all_surfAreas2;
saveMemory = 1;
writeLogfile = 1;
epsIn1  =  1;
epsIn2  =  1;
Tbase = 300; mytemp=Tbase;
staticpotential = 0.0;
KelvinOffset = 273.15;
conv_factor = 332.112;
%staticpotential = 0.0; % this only affects charged molecules;

kappa = 0.0;  % should be zero, meaning non-ionic solutions!
UsefulConstants2 = struct('epsIn1',epsIn1,'epsIn2',epsIn2,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);

% here we define the actual params for the NLBC test

curdir = pwd;
for i=1
  dir=sprintf('%s/%s',curdir,'two-region');
  chdir(dir);
  pqrData1 = loadPqr('test1.pqr');
  pqrData2 = loadPqr('test2.pqr');
  pqrAll1{i} = pqrData1;
  pqrAll2{i} = pqrData2;
  srfFile1{i} = sprintf('%s/test1_2.srf',dir);
  srfFile2{i} = sprintf('%s/test2_2.srf',dir);
  chargeDist1{i} = pqrData1.q;
  chargeDist2{i} = pqrData2.q;
  referenceData{i} = dG_list(i);
  surfArea1{i} = surfArea_list1(i);
  surfArea2{i} = surfArea_list2(i);
  chdir(curdir);
  addProblemSA2(mol_list1{i},mol_list2{i},pqrAll1{i},pqrAll2{i},srfFile1{i},srfFile2{i},...
               chargeDist1{i},chargeDist2{i},referenceData{i},surfArea1{i},surfArea2{i});
end

%[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA2(x);
%SLICFileName = sprintf('RunTwoRegionSLIC');
%save(SLICFileName,'mol_list1','mol_list2','errfinal','calcE','refE','es','np');
[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA2(x_PB);
PBFileName = sprintf('RunTwoRegionPB');
save(PBFileName,'mol_list1','mol_list2','errfinal','calcE','refE','es','np');

   
