% Path information
clear all; clear global
DataStructure = load('OptFiles.mat');
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
logfileName = 'chloroform.out';
epsOut = 2.0402; % from MNSol

fid = fopen('mnsol/chloroform.csv','r'); 
Data = textscan(fid,'%s %f %f','delimiter',',');
fclose(fid);
mol_list = Data{1};
dG_list = Data{2};

fid = fopen('mnsol/mobley_sa.csv','r');
Data = textscan(fid,'%s %f','delimiter',',');
fclose(fid);
all_solutes = Data{1};
all_surfAreas = Data{2};
[m, index] = ismember(mol_list,all_solutes);
surfArea_list = all_surfAreas(index);
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

curdir = pwd;
for i=1:length(mol_list)
  dir=sprintf('%s/Dropbox/lab/Projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',getenv('HOME'),mol_list{i});
  chdir(dir);
  pqrData = loadPqr('test.pqr');
  pqrAll{i} = pqrData;
  srfFile{i} = sprintf('%s/test_2.srf',dir);
  chargeDist{i} = pqrData.q;
  referenceData{i} = dG_list(i);
  surfArea{i} = surfArea_list(i);
  chdir(curdir);
  addProblemSA(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end

[n,m] = ismember({DataStructure.LoadData.Solvent},'Chloroform');
ParamFiles = DataStructure.LoadData(n);
for i = 1:length(ParamFiles)
    ParamInfo = load(ParamFiles(i).loadFile);
    testSets(:,i) = ParamInfo.testset;
    x = ParamInfo.x;
    [errfinal(:,i),calcE(:,i),refE(:,i),es(:,i),np(:,i)]=ObjectiveFromBEMSA(x);
end

save('RunChloroform','errfinal','calcE','refE','es','np','testSets');
