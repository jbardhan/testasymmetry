% Path information
Home = getenv('HOME');
addpath(sprintf('%s/repos/pointbem',Home));
addpath(sprintf('%s/repos/panelbem',Home));
addpath(sprintf('%s/repos/testasymmetry',Home));
addpath(sprintf('%s/repos/testasymmetry/functions',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley',Home));
addpath(sprintf('%s/repos/testasymmetry/born',Home));
addpath(sprintf('%s/repos/testasymmetry/two-region',Home));

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet2" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants3 ProblemSet3 saveMemory writeLogfile logfileName
logfileName = '2D1S.out';
epsOut = 78.36;

curdir = pwd;

ParamWatInfo = load('OptWater_thermo');
x = ParamWatInfo.xvec(3,:);
x_PB = [0 0 0 0 0 ParamWatInfo.xvec(3,6:7)];

fid = fopen('two-region.csv','r');
Data = textscan(fid,'%s %s %f %f %f','delimiter',',');
fclose(fid);
mol_list1 = Data{1};
mol_list2 = Data{2};
all_surfAreas1 = Data{3};%%%FIX
dG_list = Data{5};
surfArea_list1 = all_surfAreas1;
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
UsefulConstants3 = struct('epsIn1',epsIn1,'epsIn2',epsIn2,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);

% here we define the actual params for the NLBC test
% assuming that meshes are generated using Python3 command
% twoRegMeshGen(d0,spacing,df) (after running meshGen2D1S.py)
% this test is for d0 = 2.6 A, spacing = 0.3A and df = 6.5A
i = 1; %problem index
d0 = 3;
spacing = 0.5;
df = 6.5;

for j=d0:spacing:df
    
    dir=sprintf('%s/mesh-2D1S/%4.2f',curdir,j);
    chdir(dir);
    pqrData1 = loadPqr('mol1.pqr');
    pqrData2 = loadPqr('mol2.pqr');
    pqrAll1{i} = pqrData1;
    pqrAll2{i} = pqrData2;
    mol_list1{i} = 'Na-';
    mol_list2{i} = 'Na+';
    srfFile1{i} = sprintf('%s/mol1.srf',dir);
    srfFile2{i} = sprintf('%s/mol2.srf',dir);
    srfFile3{i} = sprintf('%s/mol12.srf',dir);
    chargeDist1{i} = pqrData1.q;
    chargeDist2{i} = pqrData2.q;
    referenceData{i} = -180; %%%%%%%%%%% fix later
    surfArea12{i} = surfArea_list1(1);
    chdir(curdir);
    addProblemSA3(mol_list1{i},mol_list2{i},pqrAll1{i},pqrAll2{i},srfFile1{i},srfFile2{i},srfFile3{i},...
               chargeDist1{i},chargeDist2{i},referenceData{i},surfArea12{i});
    i = i + 1;
end
chdir(curdir);
distance=(d0:spacing:df).';
[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA3(x);
SLICFileName = sprintf('Run2D1S_SLIC_Picard20_Na(-1)_Na(+1)');
save(SLICFileName,'distance','mol_list1','mol_list2','errfinal','calcE','refE','es','np');
[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA3(x_PB);
PBFileName = sprintf('Run2D1S_PB_Picard20_Na(-1)_Na(+1)');
save(PBFileName,'distance','mol_list1','mol_list2','errfinal','calcE','refE','es','np');

