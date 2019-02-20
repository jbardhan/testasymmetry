% Path information
Home = getenv('HOME');
addpath(sprintf('%s/repos/pointbem',Home));
addpath(sprintf('%s/repos/panelbem',Home));
addpath(sprintf('%s/repos/testasymmetry',Home));
addpath(sprintf('%s/repos/testasymmetry/functions',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley',Home));
addpath(sprintf('%s/repos/testasymmetry/born',Home));
addpath(sprintf('%s/repos/testasymmetry/two-region',Home));

%This is the original SLIC code (one dielectric surface and its
%corresponding stern layer in solvent) which is used to validate the other
%extensions of SLIC (two-region, membrane etc.)
% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet2" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
logfileName = '1D1S.out';
epsOut = 78.36;

curdir = pwd;

ParamWatInfo = load('OptWater_thermo');
x = ParamWatInfo.xvec(3,:);
x(5) =0 ;
x_PB = [0 0 0 0 0 ParamWatInfo.xvec(3,6:7)];

fid = fopen('two-region.csv','r');
Data = textscan(fid,'%s %s %f %f %f','delimiter',',');
fclose(fid);
mol_list1 = Data{1};
%mol_list2 = Data{2};
all_surfAreas1 = Data{3};
%all_surfAreas2 = Data{4};
dG_list = Data{5};
surfArea_list1 = all_surfAreas1;
%surfArea_list2 = all_surfAreas2;
saveMemory = 1;
writeLogfile = 1;
epsIn  =  1;
%epsIn2  =  1;
Tbase = 300; mytemp=Tbase;
staticpotential = 0.0;
KelvinOffset = 273.15;
conv_factor = 332.112;
%staticpotential = 0.0; % this only affects charged molecules;

kappa = 0.0;  % should be zero, meaning non-ionic solutions!
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);

% here we define the actual params for the NLBC test
% assuming that meshes are generated using Python3 command
% twoRegMeshGen(d0,spacing,df) (after running twoRegMeshGen.py)
% this test is for d0 = 3.0 A, spacing = 0.5A and df = 12.0A
i = 1; %problem index
%d0 = 10.0;
%spacing = 10.0;
%df = 30.0;

for j=1
    
    dir=sprintf('%s/%s',curdir,'Na');
    chdir(dir);
    pqrData1 = loadPqr('test.pqr');
    %pqrData2 = loadPqr('mol2.pqr');
    pqrAll1{i} = pqrData1;
    %pqrAll2{i} = pqrData2;
    mol_list1{i} = 'Na+';
    %mol_list2{i} = 'Na-';
    srfFile1{i} = sprintf('%s/test_2.srf',dir);
    %srfFile2{i} = sprintf('%s/mol2.srf',dir);
    chargeDist1{i} = pqrData1.q;
    %chargeDist2{i} = pqrData2.q;
    referenceData{i} = ParamWatInfo.refvec(3,14); %Na
    surfArea1{i} = surfArea_list1(1);
    %surfArea2{i} = surfArea_list2(1);
    chdir(curdir);
    addProblemSA(mol_list1{i},pqrAll1{i},srfFile1{i},...
               chargeDist1{i},referenceData{i},surfArea1{i});
    i = i + 1;
end
chdir(curdir);
%distance=(d0:spacing:df).';
[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA(x);
SLICFileName = sprintf('Run1D1S_SLIC_Picard100_Na(-1)_test');
save(SLICFileName,'mol_list1','errfinal','calcE','refE','es','np');
%[errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA(x_PB);
%PBFileName = sprintf('Run1D1S_PB_Na(-1)');
%save(PBFileName,'mol_list1','errfinal','calcE','refE','es','np');

