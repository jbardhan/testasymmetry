% Path information
Home = getenv('HOME');
addpath(sprintf('%s/repos/pointbem',Home));
addpath(sprintf('%s/repos/panelbem',Home));
addpath(sprintf('%s/repos/testasymmetry',Home));
addpath(sprintf('%s/repos/testasymmetry/functions',Home));
addpath(sprintf('%s/repos/testasymmetry/mobley',Home));
addpath(sprintf('%s/repos/testasymmetry/born',Home));
addpath(sprintf('%s/repos/testasymmetry/membrane',Home));

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
logfileName = 'membrane.out';
epsOut = 78.36;

curdir = pwd;
chdir('../mobley');
ParamWatInfo = load('OptWater_thermo_Picard_20');
x = ParamWatInfo.xvec(3,:);
x_PB = [0 0 0 0 0 ParamWatInfo.xvec(3,6:7)];

chdir(curdir);
lambdaFEP = linspace (0,1,11);

%fid = fopen('membrane_exp.csv','r'); 
%Data = textscan(fid,'%s %f %f','delimiter',',');
%fclose(fid);
%dG_list = Data{2};
% No reference data for membrane yet!
dG_list = zeros(length(lambdaFEP),1);

%fid = fopen('SLIC_Membrane.csv','r');
%Data = textscan(fid,'%s %f','delimiter',',');
%fclose(fid);
%all_solutes = Data{1};
%all_surfAreas = Data{2};
%[m, index] = ismember(mol_list,all_solutes);
%surfArea_list = all_surfAreas(index);
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


lambdaFEP = linspace (0,10,11);
for k=1:length(lambdaFEP)
    mol_list{k}=sprintf('test_lambda_%3.1f',lambdaFEP(k));
end

% Radius of the biggest cylinder to run SLIC (starting from R=1A with 1A spacing) 
%lastCylRadiusToCalc = 10;

%Cylinder height (Angstrom)
h_Cylinder = 50;
r0 = 15;
spacing = 3;
rf = 20;
radii = linspace(r0,rf,spacing-1);
Radii = 20.0;%[5.0,7.5,10.0,15.0,20.0,30.0,40.0];
% Cylinder (membrane) radii (Angstrom)
%cylinderRadius = linspace(1,lastCylRadiusToCalc,lastCylRadiusToCalc);

for j=1%:length(Radii)
    dir=sprintf('%s/mesh-membrane/%3.1f',curdir,Radii(j));
    chdir(dir);
    for i=1%:11
        pqrFile=sprintf('%s.pqr',mol_list{i});
        pqrData = loadPqr(pqrFile);
        pqrAll{i} = pqrData;
        srfFile{i} = sprintf('%s/test_2.srf',dir);
        chargeDist{i} = pqrData.q;
        referenceData{i} = 0;
        surfArea_list = pi*(2*Radii(j)^2 + 2 * Radii(j) * h_Cylinder);
        addProblemSA(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea_list);
    end
    distance=radii.';
    [errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA(x);
    SLIC_membraneFileName = sprintf('RunMembrane_SLIC_po4_%3.1f.mat',Radii(j));
    save(SLIC_membraneFileName,'distance','mol_list','errfinal','calcE','refE','es','np');
    [errfinal,calcE,refE,es,np]=ObjectiveFromBEMSA(x_PB);
    PB_membraneFileName = sprintf('RunMembrane_PB_po4_%3.1f.mat',Radii(j));
    save(PB_membraneFileName,'distance','mol_list','errfinal','calcE','refE','es','np');
    
    chdir(curdir)
    ProblemSet = [];
end







