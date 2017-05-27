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
epsOut = 1.8819; % from MNSol

mytemp=Tbase;
KelvinOffset = 273.15;
conv_factor = 332.112;
staticpotential = 0.0; % this only affects charged molecules;
kappa = 0.0;  % should be zero, meaning non-ionic solutions!


% the staticpotential below should not be used any more, please check
UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
			 kappa,'conv_factor',conv_factor,...
			 'staticpotential',staticpotential);
     
fid = fopen('mnsol/hexane.csv','r'); 
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

testset  = {'33_dimethylbutan_2_one','ethanol','ethyl_acetate','ethylamine','heptan_1_ol','hexan_2_one','methanol','methyl_pentanoate','n_propyl_acetate','propan_1_ol'};


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
  referenceData{i} = dG_list(index);
  surfArea{i} = surfArea_list(index);
  chdir(curdir);
  addProblemSA(testset{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0.5 -60 -0.5   -0.5*tanh(- -0.5)     0 -0.03 0];
lb = [0.4 -200 -100 -1  -0.1  -0.1  -4];
ub = [0.6 +200 +100 +1  +0.1  +0.1  +4];

options = optimoptions('lsqnonlin','MaxIter',6);
options = optimoptions(options,'Display', 'iter');

y = @(x)ObjectiveFromBEMSA(x);
[x,resnorm,residual,exitflag,output,] = lsqnonlin(y,x0,lb,ub,options);
[err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
[err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);

save('OptHexane','x','ref','calc','es','np','x0','calc0','es0','np0');
