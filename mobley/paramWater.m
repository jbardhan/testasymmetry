
clear all


% Path information
Home = getenv('HOME');
addpath(sprintf('%s/Research/pointbem',Home));
addpath(sprintf('%s/Research/panelbem',Home));
addpath(sprintf('%s/Research/testasymmetry',Home));
addpath(sprintf('%s/Research/testasymmetry/functions',Home));
addpath(sprintf('%s/Research/testasymmetry/mobley',Home));
addpath(sprintf('%s/Research/testasymmetry/born',Home));



tempdiv=5;
TEMP=linspace(0,100,tempdiv);

% tempdiv=1;
% TEMP=24.85;

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.

for kk=1:tempdiv
    clear global
    loadConstants
    convertKJtoKcal = 1/joulesPerCalorie;
    global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
    saveMemory = 0;
    writeLogfile = 0;
    logfileName = 'junklogfile';
    temp=TEMP(kk);
    epsIn  =  1;
    Tbase = 300; 
    %epsOut = 78.36; % from MNSol
    epsOut = (-1.410e-6)*TEMP(kk)^3+(9.398e-4)*TEMP(kk)^2-0.40008*TEMP(kk)+87.740;
    mytemp=Tbase;
    KelvinOffset = 273.15;
    conv_factor = 332.112;
    staticpotential = 0.0; % this only affects charged molecules;
    kappa = 0.0;  % should be zero, meaning non-ionic solutions!


    % the staticpotential below should not be used any more, please check
    UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
                 kappa,'conv_factor',conv_factor,...
                 'staticpotential',staticpotential);

%     fid = fopen('mnsol/water.csv','r'); 
%     Data = textscan(fid,'%s %f %f','delimiter',',');
%     fclose(fid);
%     mol_list = Data{1};
%     dG_list = Data{2};
% 
%     fid = fopen('mnsol/mobley_sa.csv','r');
%     Data = textscan(fid,'%s %f','delimiter',',');
%     fclose(fid);
%     all_solutes = Data{1};
%     all_surfAreas = Data{2};
%     [m, index] = ismember(mol_list,all_solutes);
%     surfArea_list = all_surfAreas(index);

     testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane'};

%testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'toluene', 'methanol', 'ethanol', 'p_cresol', 'propane'};


    fid = fopen('~/Research/testasymmetry/mobley/mnsol/mobley_sa.csv','r');
    Data = textscan(fid,'%s %f','delimiter',',');
    fclose(fid);
    all_solutes = Data{1};
    all_surfAreas = Data{2};
    [m, index] = ismember(testset,all_solutes);
    surfArea_list = all_surfAreas(index);
    dG_list_ref_at_298=[8.1,-40.5,-5.2,9,9.5,-6.2,-3.2,-21.2,-20.4,-24.6,-25.6,8.3]'./ 4.184;%Hess
    
   % dG_list_ref_at_298=[2,-9.71,-1.24,2.08,2.32,-0.89,-5.11,-5.01,-6.14,1.96]';%MNSol

    
    
    H_list_ref_at_298=[-8.3,-67.0,-23.9,-17.1,-17.1,-34.6,-25.3,-43.0,-45,-58.8,-57.4,-13.7]'./ 4.184;
    dS_list_ref_at_298=(H_list_ref_at_298-dG_list_ref_at_298)/298;
    t_ref=24.85;
    
    dG_list=dG_list_ref_at_298-dS_list_ref_at_298*(TEMP(kk)-t_ref);
    
    
   

    % all octanol available side chain analogues 
    %testset = {'2_methylpropane', 'acetic_acid', 'ethanol', 'methane', 'methanol',...
    % 'n_butane', 'n_butylamine', 'p_cresol', 'propane', 'propanoic_acid','toluene'};

    % complete list of side chain analogues. not available for all solvents
    %testset = {'1_methyl_imidazole','2_methylpropane', ...
    %	   '3_methyl_1h_indole','acetic_acid','ethanamide', ...
    %	   'ethanol','methane','methanethiol','methanol', ...
    %	   'methyl_ethyl_sulfide','n_butane','n_butylamine', ...
    %	   'p_cresol','propane','propanoic_acid','toluene'};
    curdir=pwd;
    for i=1:length(testset)
      dir=sprintf('%s/Dropbox-NEU/Dropbox/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',getenv('HOME'),testset{i});
      chdir(dir);
      pqrData = loadPqr('test.pqr');
      pqrAll{i} = pqrData;
      srfFile{i} = sprintf('%s/test_2.srf',dir);
      chargeDist{i} = pqrData.q;%chargeDistribution;
      foo = strcmp(all_solutes,testset{i});
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
    % alpha beta gamma mu phi_stat np_a np_b
    x0 = [0.5 -60 -0.5   -0.5*tanh(- -0.5)     0 -0.03 1.6];
    lb = [-2 -200 -100 -20  -0.1  -0.1  -2];
    ub = [+2 +200 +100 +20  +0.1  +0.1  +2];

    options = optimoptions('lsqnonlin','MaxIter',8);
    options = optimoptions(options,'Display', 'iter');

    y = @(x)ObjectiveFromBEMSA(x);
    [x,resnorm,residual,exitflag,output] = lsqnonlin(y,x0,lb,ub,options);
    [err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
    [err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);

    xvec(kk,:)=x;
    refvec(kk,:)=ref;
    calcvec(kk,:)=calc;
    esvec(kk,:)=es;
    npvec(kk,:)=np;
    x0vec(kk,:)=x0;
    calc0vec(kk,:)=calc0;
    es0vec(kk,:)=es0;
    np0vec(kk,:)=np0;
    tempvec(kk,:)=temp;
end

save('OptWater','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','tempvec');