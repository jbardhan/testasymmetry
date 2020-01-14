
clear all
clc

% Path information
Home = getenv('HOME');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%                                Set these values before running the code                      %%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

repo_path=sprintf('%s/repos',Home);
dropbox_path=sprintf('%s/Dropbox',Home);

ionflag=0;          % if ionflag=0, ions data are not included in the testset 
                    % if ionflag=1, ions data are included in the testset
paramboundflag=1;   % if paramboundflag=0 there is no bound for parameters in the optimization process 
                    % if paramboundflag=1 there is abound
                   
temp_min=24.85;     % lower bound of the temperature interval 
temp_max=84.85;    % upper bound in the temperature interval
tempdiv=5;      % number of divisions in the temperature interval                     
                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Note: we calculate everything at 298 K which is equal to 24.85C, We
%%%% will calculate \Delta G of ions which is available at 25C at 24.85 C
%%%% to use them in our objective function. But in the postprocessing
%%%% script (whaterthermo.m) we will compare 
                    
addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley/reference-data',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));                    

TEMP=linspace(temp_min,temp_max,tempdiv); % create the temperature vector

% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.

for j=1:tempdiv
    clear global
    loadConstants
    convertKJtoKcal = 1/joulesPerCalorie;
    global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
    saveMemory = 0;
    writeLogfile = 0;
    logfileName = 'junklogfile';

    epsIn  =  1;
    KelvinOffset = 273.15;
    conv_factor = 332.112;
    staticpotential = 0.0; % this only affects charged molecules;
    kappa = 0.0;  % should be zero, meaning non-ionic solutions!
    tempKelvin=TEMP(j)+KelvinOffset;
    epsOut = 12.85; 
    % See Choi13 (Fig 8) or Krossing06

    
    % the staticpotential below should not be used any more, please check
    UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
                 kappa,'conv_factor',conv_factor,...
                 'staticpotential',staticpotential);

    if ionflag==1
        fprintf('Error: paramOctanol does not treat ions yet.\n');
    elseif ionflag==0
        testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane',...
                    '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol',...
                    'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane'};
    end

    fid = fopen('~/repos/testasymmetry/mobley/mnsol/mobley_dG_AND_sa_and_vol_fixed.csv','r');
    Data = textscan(fid,'%s %f  %f  %f  %f  %f  %f  %f','delimiter',',');
    fclose(fid);
    all_solutes = Data{1};
    all_dG = Data{2};
    all_surfAreas = Data{3};
    [m, index] = ismember(testset,all_solutes);
    surfArea_list = all_surfAreas(index);
    t_ref_aca=24.85; %reference tempereture for amino acid analogues
    
    [dGs,error_code] = determinePaluchSolvationFreeEnergy('reference-data/mobley_paluch_omim_tf2n.csv',testset,tempKelvin);
    dG_list = dGs*kB*tempKelvin;
    %return;
    if error_code > 0
      fprintf(['Halting optimization because' ...
         ' determinePaluchSolvationFreeEnergy failed.\n']);
    end
    curdir=pwd;
    for i=1:length(testset)
      dir=sprintf('%s/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',dropbox_path,testset{i});
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
    x0 = [1.333	-77.589	-2.229	-1.538	0.078	-0.018	3.382]; 
    %from SLIC w/o ion for octanol 298
    
    if ionflag==0
        if paramboundflag==1
            lb = [0  -200 -100 -20  -0.1  -0.1  -4];
            ub = [+2 +200 +100 +20  +0.1  +0.1  +4];
        elseif paramboundflag==0
            lb = [-inf -inf -inf -inf  -inf  -inf  -inf];
            ub = [ inf  inf  inf  inf   inf   inf   inf];
        end
    elseif ionflag==1
        if paramboundflag==1
            lb = [0  -200 -100 -20  -20  -0.1  -4];
            ub = [+2 +200 +100 +20  +20  +0.1  +4];
        elseif paramboundflag==0
            lb = [-inf -inf -inf -inf  -inf  -inf  -inf];
            ub = [ inf  inf  inf  inf   inf   inf   inf];
        end
    end

    options = optimoptions('lsqnonlin','MaxIter',8);
    options = optimoptions(options,'Display', 'iter');

    y = @(x)ObjectiveFromBEMSA(x);
    [x,resnorm,residual,exitflag,output] = lsqnonlin(y,x0,lb,ub,options);
    [err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
    [err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);

    xvec(j,:)=x;
    refvec(j,:)=ref;
    calcvec(j,:)=calc;
    esvec(j,:)=es;
    npvec(j,:)=np;
    x0vec(j,:)=x0;
    calc0vec(j,:)=calc0;
    es0vec(j,:)=es0;
    np0vec(j,:)=np0;
    tempvec(j,:)=tempKelvin;
    epsvec(j,:)=epsOut;
end
if ionflag==0
    if paramboundflag==1
        save('OptOctanolPaluch_wo_ion','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','tempvec','epsvec');
    else
        save('OptOctanolPaluch_wo_ion_wo_bound','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','tempvec','epsvec');
    end
elseif ionflag==1
    if paramboundflag==1
        save('OptOctanolPaluch_w_ion','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','tempvec','epsvec');
    else
        save('OptOctanolPaluch_w_ion_wo_bound','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','tempvec','epsvec');
    end
end
