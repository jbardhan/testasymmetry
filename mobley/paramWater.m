
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

repo_path=sprintf('%s/Research',Home);
dropbox_path=sprintf('%s/Dropbox-NEU/Dropbox',Home);


ionflag=1;          % if ionflag=0, ions data are not included in the testset 
                    % if ionflag=1, ions data are included in the testset

temp_min=4.85;     % lower bound of the temperature interval 
temp_max=44.85;    % upper bound in the temperature interval
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
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));                    
 


TEMP=linspace(temp_min,temp_max,tempdiv);        % create the temperature vector


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

    epsIn  =  1;
    Tbase = 300; 
    %epsOut = 78.36; % from MNSol
    
    mytemp=Tbase;
    KelvinOffset = 273.15;
    conv_factor = 332.112;
    staticpotential = 0.0; % this only affects charged molecules;
    kappa = 0.0;  % should be zero, meaning non-ionic solutions!
    temp=TEMP(kk);
    epsOut = (-1.410e-6)*temp^3+(9.398e-4)*temp^2-0.40008*temp+87.740;  % function for dielectric constant
                                             % of water as a function of
                                             % temperature in C


    % the staticpotential below should not be used any more, please check
    UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
                 kappa,'conv_factor',conv_factor,...
                 'staticpotential',staticpotential);


    if ionflag==1
        testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane','Li','Na','K','Rb','Cs','Cl','Br','I'};  % test set without florine

    elseif ionflag==0
        testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane'};
    end


    fid = fopen('~/Research/testasymmetry/mobley/mnsol/mobley_dG_AND_sa_and_vol.csv','r');
    Data = textscan(fid,'%s %f  %f  %f  %f  %f  %f  %f','delimiter',',');
    fclose(fid);
    all_solutes = Data{1};
    all_dG = Data{2};
    all_surfAreas = Data{3};
    [m, index] = ismember(testset,all_solutes);
    surfArea_list = all_surfAreas(index);
    t_ref_aca=24.85; %reference tempereture for amino acid analogues
    t_ref_ion=25;  %reference tempereture for ions
    
    dG_list_ref_at_298=[8.1,-40.5,-5.2,9,9.5,-6.2,-3.2,-21.2,-20.4,-24.6,-25.6,8.3]'./joulesPerCalorie; %Hess in kcal/mol
    H_list_ref_at_298=[-8.3,-67.0,-23.9,-17.1,-17.1,-34.6,-25.3,-43.0,-45,-58.8,-57.4,-13.7]'./joulesPerCalorie;  %Hess . in kcal/mol
    dS_list_ref_at_298=(H_list_ref_at_298-dG_list_ref_at_298(1:12))/298;  % in kcal/mol/K
    dG_list_aca=dG_list_ref_at_298-dS_list_ref_at_298*(TEMP(kk)-t_ref_aca);
    
    aca_num=length(dG_list_ref_at_298);
    ion_num=0;
    
    if ionflag==1
      
        dG_list_ref_ion_at_298_15=[-529;-424;-352;-329;-306;-304;-278;-243]./joulesPerCalorie;         % with out florine Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference
        dS_list_ref_ion_at_298_15=[-0.164;-0.133;-0.096;-0.087;-0.081;-0.053;-0.037;-0.014]./joulesPerCalorie;   % with out florine Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference
        dG_list_ion=dG_list_ref_ion_at_298_15-dS_list_ref_ion_at_298_15*(TEMP(kk)-t_ref_ion);
        
        dG_list=[dG_list_aca;dG_list_ion];
        dS_list=[dS_list_ref_at_298;dS_list_ref_ion_at_298_15];
        
        ion_num=length(dG_list_ref_ion_at_298_15);
        
    elseif ionflag==0
        
        dG_list=dG_list_aca;
        dS_list=dS_list_ref_at_298;
        
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
    x0 = [0.5 -60 -0.5   -0.5*tanh(- -0.5)  0 -0.03 1.6];
    if ionflag==0
        lb = [-2 -200 -100 -20  -0.1  -0.1  -2];
        ub = [+2 +200 +100 +20  +0.1  +0.1  +2];
    elseif ionflag==1
            lb = [-2 -200 -100 -20  -20  -0.1  -2];
            ub = [+2 +200 +100 +20  +20  +0.1  +2];
    end

    options = optimoptions('lsqnonlin','MaxIter',8);
    options = optimoptions(options,'Display', 'iter');

    y = @(x)ObjectiveFromBEMSA(x);
    [x,resnorm,residual,exitflag,output] = lsqnonlin(y,x0,lb,ub,options);
    [err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
    [err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);

    xvec(kk,:)=x;refvec(kk,:)=ref;calcvec(kk,:)=calc;
    esvec(kk,:)=es;npvec(kk,:)=np;x0vec(kk,:)=x0;
    calc0vec(kk,:)=calc0;es0vec(kk,:)=es0;np0vec(kk,:)=np0;
    tempvec(kk,:)=temp;
end
if ionflag==0
    save('OptWater_thermo_wo_ion','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','testset','dS_list','tempvec','ionflag','aca_num','ion_num','t_ref_aca','t_ref_ion');
elseif ionflag==1
    save('OptWater_thermo','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','testset','dS_list','tempvec','ionflag','aca_num','ion_num','t_ref_aca','t_ref_ion');
end
