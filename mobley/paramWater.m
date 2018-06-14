
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


    testset  = {'Li','Na','K','Rb','Cs','Cl','Br','I'};  % test set without florine



    fid = fopen('~/repos/testasymmetry/mobley/mnsol/mobley_dG_AND_sa_and_vol_fixed_ions.csv','r');
    Data = textscan(fid,'%s %f  %f  %f  %f  %f  %f  %f','delimiter',',');
    fclose(fid);
    
    if ionflag==1
        
        all_solutes = Data{1};
        all_dG = Data{2};
        all_surfAreas = Data{3};
        dG_np = Data{6};% nonpolar energies
        dG_es = Data{7};% electrotatic energies
        
    elseif ionflag==0
        
        all_solutes = Data{1}(1:502);
        all_dG = Data{2}(1:502);
        all_surfAreas = Data{3}(1:502);
        dG_np = Data{6}(1:502);% nonpolar energies (no ion)
        dG_es = Data{7}(1:502);% electrotatic energies
        
    end

    [m, index] = ismember(testset,all_solutes);
    surfArea_list = all_surfAreas(index);
    dG_es_list = dG_es(index);
    t_ref_aca=24.85; %reference tempereture for amino acid analogues
    t_ref_ion=25;  %reference tempereture for ions
    
    dG_list_ref_at_298=[8.1,-40.5,-5.2,9,9.5,-6.2,-3.2,-21.2,-20.4,-24.6,-25.6,8.3]'./joulesPerCalorie; %Hess in kcal/mol
    H_list_ref_at_298=[-8.3,-67.0,-23.9,-17.1,-17.1,-34.6,-25.3,-43.0,-45,-58.8,-57.4,-13.7]'./joulesPerCalorie;  %Hess . in kcal/mol
    CP_list_ref_at_298=1e-3*[142,25,213,310,290,86,285,44,122,336,174,246]'./joulesPerCalorie;
    dS_list_ref_at_298=(H_list_ref_at_298-dG_list_ref_at_298(1:12))/298;  % in kcal/mol/K
    dG_list_aca=dG_list_ref_at_298-dS_list_ref_at_298*(TEMP(kk)-t_ref_aca)+CP_list_ref_at_298*((TEMP(kk)-t_ref_aca)-(TEMP(kk)+KelvinOffset)*log(((TEMP(kk)+KelvinOffset))/((t_ref_aca+KelvinOffset))));
    
    aca_num=length(dG_list_ref_at_298);
    ion_num=0;
    
    if ionflag==1
      
        dG_list_ref_ion_at_298_15=[-529;-424;-352;-329;-306;-304;-278;-243]./joulesPerCalorie;         % with out florine Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference
        dS_list_ref_ion_at_298_15=[-0.164;-0.133;-0.096;-0.087;-0.081;-0.053;-0.037;-0.014]./joulesPerCalorie;   % with out florine Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference
        CP_list_ref_ion_at_298_15=1e-3*[-23;-42;-72;-94;-108;-70;-74;-64]./joulesPerCalorie;
        dG_list_ion=dG_list_ref_ion_at_298_15-dS_list_ref_ion_at_298_15*(TEMP(kk)-t_ref_ion)+CP_list_ref_ion_at_298_15*((TEMP(kk)-t_ref_ion)-(TEMP(kk)+KelvinOffset)*log(((TEMP(kk)+KelvinOffset))/((t_ref_ion+KelvinOffset))));
    
        
        dG_list=[dG_list_aca;dG_list_ion];
        dS_list=[dS_list_ref_at_298;dS_list_ref_ion_at_298_15];
        CP_list=[CP_list_ref_at_298;CP_list_ref_ion_at_298_15];
        
        ion_num=length(dG_list_ref_ion_at_298_15);
        
    elseif ionflag==0
        
        dG_list=dG_list_aca;
        dS_list=dS_list_ref_at_298;
        CP_list=CP_list_ref_at_298;
        
    end
    
    NpInfo=load('Np_opt');
    xnp = NpInfo.x_np;
    dG_list_ion = dG_list_ion - surfArea_list.*xnp(kk,1) - xnp(kk,2);
    testset  = {'Li','Na','K','Rb','Cs','Cl','Br','I'};
    
    
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
      referenceData{i} = dG_list_ion(i);
      surfArea{i} = surfArea_list(i);
      
      chdir(curdir);
      addProblemSA(testset{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
    end


    % The following script is specialized to this example.  We'll
    % handle generating others.  Not complicated, but it's not self-explanatory.


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % alpha beta gamma mu phi_stat np_a np_b
    x0 = [0.5 -60 -0.5   -0.5*tanh(- -0.5)  0 0 0];
    
    lb = [-2 -200 -100 -20  -20  -0.00001  -0.0001];
    ub = [+2 +200 +100 +20  +20  +0.00001  +0.0001];

    options = optimoptions('lsqnonlin','MaxIter',8);
    options = optimoptions(options,'Display', 'iter');

    y = @(x)ObjectiveFromBEMSA(x);
    [x,resnorm,residual,exitflag,output] = lsqnonlin(y,x0,lb,ub,options);
    x = [x(1) x(2) x(3) x(4) x(5) xnp(kk,1) x(kk,2)];
    [err,calc,ref,es,np]=ObjectiveFromBEMSA(x);

    x0 = [0.5 -60 -0.5   -0.5*tanh(- -0.5)  0 0 0];
    [err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);

    xvec(kk,:)=x;refvec(kk,:)=ref;calcvec(kk,:)=calc;
    esvec(kk,:)=es;npvec(kk,:)=np;x0vec(kk,:)=x0;
    calc0vec(kk,:)=calc0;es0vec(kk,:)=es0;np0vec(kk,:)=np0;
    tempvec(kk,:)=temp;
end
if ionflag==0
    save('OptWater_thermo_wo_ion','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','testset','dS_list','CP_list','tempvec','ionflag','aca_num','ion_num','t_ref_aca','t_ref_ion');
elseif ionflag==1
    save('OptWater_thermo_ion_only','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','testset','dS_list','CP_list','tempvec','ionflag','aca_num','ion_num','t_ref_aca','t_ref_ion');
end
