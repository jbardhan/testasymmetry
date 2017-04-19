clear all
close all
clc

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
paramboundflag=1;   % if paramboundflag=0 there is no bound for parameters in the optimization process 
                    % if paramboundflag=1 there is abound
                    
calcflag=0;     % if calcflag=1 the code actually calculate the /delta G 's using BEM if it is zero, it means 
                % delta G's has been calculated before and all we need is
                % to load the data. 
                    
temp_min=5;     % lower bound of the temperature interval 
temp_max=45;    % upper bound in the temperature interval
tempdiv=5;      % number of divisions in the temperature interval                     
                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));


if calcflag==1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load data from parameterization


    %%% define the new tempereture vector 
    new_temp=linspace(temp_min,temp_max,9);

    %%% Finding the values of parameters from the quadratic fit at the new
    %%% temperatures
    load('param_func')
    param_num=length(paramfunc(1,:));

    for i=1:length(new_temp)
        for j=1:param_num
            x_new(i,j)=paramfunc(j).func(new_temp(i));
        end
    end

    for j=1:length(new_temp)
        clear global
        loadConstants

        convertKJtoKcal = 1/joulesPerCalorie;
        global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
        saveMemory = 0;
        writeLogfile = 0;
        logfileName = 'junklogfile';
        temp=new_temp(j);
        epsIn  =  1;
        Tbase = 300; 
        %epsOut = 78.36; % from MNSol
        epsOut = (-1.410e-6)*temp^3+(9.398e-4)*temp^2-0.40008*temp+87.740;
        mytemp=Tbase;
        KelvinOffset = 273.15;
        conv_factor = 332.112;
        staticpotential = 0.0; % this only affects charged molecules;
        kappa = 0.0;  % should be zero, meaning non-ionic solutions!


        % the staticpotential below should not be used any more, please check
        UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
                     kappa,'conv_factor',conv_factor,...
                     'staticpotential',staticpotential);

        if ionflag==1
        %    testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane','Li','Na','K','Rb','Cs','F','Cl','Br','I'};
            testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane','Li','Na','K','Rb','Cs','Cl','Br','I'};  % test set without florine

        elseif ionflag==0
            testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane'};
        end         

        fid = fopen('~/Research/testasymmetry/mobley/mnsol/mobley_sa.csv','r');
        Data = textscan(fid,'%s %f','delimiter',',');
        fclose(fid);
        all_solutes = Data{1};
        all_surfAreas = Data{2};
        [m, index] = ismember(testset,all_solutes);
        surfArea_list = all_surfAreas(index);
        t_ref_aca=24.85; %reference tempereture for amino acid analogouses
        t_ref_ion=25;  %reference tempereture for ions

        if ionflag==1
            dG_list_ref_at_298=[8.1,-40.5,-5.2,9,9.5,-6.2,-3.2,-21.2,-20.4,-24.6,-25.6,8.3]'./joulesPerCalorie; %Hess in kcal/mol
            H_list_ref_at_298=[-8.3,-67.0,-23.9,-17.1,-17.1,-34.6,-25.3,-43.0,-45,-58.8,-57.4,-13.7]'./joulesPerCalorie;  %Hess . in kcal/mol
            dS_list_ref_at_298=(H_list_ref_at_298-dG_list_ref_at_298(1:12))/298;  % in kcal/mol/K

          %  dG_list_ref_ion_at_298=[-529;-424;-352;-329;-306;-429;-304;-278;-243]./joulesPerCalorie;         %Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference
          %  dS_list_ref_ion_at_298=[-0.164;-0.133;-0.096;-0.087;-0.081;-0.115;-0.053;-0.037;-0.014]./joulesPerCalorie;   %Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference

            dG_list_ref_ion_at_298_15=[-529;-424;-352;-329;-306;-304;-278;-243]./joulesPerCalorie;         % with out florine Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference
            dS_list_ref_ion_at_298_15=[-0.164;-0.133;-0.096;-0.087;-0.081;-0.053;-0.037;-0.014]./joulesPerCalorie;   % with out florine Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference

            dG_list_aca=dG_list_ref_at_298-dS_list_ref_at_298*(temp-t_ref_aca);
            dG_list_ion=dG_list_ref_ion_at_298_15-dS_list_ref_ion_at_298_15*(temp-t_ref_ion);

            dG_list=[dG_list_aca;dG_list_ion];
            dS_list=[dS_list_ref_at_298;dS_list_ref_ion_at_298_15]; % reference list of Entropy, note that teh entropy for amino acids analogouses are at 298K but the entropy for ions are at 298.15
                               
    %         
    %         dG_list_ref_at_298=[dG_list_ref_at_298;dG_list_ref_ion_at_298];   % in kcal/mol
    %         dS_list_ref_at_298=[dS_list_ref_at_298;dS_list_ref_ion_at_298];   % in kcal/mol/K    
    %         dG_list=dG_list_ref_at_298-dS_list_ref_at_298*(TEMP(kk)-t_ref);

        elseif ionflag==0
            dG_list_ref_at_298=[8.1,-40.5,-5.2,9,9.5,-6.2,-3.2,-21.2,-20.4,-24.6,-25.6,8.3]'./joulesPerCalorie;%Hess .  in kcal/mol
            H_list_ref_at_298=[-8.3,-67.0,-23.9,-17.1,-17.1,-34.6,-25.3,-43.0,-45,-58.8,-57.4,-13.7]'./joulesPerCalorie;  % Hess in kcal/mol
            dS_list_ref_at_298=(H_list_ref_at_298-dG_list_ref_at_298)/298;  % in kcal/mol/K
            dG_list=dG_list_ref_at_298-dS_list_ref_at_298*(TEMP(kk)-t_ref_aca);  % in kcal/mol
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
        [errfinal(j,:),calcE(j,:),refE(j,:),es(j,:),np(j, :)]=ObjectiveFromBEMSA(x_new(j,:));
    end

    save('RunWater_param','errfinal','calcE','refE','es','np','new_temp','x_new','testset','dS_list');
    
end

RunWater_param=load('RunWater_param.mat');

dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures
paramfunc=struct();  % structure that has the information for the quadratic fit of the parameters as a function of tempereture

loadConstants
convertKJtoKcal = 1/joulesPerCalorie;

x = RunWater_param.x_new;   % parameters at different temperetures
calcE = RunWater_param.calcE; %calculated values for \Delta G at different temperatures
refE = RunWater_param.refE; %calculated values for \Delta G at different temperatures
refS = RunWater_param.dS_list; %calculated values for \Delta G at different temperatures
TEMP=RunWater_param.new_temp;        % create the temperature vector 
testset=RunWater_param.testset;

for i=1:length(testset)
    dGfunc(i).name=testset(i);
    dGfunc(i).func=fit(TEMP',calcE(:,i),'poly1'); % the linear function that fits to the calculated dG  at differet temperatures
end

for i=1:length(testset)
    p=[dGfunc(i).func.p1,dGfunc(i).func.p2];
    pder=polyder(p);  % derivative of the linear function dG
    if i<13
        dsvec(i)=-polyval(pder,24.85)*1000;    % Evaluationg entropy, dS at 298K = 24.85C in cal/mol/K
    else
        dsvec(i)=-polyval(pder,25)*1000;    % Evaluationg entropy, dS at 298K = 24.85C in cal/mol/K
%     figure(i);
%     plot(dGfunc(i).func,TEMP',calcvec(:,i),'o')   
    end
end

Reference_dG =refE
Calculated_dG=calcE

dG_err=abs(Reference_dG-Calculated_dG)

Reference_dS =refS'*1000
Calculated_ds=dsvec

dS_err=abs(Reference_dS-Calculated_ds)


    
    


    