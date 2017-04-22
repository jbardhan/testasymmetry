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
                    
temp_min=278;     % lower bound of the temperature interval 
temp_max=318;    % upper bound in the temperature interval
tempdiv=5;      % number of divisions in the temperature interval                     
                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));


               

if ionflag==1
  %  testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane','Li','Na','K','Rb','Cs','F','Cl','Br','I'};
    testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane','Li','Na','K','Rb','Cs','Cl','Br','I'}; % test set with out florine

elseif ionflag==0
    testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane'};
end


dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures
paramfunc=struct();  % structure that has the information for the quadratic fit of the parameters as a function of tempereture

loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
KelvinOffset = 273.15;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data from parameterization

if ionflag==1
    if paramboundflag==0
        ParamWatInfo = load('OptWater_w_ion_wo_bound');
    elseif paramboundflag==1
        %ParamWatInfo = load('OptWater_w_ion');
        ParamWatInfo = load('OptWater_w_ion_wo_florine');
    end
elseif ionflag==0
    if paramboundflag==0
        ParamWatInfo = load('OptWater_wo_ion_wo_bound');
    elseif paramboundflag==1
        ParamWatInfo = load('OptWater_wo_ion');
    end
end
      
x = ParamWatInfo.xvec;   % parameters at different temperetures
calcvec = ParamWatInfo.calcvec; %calculated values for \Delta G at different temperatures
TEMP=linspace(temp_min,temp_max,tempdiv);        % create the temperature vector  
TEMP=TEMP-KelvinOffset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_num=length(x(1,:));  % number of parameters
param_list={'alpha','beta','gamma','mu','phi_stat','np_a','np_b'};   % list of parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fitting the calculated \Delta G's in optimized tempereture as a linear function of temperature

for i=1:length(testset)
    dGfunc(i).name=testset(i);
    dGfunc(i).func=fit(TEMP',calcvec(:,i),'poly1'); % the linear function that fits to the calculated dG  at differet temperatures
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fitting the parameters as a quadratic function of temperature

for i=1:param_num
    paramfunc(i).name=param_list(i);
    paramfunc(i).func=fit(TEMP',x(:,i),'poly2'); % the quadratic function that fits to the calculated dG  at differet temperatures
%     figure
%     plot(paramfunc(i).func,TEMP',x(:,i),'o')
end

save ('param_func','paramfunc');



    
    



