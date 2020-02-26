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

repo_path=sprintf('%s/repos',Home);
dropbox_path=sprintf('%s/Dropbox',Home);
                   
temp_min=24.85;     % lower bound of the temperature interval 
temp_max=84.85;    % upper bound in the temperature interval
tempdiv=5;      % number of divisions in the temperature interval                                                   11      
                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));


ParamOctanolInfo=load('OptOctanolPaluch_wo_ion');
x = ParamOctanolInfo.xvec;
dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures
loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
calcE = ParamOctanolInfo.calcvec; %calculated values for \Delta G at different temperatures
refE = ParamOctanolInfo.refvec; %calculated values for \Delta G at different temperatures
es=ParamOctanolInfo.esvec;
np=ParamOctanolInfo.npvec;
mol_list=ParamOctanolInfo.testset;
TEMP_K=ParamOctanolInfo.tempvec';     
[m,index]=ismember(328,TEMP_K);
%mol_list={'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide',...
%            'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane'}';
testset=mol_list;
for i=1:length(mol_list)
    f = @(R) (R(1)-R(2)*(TEMP_K-328)+R(3)*((TEMP_K-328)-TEMP_K.*log(TEMP_K./328)))-calcE(:,i);
    R0=[refE(index,i),1,1];
    options=optimoptions('lsqnonlin','StepTolerance',1e-6);
    options=optimoptions(options,'OptimalityTolerance',1e-6);
    options=optimoptions(options,'FunctionTolerance',1e-6);
    [R,resnorm,residual,exitflag,output]=lsqnonlin(f,R0,[],[],options);
    dGfunc(i).name=mol_list(i); 
    dGfunc(i).dg=R(1);
    dGfunc(i).ds=R(2);
    dGfunc(i).cp=R(3); 
    dsvec(i)=dGfunc(i).ds*1000;
    cpvec(i)=dGfunc(i).cp*1000;
    resnorm(i)=resnorm;
    exitflag(i)=exitflag;
    output(i)=output;
end

% TT = TEMP_K'-328;
% 
% for i=1:length(testset)
%     p = polyfit(TT(1:2:5),calcE(1:2:5,i),1);
%     pder=polyder(p);  % derivative of the linear function dG
%     tdsvec(i)=-polyval(pder,0)*328;    % Evaluationg entropy, dS at 298K = 24.85C in cal/mol/K
% end

dg_rms=rms(calcE-refE);



%output_name='Run_Oct_Paluch_thermo_new';
%save(output_name,'mol_list','calcE','refE','es','np','TEMP_K','x','dg_rms','dGfunc','dsvec','cpvec','index','resnorm','residual','output','exitflag');



    