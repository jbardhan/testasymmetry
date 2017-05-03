clear all
close all
clc

Home = getenv('HOME');

repo_path=sprintf('%s/Research',Home);
dropbox_path=sprintf('%s/Dropbox',Home);

addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path)); 


loadConstants
convertKJtoKcal = 1/joulesPerCalorie;

ParamWatInfo = load('Optwater_thermo_50_percent_1');
TEMP= ParamWatInfo.tempvec;
[m,index]=ismember(24.85,TEMP);

x = ParamWatInfo.xvec;
calcE = ParamWatInfo.calcvec;
refE = ParamWatInfo.dG_list_ref;
randE = ParamWatInfo.dG_list_rand;
testset = ParamWatInfo.testset;
refS = ParamWatInfo.dS_list_ref'*1000;
randS= ParamWatInfo.dS_list_rand'*1000;

refCP= ParamWatInfo.CP_list_ref'*1000;
randCP= ParamWatInfo.CP_list_rand'*1000;

ds_rand=[ParamWatInfo.ds_rand_coef_aca;ParamWatInfo.ds_rand_coef_ion]';

ionflag = ParamWatInfo.ionflag;
aca_num = ParamWatInfo.aca_num;
ion_num = ParamWatInfo.ion_num;
t_ref_ion=ParamWatInfo.t_ref_ion;
t_ref_aca=ParamWatInfo.t_ref_aca;

TEMP_K=TEMP+273.15;



dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures



for i=1:length(testset)
    f = @(R) (R(1)-R(2)*(TEMP_K-298)+R(3)*((TEMP_K-298)-TEMP_K.*log(TEMP_K./298)))-calcE(:,i);
    R0=[refE(i),refS(i)/1000,refCP(i)/1000];
    options=optimoptions('lsqnonlin','StepTolerance',1e-12);
    options=optimoptions(options,'OptimalityTolerance',1e-12);
    options=optimoptions(options,'FunctionTolerance',1e-12);
    [R,resnorm,residual,exitflag,output]=lsqnonlin(f,R0,[],[],options);
%     figure()
%     plot(TEMP_K,calcE(:,i),'ko',TEMP_K,R(1)-R(2)*(TEMP_K-298)+R(3)*((TEMP_K-298)-TEMP_K.*log(TEMP_K./298)),'b-')
    dGfunc(i).name=testset(i); 
    dGfunc(i).dg=R(1);
    dGfunc(i).ds=R(2);
    dGfunc(i).cp=R(3); 
    calcS(i)=dGfunc(i).ds*1000;
    calcCP(i)=dGfunc(i).cp*1000;
end



dg_rms_298=rms(refE(index,:)-calcE(index,:));
ds_rms_298=rms(refS-calcS);
cp_rms_298=rms(refCP-calcCP);


save('RunWater_training_thermo_50_percent_1','x','calcE','refE','refCP','randCP','testset','refS','randS','TEMP','ionflag','aca_num','ion_num','dGfunc','calcS','calcCP','t_ref_ion','t_ref_aca','dg_rms_298','ds_rms_298','cp_rms_298','index','ds_rand');

   
    
    
    