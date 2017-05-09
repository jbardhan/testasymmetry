clear all
close all
clc

loadConstants
convertKJtoKcal = 1/joulesPerCalorie;

ParamWatInfo = load('Optwater_thermo');
TEMP= ParamWatInfo.tempvec;
[m,index]=ismember(24.85,TEMP);

x = ParamWatInfo.xvec;
calcE = ParamWatInfo.calcvec;
refE = ParamWatInfo.refvec;
testset = ParamWatInfo.testset;
refS = ParamWatInfo.dS_list'*1000;

ionflag = ParamWatInfo.ionflag;
aca_num = ParamWatInfo.aca_num;
ion_num = ParamWatInfo.ion_num;
t_ref_ion=ParamWatInfo.t_ref_ion;
t_ref_aca=ParamWatInfo.t_ref_aca;




dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures



for i=1:length(testset)
    f=fit(TEMP,calcE(:,i),'poly1'); % the linear function that fits to the calculated dG  at differet temperatures
    dGfunc(i).name=testset(i); 
    dGfunc(i).func=f;
end

for i=1:length(testset)
    p=[dGfunc(i).func.p1,dGfunc(i).func.p2];
    pder=polyder(p);  % derivative of the linear function dG
    if i <= aca_num
        calcS(i)=-polyval(pder,t_ref_aca)*1000;    % Evaluationg entropy, dS at 298K = 24.85C in cal/mol/K
    else
        calcS(i)=-polyval(pder,t_ref_ion)*1000;    % Evaluationg entropy, dS at 298K = 24.85C in cal/mol/K
    end
%     figure(i);
%     plot(dGfunc(i).func,TEMP',calcvec(:,i),'o')   
end

dg_rms_298=rms(refE(index,:)-calcE(index,:));
ds_rms_298=rms(refS-calcS);

save('RunWater_training_thermo','x','calcE','refE','testset','refS','TEMP','ionflag','aca_num','ion_num','dGfunc','calcS','t_ref_ion','t_ref_aca','dg_rms_298','ds_rms_298','index');

   
    
    
    