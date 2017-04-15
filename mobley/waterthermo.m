clear all
close all

Home = getenv('HOME');
addpath(sprintf('%s/Research/pointbem',Home));
addpath(sprintf('%s/Research/panelbem',Home));
addpath(sprintf('%s/Research/testasymmetry',Home));
addpath(sprintf('%s/Research/testasymmetry/functions',Home));
addpath(sprintf('%s/Research/testasymmetry/mobley',Home));
addpath(sprintf('%s/Research/testasymmetry/born',Home));

     testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane','Li','Na','K','Rb','Cs','F','Cl','Br','I'};


tempdiv=5;
TEMP=linspace(5,45,tempdiv);
thermofunc=struct();

loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
logfileName = 'water.out';
epsOutfunc = @(x) (-1.410e-6)*x+(9.398e-4)*x^2-0.40008*x+87.740;

epsOut = epsOutfunc(TEMP(1));

ParamWatInfo = load('OptWater');
x = ParamWatInfo.xvec;
calcvec = ParamWatInfo.calcvec;

for i=1:length(testset)
    f=fit(TEMP',calcvec(:,i),'poly1');
    thermofunc(i).name=testset(i);
    thermofunc(i).func=f;
end

for i=1:length(testset)
    dgvec(i)=thermofunc(i).func(24.85);
    p=[thermofunc(i).func.p1,thermofunc(i).func.p2];
%     p=[thermofunc(i).func.p1,thermofunc(i).func.p2,thermofunc(i).func.p3,thermofunc(i).func.p4,thermofunc(i).func.p5];
    pder=polyder(p);
    dsvec(i)=-polyval(pder,24.85)*1000;
%     figure(i);
%     plot(thermofunc(i).func,TEMP',calcvec(:,i))   
end
dG_list_ref_at_298=[8.1,-40.5,-5.2,9,9.5,-6.2,-3.2,-21.2,-20.4,-24.6,-25.6,8.3,-529,-424,-352,-329,-306,-429,-304,-278,-243];%Hess
H_list_ref_at_298=[-8.3,-67.0,-23.9,-17.1,-17.1,-34.6,-25.3,-43.0,-45,-58.8,-57.4,-13.7];
dS_list_ref_at_298=(H_list_ref_at_298-dG_list_ref_at_298(1:12))/298;
dS_list_ref_ion_at_298=[-0.164;-0.133;-0.096;-0.087;-0.081;-0.115;-0.053;-0.037;-0.014];
dS_list_ref_at_298=[[dS_list_ref_at_298],dS_list_ref_ion_at_298'];
    
dG_list_ref_at_298
dgvec*4.184
dgerr=abs((dG_list_ref_at_298-dgvec*4.184))
dS_list_ref_at_298*1000
dsvec*4.184
dserr=abs((dS_list_ref_at_298*1000-dsvec*4.184))

figure(1)
scatter(dG_list_ref_at_298,dS_list_ref_at_298*1000,100*ones(1,length(testset)),'x','linewidth',2);
hold on
scatter(dgvec*4.184,dsvec*4.184,100*ones(1,length(testset)),'o','linewidth',2);
hold on
line([-600 100],[0 0])
ax=gca;
ax.LineWidth=1;
ax.Box='on';
set(gca,'fontsize',30);

    
    
    