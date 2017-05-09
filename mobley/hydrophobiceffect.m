clear all
clc
close all
% for cas=1:2
cas=2;
    if cas==1
        RunWatInfo=load('RunWater_MNSol_MobleySurf_184_allthermo.mat');
    elseif cas==2
        RunWatInfo=load('RunWater_mobley_MobleySurf_513_allthermo.mat');
    end
    calcE=RunWatInfo.calcE;
    refE=RunWatInfo.refE;
    TEMP=RunWatInfo.TEMP;
    mol_list=RunWatInfo.mol_list;
    dGfunc=RunWatInfo.dGfunc;
    dsvec=RunWatInfo.dsvec;

    fid = fopen('mnsol/mobley_dG_AND_sa_and_vol.csv','r');
    Data = textscan(fid,'%s %f  %f  %f  %f  %f  %f  %f','delimiter',',');
    fclose(fid);
    all_solutes = Data{1};
    all_surfAreas = Data{3};
    all_volume=Data{4};
    all_mobley = Data{8};
    all_experiment = Data{2};
    [m, index] = ismember(mol_list,all_solutes);
    surfArea_list = all_surfAreas(index);
    volume_list = all_volume(index);
%%
%     figure()
%     scatter(surfArea_list(1:502),-dsvec(1:502)*298,'o')
%     ax=gca;
%     ax.LineWidth=1;
%     ax.Box='on';
%     set(gca,'fontsize',30);
%     xlabel('solute surface area');
%     ylabel('-T \DeltaS (cal/mol)');
%  %%   
%     figure()
%     scatter(volume_list(1:502),-dsvec(1:502)*298,'x')
%     ax=gca;
%     ax.LineWidth=1;
%     ax.Box='on';
%     set(gca,'fontsize',30);
%     xlabel('solute volume ');
%     ylabel('-T \DeltaS (cal/mol)');
% %%    
%     figure()
%     scatter(surfArea_list(1:502),calcE(2,1:502),'o')
%     ax=gca;
%     ax.LineWidth=1;
%     ax.Box='on';
%     set(gca,'fontsize',30);
%     xlabel('solute surface area');
%     ylabel('\DeltaG (kcal/mol)');
% %%
% 
%     figure()
%     scatter(volume_list(1:502),calcE(2,1:502),'x')
%     ax=gca;
%     ax.LineWidth=1;
%     ax.Box='on';
%     set(gca,'fontsize',30);
%     xlabel('solute volume ');
%     ylabel('\DeltaG (kcal/mol)');

%%
    figure()
    plot(refE(2,1:502),calcE(2,1:502),'*')
    hold on
    plot(refE(2,1:502),all_mobley(1:502),'o')
    ax=gca;
    ax.LineWidth=1;
    ax.Box='on';
    set(gca,'fontsize',30);
    xlabel('Experimental ');
    ylabel('Calculated \DeltaG (kcal/mol)'); 
    hold off
    legend('SLIC','Mobley','location','southeast');
    axis([min(refE(2,1:502))-3 max(refE(2,1:502))+3 min(refE(2,1:502))-3 max(refE(2,1:502))+3])
    hold on
    line([min(refE(2,1:502))-3, max(refE(2,1:502))+3],[min(refE(2,1:502))-3 ,max(refE(2,1:502))+3])

    
dg_rms_slic_vs_exp=rms(calcE(2,1:502)-refE(2,1:502))
dg_rms_slic_vs_mobley=rms(calcE(2,1:502)-all_mobley(1:502)')
dg_rms_mobley_vs_exp=rms(refE(2,1:502)-all_mobley(1:502)')




