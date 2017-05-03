clear all
%clc
% close all

%Water
%%
trainingplot_flag=1;
mnsolplot_flag=0;
mobleyplot_flag=0;

outlier_flag=0;
%%
run_water_training=load('RunWater_training_thermo_100_percent_1.mat');

index_training=run_water_training.index;  % index of 298K =24.85C in the temp vector
dg_ref_training=run_water_training.refE(index_training,:);   % expaerimental Delta_G of the training set in kcal/mol
ds_ref_training=run_water_training.refS;  % expaerimentalDelta S of the training set in cal/mol/K
cp_ref_training=run_water_training.refCP;  % expaerimentalDelta S of the training set in cal/mol/K
dg_calc_training=run_water_training.calcE(index_training,:); % calculated Delta_G of the training set in kcal/mol
ds_calc_training=run_water_training.calcS; % calculated S of the training set in cal/mol/K
cp_calc_training=run_water_training.calcCP; % calculated S of the training set in cal/mol/K
dg_rms_298_training=run_water_training.dg_rms_298;
ds_rms_298_training=run_water_training.ds_rms_298;
cp_rms_298_training=run_water_training.cp_rms_298;
ds_rand_training=run_water_training.randS; % calculated S of the training set in cal/mol/K
cp_rand_training=run_water_training.randCP; % calculated S of the training set in cal/mol/K

ds_rand_num=run_water_training.ds_rand;

if outlier_flag==1

    dg_err_training=abs(dg_ref_training-dg_calc_training);
    ds_err_training=abs(ds_ref_training-ds_calc_training);
    cp_err_training=abs(cp_ref_training-cp_calc_training);

    dg_err_training=sort(dg_err_training);
    ds_err_training=sort(ds_err_training);
    cp_err_training=sort(cp_err_training);

    dg_err_training_outlier=dg_err_training(end-1:end);
    ds_err_training_outlier=ds_err_training(end-1:end);
    cp_err_training_outlier=cp_err_training(end-1:end);
    
end

    dgcorrcoef=corrcoef(dg_ref_training,dg_calc_training);
    dscorrcoef=corrcoef(ds_ref_training,ds_calc_training);
    cpcorrcoef=corrcoef(cp_ref_training,cp_calc_training);
    





%%
run_water_MNSol=load('RunWater_mnsol_mobleysurf_thermo.mat');

index_mnsol=run_water_MNSol.index;
dg_ref_mnsol=run_water_MNSol.refE(index_mnsol,:);
dg_calc_slic_mnsol=run_water_MNSol.calcE(index_mnsol,:);
ds_calc_slic_mnsol=run_water_MNSol.dsvec;
dg_rms_slic_mnsol=run_water_MNSol.dg_rms_298;

%%
run_water_Mobley=load('RunWater_mobley_thermo.mat');

index_mobley=run_water_Mobley.index;
dg_ref_mobley=run_water_Mobley.refE(index_mobley,:);
dg_calc_slic_mobley=run_water_Mobley.calcE(index_mobley,:);
ds_calc_slic_mobley=run_water_Mobley.dsvec;
dg_calc_mobley=run_water_Mobley.calc_mobley;
dg_rms_slic_mobley=run_water_Mobley.dg_rms_298;
dg_rms_mobley=run_water_Mobley.dg_rms_298_MD;



%%
if trainingplot_flag==1
    figure()
    p=plot(dg_ref_training,dg_calc_training,'o');
    plotsettings(p,'b'); % b is the color of markers
    xlabel('Experimental');
    ylabel('Computed');
    ax=gca;
    ax.XLim=[-120 20];
    ax.YLim=[-120 20];
    diagline=refline(1,0);
    set(diagline,'LineWidth',2);
    set(diagline,'Color','k');
    leg=legend(['SLIC \DeltaG (kcal/mol);  Training set; RMS = ',num2str(dg_rms_298_training)],'location','southeast');
    leg.Location='southeast';
    leg.FontSize=20;
    hold off

    figure()
    p=plot(ds_ref_training,ds_calc_training,'o');
    %p=plot(ds_rand_training,ds_calc_training,'x');
    plotsettings(p,'b') % b is the color of markers
    xlabel('Experimental \DeltaS (cal/mol K)');
    ylabel('Computed \DeltaS (cal/mol K)'); 
    diagline=refline(1,0);
    set(diagline,'LineWidth',2);
    set(diagline,'Color','k');
    leg=legend(['SLIC \DeltaS (cal/mol^\circK); training set; RMS = ',num2str(ds_rms_298_training)]);
    leg.Location='southeast';
    leg.FontSize=20;
    hold off
    
    
    figure()
    p=plot(cp_ref_training,cp_calc_training,'o');
    %p=plot(cp_rand_training,cp_calc_training,'x');
    plotsettings(p,'b') % b is the color of markers
    xlabel('Experimental CP (cal/mol K)');
    ylabel('Computed CP (cal/mol K)'); 
    diagline=refline(1,0);
    set(diagline,'LineWidth',2);
    set(diagline,'Color','k');
    leg=legend(['SLIC CP (cal/mol^\circK); training set; RMS = ',num2str(cp_rms_298_training)]);
    leg.Location='southeast';
    leg.FontSize=20;
    hold off
end

%%
if mnsolplot_flag==1
    figure()
    p=plot(dg_ref_mnsol,dg_calc_slic_mnsol,'o');
    plotsettings(p,'b'); % b is the color of markers
    xlabel('Experimental');
    ylabel('Computed'); 
    hold off
    diagline=refline(1,0);
    set(diagline,'LineWidth',2);
    set(diagline,'Color','k');
    leg=legend(['SLIC \DeltaG (kcal/mol) for the MNSol database RMS error = ',num2str(dg_rms_slic_mnsol)]);
    leg.Location='southeast';
    leg.FontSize=20;
end

%%
if mobleyplot_flag==1
    figure()
    p=plot(dg_ref_mobley,dg_calc_slic_mobley,'o');
    plotsettings(p,'b'); % b is the color of markers
    hold on
    q=plot(dg_ref_mobley,dg_calc_mobley,'s');
    plotsettings(q,'r'); % b is the color of markers
    hold on
    leg=legend( ['SLIC \DeltaGs (kcal/mol); Mobley database; RMS = ',num2str(dg_rms_slic_mobley)],...
                ['MD \DeltaGs (kcal/mol); Mobley database; RMS =',num2str(dg_rms_mobley)]);
    leg.Location='southeast';
    leg.FontSize=20;
    
    diagline=refline(1,0);
    set(diagline,'LineWidth',2);
    set(diagline,'Color','k');
    xlabel('Experimental');
    ylabel('Computed'); 
    hold off
end

% axis([min(dg_ref_training)-3 max(dg_ref_training)+3 min(dg_ref_training)-3 max(dg_ref_training)+3])
% hold on
%plot([min(dg_ref_training)-3, max(dg_ref_training)+3],[min(dg_ref_training)-3 ,max(dg_ref_training)+3],'k','linewidth',2)



