title_font = 'Helvetica';
legend_font = 'Helvetica';
label_font = 'Helvetica';
title_font_size = 18;
legend_font_size = 14;
label_font_size = 18;
data = load('RunCosmoFixed.mat');

[~, index] = ismember(data.training_set,data.mol_list);
close all
%% TOTAL SOLVATION ENERGY
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
plot(data.ref,data.calc,'bo','markers',12,'linewidth',2);
hold on
plot(data.ref(index),data.calc(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.ref) max(data.calc)]);
xmin = min([min(data.ref) min(data.calc)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.training_set);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,expt}^{\,solv}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"Experimental vs. COSMO predicted solvation energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.ref-data.calc);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.ref(index)-data.calc(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Training-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_expt.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% NONPOLAR
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
np_cosmo = data.np;
plot(data.np_mob,np_cosmo,'bo','markers',12,'linewidth',2);
hold on
plot(data.np_mob(index),np_cosmo(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.np_mob) max(np_cosmo)]);
xmin = min([min(data.np_mob) min(np_cosmo)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.training_set);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,MD}^{\,solv,np}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,np}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. COSMO predicted nonpolar energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.np_mob-np_cosmo);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.np_mob(index)-np_cosmo(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_mob_np.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')

%% Dispersion
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
disp_cosmo = data.disp_svsv - 2*data.disp_svsl;
plot(data.disp_mob,disp_cosmo,'bo','markers',12,'linewidth',2);
hold on
plot(data.disp_mob(index),disp_cosmo(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.disp_mob) max(disp_cosmo)]);
xmin = min([min(data.disp_mob) min(disp_cosmo)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.training_set);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,MD}^{\,solv,disp}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,disp}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. COSMO predicted dispersion energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.disp_mob-disp_cosmo);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.disp_mob(index)-disp_cosmo(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Training-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_mob_disp.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% Cavity
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
cav_cosmo = (data.comb-data.disp_slsl);
plot(data.cav_mob,cav_cosmo,'bo','markers',12,'linewidth',2);
hold on
plot(data.cav_mob(index),cav_cosmo(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.cav_mob) max(cav_cosmo)]);
xmin = min([min(data.cav_mob) min(cav_cosmo)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.training_set);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,MD}^{\,solv,cav}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,cav}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. COSMO predicted cavity energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.cav_mob-cav_cosmo);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.cav_mob(index)-cav_cosmo(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Training-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_mob_cav.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% ES
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
plot(data.es_mob,data.es,'bo','markers',12,'linewidth',2);
hold on
plot(data.es_mob(index),data.es(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.es_mob) max(data.es)]);
xmin = min([min(data.es_mob) min(data.es)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.training_set);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,MD}^{\,solv,es}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,es}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. COSMO predicted electrostatic energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.es_mob-data.es);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.es_mob(index)-data.es(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Training-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_es_np.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% ES + HB
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
es = data.es+data.hb;
plot(data.es_mob,es,'bo','markers',12,'linewidth',2);
hold on
plot(data.es_mob(index),es(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.es_mob) max(es)]);
xmin = min([min(data.es_mob) min(es)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.training_set);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,MD}^{\,solv,es}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,es+HB}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. COSMO predicted electrostatic (+HB) energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.es_mob-es);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.es_mob(index)-es(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Training-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_eshb_np.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% ES_original_SLIC vs Mobley
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
es = data.es_SLIC;
plot(data.es_mob,es,'bo','markers',12,'linewidth',2);
xmax = max([max(data.es_mob) max(es)]);
xmin = min([min(data.es_mob) min(es)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.training_set);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,MD}^{\,solv,es}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{SLIC}}^{\,solv,es}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. SLIC$ predicted electrostatic energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.es_mob-es);
rmsErr = sprintf('%4.2f',rmse);
legend({"RMSE Training-set = "+rmsErr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('SLIC_vs_mob_es.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% ES_original_SLIC vs ES_cosmo_SLIC
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
es = data.es_SLIC;
plot(data.es,es,'bo','markers',12,'linewidth',2);
xmax = max([max(data.es) max(es)]);
xmin = min([min(data.es) min(es)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.training_set);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,es}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{SLIC}}^{\,solv,es}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. SLIC$ predicted electrostatic energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.es-es);
rmsErr = sprintf('%4.2f',rmse);
legend({"RMSE Training-set = "+rmsErr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('SLIC_vs_cosmo_es.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
% %% Dispersion slsl vs slsv
% fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
% plot(data.disp_slsl,data.disp_slsv,'bo','markers',12,'linewidth',2);
% hold on
% plot(data.disp_slsl(index),data.disp_slsv(index),'r*','markers',10,'linewidth',3)
% xmax = max([max(data.disp_slsl) max(data.disp_slsv)]);
% xmin = min([min(data.disp_slsl) min(data.disp_slsv)]);
% axis([xmin-1 xmax+1 xmin-1 xmax+1]);
% set(gca,'FontSize',16)
% foo = refline(1,0);
% set(foo,'Linewidth',2,'color','k');
% num_tr_data = length(data.training_set);
% num_data = length(data.mol_list);
% xlabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,disp-slsl}$ " + "water" + ...
%            " (kcal/mol)","Interpreter","latex",...
%            'FontName',label_font,'FontSize',label_font_size)
% ylabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,disp-slsv}$ " + "water" + ...
%        " (kcal/mol)","Interpreter","latex",...
%        'FontName',label_font,'FontSize',label_font_size)
% 
% title({"MD (Mobley) vs. COSMO predicted dispersion energies of "+num_data...
%            "solutes in water"+" using a training set of length "+num_tr_data},...
%            'FontName',label_font,'FontSize',title_font_size)
% rmse = rms(data.disp_slsl-data.disp_slsv);
% rmsErr = sprintf('%4.2f',rmse);
% rmse_tr = rms(data.disp_slsl(index)-data.disp_slsv(index));
% rmsErrTr = sprintf('%4.2f',rmse_tr);
% legend({"RMSE Training-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
%       'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
% set(legend,'location','northwest')
% fileName= sprintf('cosmo_disp.pdf');
% orient(fig,'landscape')
% print(fig,fileName,'-dpdf')