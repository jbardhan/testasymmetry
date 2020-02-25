title_font = 'Helvetica';
legend_font = 'Helvetica';
label_font = 'Helvetica';
title_font_size = 18;
legend_font_size = 14;
label_font_size = 18;
data = load('Run_Cosmo_Mob_SLIC.mat');
[~, index] = ismember(data.testset,data.mol_list);
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
num_tr_data = length(data.testset);
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
legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_expt.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% NONPOLAR
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
np_cosmo = data.comb+data.disp_slsl-data.disp_slsv;
plot(data.np_mobley,np_cosmo,'bo','markers',12,'linewidth',2);
hold on
plot(data.np_mobley(index),data.np(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.np_mobley) max(data.np)]);
xmin = min([min(data.np_mobley) min(data.np)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.testset);
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
rmse = rms(data.np_mobley-data.np);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.np_mobley(index)-data.np(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_mob_np.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')

%% Dispersion
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
plot(data.disp_mobley,data.disp_slsv,'bo','markers',12,'linewidth',2);
hold on
plot(data.disp_mobley(index),data.disp_slsv(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.disp_mobley) max(data.disp_slsv)]);
xmin = min([min(data.disp_mobley) min(data.disp_slsv)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.testset);
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
rmse = rms(data.disp_mobley-data.disp_slsv);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.disp_mobley(index)-data.disp_slsv(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_mob_disp.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% Cavity
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
cav_cosmo = -(data.comb+data.disp_slsl);
plot(data.cav_mobley,cav_cosmo,'bo','markers',12,'linewidth',2);
hold on
plot(data.cav_mobley(index),cav_cosmo(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.cav_mobley) max(cav_cosmo)]);
xmin = min([min(data.cav_mobley) min(cav_cosmo)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.testset);
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
rmse = rms(data.cav_mobley-cav_cosmo);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.cav_mobley(index)-cav_cosmo(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_mob_cav.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% ES
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
plot(data.es_mobley,data.es,'bo','markers',12,'linewidth',2);
hold on
plot(data.es_mobley(index),data.es(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.es_mobley) max(data.es)]);
xmin = min([min(data.es_mobley) min(data.es)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.testset);
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
rmse = rms(data.es_mobley-data.es);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.es_mobley(index)-data.es(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_es_np.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% ES + HB
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
es = data.es+data.hb;
plot(data.es_mobley,es,'bo','markers',12,'linewidth',2);
hold on
plot(data.es_mobley(index),es(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.es_mobley) max(es)]);
xmin = min([min(data.es_mobley) min(es)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.testset);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,MD}^{\,solv,es+HB}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{COSMO}}^{\,solv,es+HB}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. COSMO predicted electrostatic (+HB) energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.es_mobley-es);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.es_mobley(index)-es(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('cosmo_vs_eshb_np.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
%% ES_original_SLIC vs Mobley
fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
es = data.es_SLIC;
plot(data.es_mobley,es,'bo','markers',12,'linewidth',2);
%hold on
%plot(data.es_mobley(index),es(index),'r*','markers',10,'linewidth',3)
xmax = max([max(data.es_mobley) max(es)]);
xmin = min([min(data.es_mobley) min(es)]);
axis([xmin-1 xmax+1 xmin-1 xmax+1]);
set(gca,'FontSize',16)
foo = refline(1,0);
set(foo,'Linewidth',2,'color','k');
num_tr_data = length(data.testset);
num_data = length(data.mol_list);
xlabel("$\Delta\,G_{\,MD}^{\,solv}$ " + "water" + ...
           " (kcal/mol)","Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
ylabel("$\Delta\,G_{\,calc,\mathrm{SLIC}}^{\,solv,es}$ " + "water" + ...
       " (kcal/mol)","Interpreter","latex",...
       'FontName',label_font,'FontSize',label_font_size)

title({"MD (Mobley) vs. SLIC$ predicted electrostatic energies of "+num_data...
           "solutes in water"+" using a training set of length "+num_tr_data},...
           'FontName',label_font,'FontSize',title_font_size)
rmse = rms(data.es_mobley-es);
rmsErr = sprintf('%4.2f',rmse);
rmse_tr = rms(data.es_mobley(index)-es(index));
rmsErrTr = sprintf('%4.2f',rmse_tr);
legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},...
      'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
set(legend,'location','northwest')
fileName= sprintf('SLIC_vs_mob_es.pdf');
orient(fig,'landscape')
print(fig,fileName,'-dpdf')
