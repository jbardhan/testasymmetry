clear all
Home = getenv('HOME');
resultsFolderPath = sprintf('%s/repos/testasymmetry/mobley/thermo_and_IL_results',Home);
figDirPath = sprintf('%s/figures',resultsFolderPath);
Data_all_thermo_vars_109 = readtable(sprintf('%s/all_109_solutes_calc_vs_expt.csv',resultsFolderPath));
Data_entropy_all = readtable(sprintf('%s/tds_calc_vs_expt.csv',resultsFolderPath));
Data_entropy_no_alcohols = readtable(sprintf('%s/tds_calc_vs_expt_slic_wo_ion_no_alkanols.csv',resultsFolderPath));
Data_cp = readtable(sprintf('%s/cp_calc_vs_expt_slic_wo_ion.csv',resultsFolderPath));
training_set = ["methane","ethanamide","methanethiol","n_butane","2_methylpropane",...
                "methyl_ethyl_sulfide","toluene","methanol","ethanol",...
                "3_methyl_1h_indole","p_cresol","propane"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ENTROPY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('Renderer', 'painters', 'Position', [12 8 1200 800]);
  plot(Data_entropy_all.Var2(1:158),Data_entropy_all.Var3(1:158),'ro','markers',12,'linewidth',2)
  xmax = max([max(Data_entropy_all.Var3(1:158)) max(Data_entropy_all.Var2(1:158))]);
  xmin = min([min(Data_entropy_all.Var3(1:158)) min(Data_entropy_all.Var2(1:158))]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel("$-T\Delta\,S_{expt}^{solv}$ (kcal/mol)","Interpreter","latex")
  ylabel("$-T\Delta\,S_{calc}^{solv}$ (kcal/mol)","Interpreter","latex")
  rmse = rms(Data_entropy_all.Var2(1:158)-Data_entropy_all.Var3(1:158));
  rmsErr = sprintf('%4.2f',rmse);
  num_data_points = length(Data_entropy_all.Var3(1:158));
  title({"Solvation Entropies - "+num_data_points+" data points"},'FontSize',18)
  legend("SLIC-RMSE = "+rmsErr+" kcal/mol")
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,"entropy_wo_ion_all");
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')
fig = figure('Renderer', 'painters', 'Position', [12 8 1200 800]);
  plot(Data_entropy_no_alcohols.Var2,Data_entropy_no_alcohols.Var3,'bs','markers',12,'linewidth',2)
  xmax = max([max(Data_entropy_no_alcohols.Var3) max(Data_entropy_no_alcohols.Var2)]);
  xmin = min([min(Data_entropy_no_alcohols.Var3) min(Data_entropy_no_alcohols.Var2)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel("$-T\Delta\,S_{expt}^{solv}$ (kcal/mol)","Interpreter","latex")
  ylabel("$-T\Delta\,S_{calc}^{solv}$ (kcal/mol)","Interpreter","latex")
  rmse = rms(Data_entropy_no_alcohols.Var2-Data_entropy_no_alcohols.Var3);
  rmsErr = sprintf('%4.2f',rmse);
  num_data_points = length(Data_entropy_no_alcohols.Var3);
  title({"Solvation Entropies (alcohols not included in the data set)- "+num_data_points+" data points"},'FontSize',18)
  legend("SLIC-RMSE = "+rmsErr+" kcal/mol")
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,"entropy_wo_ion_no_alkanols");
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')
fig = figure('Renderer', 'painters', 'Position', [12 8 1200 800]);
  plot(Data_cp.Var2,Data_cp.Var3,'bs','markers',12,'linewidth',2)
  xmax = max([max(Data_cp.Var3) max(Data_cp.Var2)]);
  xmin = min([min(Data_cp.Var3) min(Data_cp.Var2)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel("$C_{p\,\,expt}^{\,\,\,solv}$ (cal/mol/K)","Interpreter","latex")
  ylabel("$C_{p\,\,calc}^{\,\,\,solv}$ (cal/mol/K)","Interpreter","latex")
  rmse = rms(Data_cp.Var2-Data_cp.Var3);
  rmsErr = sprintf('%4.2f',rmse);
  legend("SLIC-RMSE = "+rmsErr+" kcal/mol")
  num_data_points = length(Data_cp.Var3);
  title({"Solvation Heat Capacitities - "+num_data_points+" data points"},'FontSize',18)
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,"heat_capacity");
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')
fig = figure('Renderer', 'painters', 'Position', [12 8 800 600]);
  Data = Data_all_thermo_vars_109;
  
  % Data.Var1   |   Var2   |   Var3   |   Var4   |   Var5 
  
  %     solutes    dG_ref     dH_ref    -TdS_ref    Cp_ref
  
  %             |   Var6   |   Var7   |   Var8   |   Var9 
 
  %                dG_calc    dH_calc   -TdS_calc   Cp_calc
  [m, index] = ismember(training_set,Data.Var1);
  subplot(2,2,1)
  plot(Data.Var2,Data.Var6,'bs','markers',6,'linewidth',2)
  hold on
  plot(Data.Var2(index),Data.Var6(index),'rd','markers',6,'linewidth',2)
  xmax = max([max(Data.Var6) max(Data.Var2)]);
  xmin = min([min(Data.Var6) min(Data.Var2)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  rmse = rms(Data.Var2-Data.Var6);
  rmsErr = sprintf('%4.2f',rmse);
  rmse_tr = rms(Data.Var2(index)-Data.Var6(index));
  rmsErrTr = sprintf('%4.2f',rmse_tr);
  xlabel("$\Delta\,G_{\,expt}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  ylabel("$\Delta\,G_{\,calc}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  title('Solvation Energies')
  legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},'FontSize',8)
  set(legend,'location','northwest')
  subplot(2,2,2)
  plot(Data.Var3,Data.Var7,'bs','markers',6,'linewidth',2)
  hold on
  plot(Data.Var3(index),Data.Var7(index),'rd','markers',6,'linewidth',2)
  xmax = max([max(Data.Var7) max(Data.Var3)]);
  xmin = min([min(Data.Var7) min(Data.Var3)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  rmse = rms(Data.Var3-Data.Var7);
  rmsErr = sprintf('%4.2f',rmse);
  rmse_tr = rms(Data.Var3(index)-Data.Var7(index));
  rmsErrTr = sprintf('%4.2f',rmse_tr);
  xlabel("$\Delta\,H_{\,expt}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  ylabel("$\Delta\,H_{\,calc}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  title('Solvation Enthalpies')
  legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},'FontSize',8)
  set(legend,'location','northwest')
  subplot(2,2,3)
  plot(Data.Var4,Data.Var8,'bs','markers',6,'linewidth',2)
  hold on
  plot(Data.Var4(index),Data.Var8(index),'rd','markers',6,'linewidth',2)
  xmax = max([max(Data.Var8) max(Data.Var4)]);
  xmin = min([min(Data.Var8) min(Data.Var4)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  rmse =rms(Data.Var4-Data.Var8);
  rmsErr = sprintf('%4.2f',rmse);
  rmse_tr = rms(Data.Var4(index)-Data.Var8(index));
  rmsErrTr = sprintf('%4.2f',rmse_tr);
  xlabel("$-T\Delta\,S_{\,expt}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  ylabel("$-T\Delta\,S_{\,calc}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  title('Solvation Entropies')
  legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},'FontSize',8)
  set(legend,'location','northwest')
  subplot(2,2,4)
  plot(Data.Var5,Data.Var9,'bs','markers',6,'linewidth',2)
  hold on
  plot(Data.Var5(index),Data.Var9(index),'rd','markers',6,'linewidth',2)
  xmax = max([max(Data.Var9) max(Data.Var5)]);
  xmin = min([min(Data.Var9) min(Data.Var5)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  rmse = rms(Data.Var5-Data.Var9);
  rmse_tr = rms(Data.Var5(index)-Data.Var9(index));
  rmsErr = sprintf('%4.2f',rmse);
  rmsErrTr = sprintf('%4.2f',rmse_tr);
  xlabel("$\Delta\,C_{p\,\,expt}^{\,\,\,\,solv}$ (kcal/mol/K)", ...
         "Interpreter","latex")
  ylabel("$\Delta\,C_{p\,\,calc}^{\,\,\,\,solv}$ (kcal/mol/K)", ...
         "Interpreter","latex")
  title('Solvation Heat Capacities')
  legend({"RMSE Test-set = "+rmsErr+" kcal/mol/K","RMSE Training-set = "+rmsErrTr+" kcal/mol/K"},'FontSize',8)
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,"thermo_all_109");
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')
fig = figure('Renderer', 'painters', 'Position', [12 8 1200 800]);
  plot(Data.Var2,Data.Var6,'bs','markers',12,'linewidth',2)
  hold on
  %Data.Var2(index)=NaN;
  %Data.Var6(index)=NaN;
  plot(Data.Var2(index),Data.Var6(index),'r*','markers',12,'linewidth',2)
  xmax = max([max(Data.Var6) max(Data.Var2)]);
  xmin = min([min(Data.Var6) min(Data.Var2)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  rmse = rms(Data.Var2-Data.Var6);
  rmsErr = sprintf('%4.2f',rmse);
  rmse_tr = rms(Data.Var2(index)-Data.Var6(index));
  rmsErrTr = sprintf('%4.2f',rmse_tr);
  xlabel("$\Delta\,G_{\,expt}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  ylabel("$\Delta\,G_{\,calc}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  title('Solvation Free Energies','FontSize',18)
  legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},'FontSize',16)
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,"dG_all_109");
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')
fig = figure('Renderer', 'painters', 'Position', [12 8 1200 800]);
  plot(Data.Var3,Data.Var7,'bs','markers',12,'linewidth',2)
  hold on
  plot(Data.Var3(index),Data.Var7(index),'r*','markers',12,'linewidth',2)
  xmax = max([max(Data.Var7) max(Data.Var3)]);
  xmin = min([min(Data.Var7) min(Data.Var3)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  rmse = rms(Data.Var3-Data.Var7);
  rmsErr = sprintf('%4.2f',rmse);
  rmse_tr = rms(Data.Var3(index)-Data.Var7(index));
  rmsErrTr = sprintf('%4.2f',rmse_tr);
  xlabel("$\Delta\,H_{\,expt}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  ylabel("$\Delta\,H_{\,calc}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  title('Solvation Enthalpies','FontSize',18)
  legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},'FontSize',16)
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,"dH_all_109");
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')
fig = figure('Renderer', 'painters', 'Position', [12 8 1200 800]);
  plot(Data.Var4,Data.Var8,'bs','markers',12,'linewidth',2)
  hold on
  plot(Data.Var4(index),Data.Var8(index),'r*','markers',12,'linewidth',2)
  xmax = max([max(Data.Var8) max(Data.Var4)]);
  xmin = min([min(Data.Var8) min(Data.Var4)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  rmse =rms(Data.Var4-Data.Var8);
  rmsErr = sprintf('%4.2f',rmse);
  rmse_tr = rms(Data.Var4(index)-Data.Var8(index));
  rmsErrTr = sprintf('%4.2f',rmse_tr);
  xlabel("$-T\Delta\,S_{\,expt}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  ylabel("$-T\Delta\,S_{\,calc}^{\,solv}$ (kcal/mol)", ...
         "Interpreter","latex")
  title('Solvation Entropies','FontSize',18)
  legend({"RMSE Test-set = "+rmsErr+" kcal/mol","RMSE Training-set = "+rmsErrTr+" kcal/mol"},'FontSize',16)
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,"TdS_all_109");
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')
fig = figure('Renderer', 'painters', 'Position', [12 8 1200 800]);
  plot(Data.Var5,Data.Var9,'bs','markers',12,'linewidth',2)
  hold on
  plot(Data.Var5(index),Data.Var9(index),'r*','markers',12,'linewidth',2)
  xmax = max([max(Data.Var9) max(Data.Var5)]);
  xmin = min([min(Data.Var9) min(Data.Var5)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  rmse = rms(Data.Var5-Data.Var9);
  rmse_tr = rms(Data.Var5(index)-Data.Var9(index));
  rmsErr = sprintf('%4.2f',rmse);
  rmsErrTr = sprintf('%4.2f',rmse_tr);
  xlabel("$\Delta\,C_{p\,\,expt}^{\,\,\,\,solv}$ (kcal/mol/K)", ...
         "Interpreter","latex")
  ylabel("$\Delta\,C_{p\,\,calc}^{\,\,\,\,solv}$ (kcal/mol/K)", ...
         "Interpreter","latex")
  title('Solvation Heat Capacities','FontSize',18)
  legend({"RMSE Test-set = "+rmsErr+" kcal/mol/K","RMSE Training-set = "+rmsErrTr+" kcal/mol/K"},'FontSize',16)
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,"Cp_all_109");
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')