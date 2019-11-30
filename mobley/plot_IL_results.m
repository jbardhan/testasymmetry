clear all
Home = getenv('HOME');
resultsFolderPath = sprintf('%s/repos/testasymmetry/mobley/thermo_and_IL_results',Home);
refDataPath = sprintf('%s/repos/testasymmetry/mobley/reference-data',Home);
figDirPath = sprintf('%s/figures',resultsFolderPath);
 
solvents = {'[Emim][Tf2N]_p_1','[Emim][Tf2N]_p_2','[Emim][Tf2N]_p_4',...
            '[Bmim][Tf2N]_p_2','[Bmim][Tf2N]_l_2',...
            '[Bmim][BF4]_l_1','[Bmim][BF4]_l_2','[Bmim][BF4]_l_3',...
            '[Bmim][PF6]_l_1','[Bmim][PF6]_l_2','[Bmim][PF6]_l_3',...
            '[Bmim][TfO]_l_2','[Bmim][TfO]_l_4',...
            '[Bmim][Cl]_l_2',...
            '[Omim][Tf2N]_p_2'};   
% solventName_referenceThermoData_referenceEpsilonData
% referenceThermoData: p=Paluch12, l=Latif, \varepsilon = '+eps4
% referenceEpsilonData: 1=Wakai05, 2=Rybinska14, 3=Hunger09, 4=Huang11
% See 'thermo_and_IL_results/IL_res.xls' fo, \varepsilon = r+eps more details
 
Data = readtable(sprintf('%s/IL_res.csv',resultsFolderPath));
Data_eps = readtable(sprintf('%s/eps_data.csv',refDataPath)); 

% for solvents{k}: Data.Var2k-1 = ref_data, Data.Var2k = calc_data
%
% example: solvents{6}=[Bmim][BF4]_l_1 is for [Bmim][BF4] with Latif, \varepsilon=' + eps4 as 
% the reference for dG and Wakai05 as the reference for epsilon
% 2k=12 ==> Data.Var11 = reference dG from Latif, \varepsilon='+eps% , \varepsilon =  +eps         Data.Var12 = calcula t ed dG

% Plot all data, \varepsilon = 
figure
  i = 1;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i));
  plot(Data.Var1,Data.Var2,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var1) max(Data.Var2)]);
  xmin = min([min(Data.Var1) min(Data.Var2)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = " + eps+")")
  else
      legend("SLIC-Paluch, \epsilon = " + eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 2;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i));
  plot(Data.Var3,Data.Var4,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var3) max(Data.Var4)]);
  xmin = min([min(Data.Var3) min(Data.Var4)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 3;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i));
  plot(Data.Var5,Data.Var6,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var5) max(Data.Var6)]);
  xmin = min([min(Data.Var5) min(Data.Var6)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 4;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i));
  plot(Data.Var7,Data.Var8,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var7) max(Data.Var8)]);
  xmin = min([min(Data.Var7) min(Data.Var8)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 5;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var9,Data.Var10,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var9) max(Data.Var10)]);
  xmin = min([min(Data.Var9) min(Data.Var10)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 6;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var11,Data.Var12,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var11) max(Data.Var12)]);
  xmin = min([min(Data.Var11) min(Data.Var12)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 7;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var13,Data.Var14,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var13) max(Data.Var14)]);
  xmin = min([min(Data.Var13) min(Data.Var14)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 8;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var15,Data.Var16,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var15) max(Data.Var16)]);
  xmin = min([min(Data.Var15) min(Data.Var16)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 9;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var17,Data.Var18,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var17) max(Data.Var18)]);
  xmin = min([min(Data.Var17) min(Data.Var18)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 10;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var19,Data.Var20,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var19) max(Data.Var20)]);
  xmin = min([min(Data.Var19) min(Data.Var20)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 11;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var21,Data.Var22,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var21) max(Data.Var22)]);
  xmin = min([min(Data.Var21) min(Data.Var22)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 12;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var23,Data.Var24,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var23) max(Data.Var24)]);
  xmin = min([min(Data.Var23) min(Data.Var24)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 13;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var25,Data.Var26,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var25) max(Data.Var26)]);
  xmin = min([min(Data.Var25) min(Data.Var26)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 14;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var27,Data.Var28,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var27) max(Data.Var28)]);
  xmin = min([min(Data.Var27) min(Data.Var28)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
figure
  i = 15;
  solvent=solvents{i};
  eps = sprintf('%4.2f',Data_eps.Var2(i-1));
  plot(Data.Var29,Data.Var30,'bo','markers',12,'linewidth',2)
  xmax = max([max(Data.Var29) max(Data.Var30)]);
  xmin = min([min(Data.Var29) min(Data.Var30)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  if solvent(length(solvent)-2)=='l'
      legend("SLIC-Latif, \epsilon = "+eps+")")
  else
      legend("SLIC-Paluch, \epsilon = "+eps+")")
  end
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s.pdf',figDirPath,solvent);
  saveas(gcf,fileName)
for i=1:max(size(Data_eps))
    epsilon{i} = sprintf('%4.2f',Data_eps.Var2(i));
end
figure
  solvent=solvents{1};
  plot(Data.Var1,Data.Var2,'bo','markers',12,'linewidth',2)
  hold on
  plot(Data.Var3,Data.Var4,'ro','markers',12,'linewidth',2)
  hold on
  plot(Data.Var5,Data.Var6,'go','markers',12,'linewidth',2)
  xmax = max([max(Data.Var1) max(Data.Var2)]);
  xmin = min([min(Data.Var1) min(Data.Var2)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  legend("Paluch (\epsilon = "+epsilon{1}+")",...
         "Paluch (\epsilon = "+epsilon{2}+")",...
         "Paluch (\epsilon = "+epsilon{3}+")")
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s_compare.pdf',figDirPath,solvent(1:length(solvent)-4));
  saveas(gcf,fileName)
figure
  solvent=solvents{4};
  plot(Data.Var7,Data.Var8,'bo','markers',12,'linewidth',2)
  hold on
  plot(Data.Var9,Data.Var10,'ro','markers',12,'linewidth',2)
  xmax = max([max(Data.Var7) max(Data.Var8) max(Data.Var10)]);
  xmin = min([min(Data.Var7) min(Data.Var8) min(Data.Var10)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  legend("Paluch (\epsilon = "+epsilon{4}+")",...
         "Latif (\epsilon = " +epsilon{4}+")")
  set(legend,'location','northwest')
  fileName=sprintf('%s/compare_%s_paluch_vs_latif.pdf',figDirPath,solvent(1:length(solvent)-4));
  saveas(gcf,fileName)
figure
  solvent=solvents{6};
  plot(Data.Var11,Data.Var12,'bo','markers',12,'linewidth',2)
  hold on
  plot(Data.Var13,Data.Var14,'ro','markers',12,'linewidth',2)
  hold on
  plot(Data.Var15,Data.Var16,'go','markers',12,'linewidth',2)
  xmax = max([max(Data.Var17) max(Data.Var18)]);
  xmin = min([min(Data.Var17) min(Data.Var18)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  legend("Latif (\epsilon = "+epsilon{5}+")",...
         "Latif (\epsilon = "+epsilon{6}+")",...
         "Latif (\epsilon = "+epsilon{7}+")")
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s_compare.pdf',figDirPath,solvent(1:length(solvent)-4));
  saveas(gcf,fileName)
figure
  solvent=solvents{8};
  plot(Data.Var17,Data.Var18,'bo','markers',12,'linewidth',2)
  hold on
  plot(Data.Var19,Data.Var20,'ro','markers',12,'linewidth',2)
  hold on
  plot(Data.Var21,Data.Var22,'go','markers',12,'linewidth',2)
  xmax = max([max(Data.Var19) max(Data.Var20)]);
  xmin = min([min(Data.Var19) min(Data.Var20)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  legend("Latif (\epsilon = "+epsilon{8}+")",...
         "Latif (\epsilon = "+epsilon{9}+")",...
         "Latif (\epsilon = "+epsilon{10}+")")
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s_compare.pdf',figDirPath,solvent(1:length(solvent)-4));
  saveas(gcf,fileName)
figure
  solvent=solvents{11};
  plot(Data.Var25,Data.Var26,'bo','markers',12,'linewidth',2)
  hold on
  plot(Data.Var23,Data.Var24,'ro','markers',12,'linewidth',2)
  xmax = max([max(Data.Var21) max(Data.Var22)]);
  xmin = min([min(Data.Var21) min(Data.Var22)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  set(gca,'FontSize',16)
  foo = refline(1,0);
  set(foo,'Linewidth',2,'color','k');
  xlabel(['\Delta','G_{expt}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  ylabel(['\Delta','G_{calc}^{solv, ',solvent(1:length(solvent)-4),'}',' (kcal/mol)'])
  legend("Latif (\epsilon = "+epsilon{11}+")",...
         "Latif (\epsilon = "+epsilon{12}+")")
  set(legend,'location','northwest')
  fileName=sprintf('%s/%s_compare.pdf',figDirPath,solvent(1:length(solvent)-4));
  saveas(gcf,fileName)