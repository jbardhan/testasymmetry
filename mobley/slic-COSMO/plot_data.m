function plot_data(x_data,y_data,data,solvent,thermo_var,plot_tr,hb_compounds)

  Home = getenv('HOME');
  legend_font = 'Helvetica';
  label_font = 'Helvetica';
  legend_font_size = 14;
  label_font_size = 18;
  %readtable('../hb.csv')

  switch nargin
      case 3
          solvent = 'water';
          thermo_var = 'dGexpt';
          plot_tr = 1;
          hb_compounds = 0;
      case 4
          thermo_var = 'dGexpt';
          plot_tr = 1;
          hb_compounds = 0;
      case 5
          plot_tr = 1;
          hb_compounds = 0;
      case 6
          hb_compounds = 0;
          
  end
  repo_path = sprintf('%s/repos/testasymmetry/mobley/slic-COSMO',Home);
  figures_path = sprintf('%s/figures',repo_path);
  fig = figure('Renderer', 'painters', 'Position', [12 8 800 600]);
  get_index = @(x,n) x(n);
  thermo_var_container = containers.Map; % works like a dictionary in python
  thermo_var_container('dGexpt')={"$\Delta\,G_{\,\mathrm{expt}}^{\,\mathrm{solv}}$",...
  "$\Delta\,G_{\,\mathrm{calc,\,SLIC}}^{\,\mathrm{solv}} $ "," (kcal/mol)"};
  thermo_var_container('dGexptSASA')={"$\Delta\,G_{\,\mathrm{expt}}^{\,\mathrm{solv}}$",...
  "$\Delta\,G_{\,\mathrm{calc,\,SLIC-SASA}}^{\,\mathrm{solv}} $ "," (kcal/mol)"};
  thermo_var_container('dGMD')={"$\Delta\,G_{\,\mathrm{MD}}^{\,\mathrm{solv}}$",...
  "$\Delta\,G_{\,\mathrm{calc,\,SLIC}}^{\,\mathrm{solv}} $ "," (kcal/mol)"};
  thermo_var_container('dGcav')={"$\Delta\,G_{\,\mathrm{MD,\,cav}}^{\,\mathrm{solv}}$",...
  "$\Delta\,G_{\,\mathrm{calc}}^{\,\mathrm{solv,\,cav}} $ "," (kcal/mol)"};
  thermo_var_container('dGdisp')={"$\Delta\,G_{\,\mathrm{MD,\,disp}}^{\,\mathrm{solv}}$",...
  "$\Delta\,G_{\,\mathrm{calc}}^{\,\mathrm{solv,\,disp}} $ "," (kcal/mol)"};
  thermo_var_container('dGnp')={"$\Delta\,G_{\,\mathrm{MD}}^{\,\mathrm{solv,\,np}}$",...
  "$\Delta\,G_{\,\mathrm{calc}}^{\,\mathrm{solv,\,np}} $ "," (kcal/mol)"};
  thermo_var_container('dGnpSASA')={"$\Delta\,G_{\,\mathrm{MD}}^{\,\mathrm{solv,\,np}}$",...
  "$\Delta\,G_{\,\mathrm{calc}}^{\,\mathrm{solv,\,np-SASA}} $ "," (kcal/mol)"};
  thermo_var_container('dGes')={"$\Delta\,G_{\,\mathrm{MD}}^{\,\mathrm{solv,\,es}}$",...
  "$\Delta\,G_{\,\mathrm{calc,\,SLIC}}^{\,\mathrm{solv,\,es}} $ "," (kcal/mol)"};
  thermo_var_container('TdSMD')={"$-T\,\Delta\,S_{\,\mathrm{MD}}^{\,\mathrm{solv}}$",...
  "$-T\,\Delta\,S_{\,\mathrm{calc,\,SLIC}}^{\,\mathrm{solv}} $ "," (kcal/mol)"};
  thermo_var_container('dCp')={"$\Delta\,C_{p\,\mathrm{expt}}^{\,\mathrm{solv}}$",...
  "$\Delta\,C_{p\,\mathrm{calc,\,SLIC}}^{\,\mathrm{solv}}$ "," (kcal/mol/K)"};
  plot(x_data,y_data,'bo','markers',12,'linewidth',2);
  hold on
  xmax = max([max(x_data) max(y_data)]);
  xmin = min([min(x_data) min(y_data)]);
  axis([xmin-1 xmax+1 xmin-1 xmax+1]);
  xlab = get_index(thermo_var_container(thermo_var),1);
  ylab = get_index(thermo_var_container(thermo_var),2);
  var_unit = get_index(thermo_var_container(thermo_var),3);
  var_unit = string(var_unit);
  rmse = rms(x_data-y_data);
  rmsErr = sprintf('%4.2f',rmse);
  if ~plot_tr
     
     if hb_compounds
          hb_index = find(data.hbdata);
          rmse_hb = rms(x_data(hb_index)-y_data(hb_index));
          rmsErrhb = sprintf('%4.2f',rmse_hb);
          hold on
          plot(x_data(hb_index),y_data(hb_index),'m+','markers',10,'linewidth',3)
          set(gca,'FontSize',16)
          foo = refline(1,0);
          set(foo,'Linewidth',2,'color','k');
          hold off
          legend({"RMSE Test-set = "+rmsErr+" "+var_unit,...
              "RMSE HB-compounds = "+rmsErrhb+" "+var_unit},...
          'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
          
          set(legend,'location','northwest')
     else
         set(gca,'FontSize',16)
         foo = refline(1,0);
         set(foo,'Linewidth',2,'color','k');
         hold off
         leg_1 = sprintf('$\\mathrm{RMSE}  =  %s \\ \\mathrm{%s}$',...
         rmsErr,var_unit);

         legend({leg_1} ,'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
     end
     xlabel(string(xlab) + solvent + var_unit,...
           "Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
     ylabel(string(ylab) + solvent + var_unit,...
           "Interpreter","latex",...
           'FontName',label_font,'FontSize',label_font_size)
     set(legend,'location','northwest')
     fileName=sprintf('%s/%s.pdf',figures_path,thermo_var);
     orient(fig,'landscape')
     print(fig,fileName,'-dpdf')    
  else
      [~, index] = ismember(data.training_set,data.mol_list);
      rmse_tr = rms(x_data(index)-y_data(index));
      rmsErrTr = sprintf('%4.2f',rmse_tr);
      plot(x_data(index),y_data(index),'r*','markers',10,'linewidth',3)
      
      if hb_compounds
          hb_index = find(data.hbdata);
          rmse_hb = rms(x_data(hb_index)-y_data(hb_index));
          rmsErrhb = sprintf('%4.2f',rmse_hb);
          hold on
          plot(x_data(hb_index),y_data(hb_index),'m+','markers',10,'linewidth',3)
          set(gca,'FontSize',16)
          foo = refline(1,0);
          set(foo,'Linewidth',2,'color','k');
          hold off
          legend({"RMSE Test-set = "+rmsErr+" "+var_unit,...
              "RMSE Training-set = "+rmsErrTr+" "+var_unit,...
              "RMSE HB-compounds = "+rmsErrhb+" "+var_unit},...
          'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
          
          set(legend,'location','northwest')
      
      else
          set(gca,'FontSize',16)
          foo = refline(1,0);
          set(foo,'Linewidth',2,'color','k');
          hold off
          legend({"RMSE Test-set = "+rmsErr+" "+var_unit,...
              "RMSE Training-set = "+rmsErrTr+" "+var_unit},...
          'FontName',legend_font,'FontSize',legend_font_size,"Interpreter","latex")
          set(legend,'location','northwest')   

      end
          
  end
  xlabel(string(xlab) + solvent + var_unit,...
                   "Interpreter","latex",...
                   'FontName',label_font,'FontSize',label_font_size)
  ylabel(string(ylab) + solvent + var_unit,...
                   "Interpreter","latex",...
                   'FontName',label_font,'FontSize',label_font_size)

  fileName=sprintf('%s/%s.pdf',figures_path,thermo_var);
  orient(fig,'landscape')
  print(fig,fileName,'-dpdf')
