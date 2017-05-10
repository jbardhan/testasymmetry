close all; clear all
addpath('export_fig/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Set your toggles and define the solute and solvent lists %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploton = 1;

solvents = {'Water','Octanol','Dichloroethane', 'Propanone', 'Dimethylsulfoxide',...
             'Propanol', 'Dimethylformamide', 'Ethanol', 'Methanol', 'Acetonitrile'}; 
solvents_1 = {'Octanol'};
       
%common_solutes = {'ethanol','butanone','n_octane','pyrene','cyclohexane',};
        
testset  = {'butanone','n_octane','ethanol','benzene','cyclohexane','pyrene'};
testset_ions = {'Li','Na','K','Cl','Br','I'};
test_ion = {'Rb','Cs'};


% This routine will plot all of the calulated vs. experimental results if 
% plot is set to 1.  Futhermore, it creates structures that contain all of
% the results for calculated solvation free energies as a function of each 
% solvent and again as a function of solute
for i = 1:length(solvents)
           
    fid = fopen(['mnsol/',lower(solvents{i}),'_ions.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};
    [~,m] = ismember(testset,mol_list);
    [~,n] = ismember(testset_ions,mol_list);
    
    allCases = 1:length(mol_list);
    if ~ismember(solvents{i},solvents_1)
        [~,p] = ismember(test_ion,mol_list);
        allCases(p) = 0;
    end 
    
    allCases(m) = 0;
    allCases(n) = 0;
    
    temp = load(['Run',solvents{i}]);
%     solute_errors(i,:) = temp.errfinal(n);
%    sorted_errors = sort(abs(temp.errfinal));
    results(i) = struct('Solvent', solvents{i},'Num_Solutes',round(length(mol_list),3),'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE),...
                        'Mean_Abs_error',mean(abs(temp.errfinal)));
    if ploton
        figure
        plot(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)),'bo','markers',12)
        %plot(temp.refE,temp.calcE,'bo','markers',12)
        set(gca,'FontSize',15)
        hold on
        plot(temp.refE(m),temp.calcE(m),'rs','markers',12)
        minax = round(min(min(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)))));
        maxax = round(max(max(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)))));
        axis([minax-2 maxax+2 minax-2 maxax+2]);
        foo = refline(1,0);
        set(foo,'Linewidth',2,'color','k');
        xlabel(['\Delta G_{expt}^{solv, ',solvents{i},'}'])
        ylabel(['\Delta G_{calc}^{solv, ',solvents{i},'}'])
        legend('Predictions','Training Set','Location','southeast')
      filename = sprintf('Output/Figures/DeltaG-%s.PDF',solvents{i});
      export_fig(filename,'-painters','-transparent');

    end
    
end
clear all
addpath('export_fig/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Set your toggles and define the solute and solvent lists %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploton = 1;

solvents = {'Water','Octanol','Dichloroethane', 'Propanone', 'Dimethylsulfoxide',...
             'Propanol', 'Dimethylformamide', 'Ethanol', 'Methanol', 'Acetonitrile'}; 
solvents_1 = {'Octanol'};
       
%common_solutes = {'ethanol','butanone','n_octane','pyrene','cyclohexane',};
        
testset_ions = {'Li','Na','K','Cl','Br','I'};
test_ion = {'Rb','Cs'};
for i = 1:length(solvents)
           
    fid = fopen(['mnsol/',lower(solvents{i}),'_ions.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};
%     [~,m] = ismember(testset,mol_list);
    [~,n] = ismember(testset_ions,mol_list);
    
    allCases = 1:length(mol_list);
    if ~ismember(solvents{i},solvents_1)
        [~,p] = ismember(test_ion,mol_list);
    end 
    
%     allCases(~m) = 0;
%     allCases(~n) = 0;
    
    temp = load(['Run',solvents{i}]);
%     solute_errors(i,:) = temp.errfinal(n);
%    sorted_errors = sort(abs(temp.errfinal));
    results(i) = struct('Solvent', solvents{i},'Num_Solutes',round(length(mol_list),3),'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE),...
                        'Mean_Abs_error',mean(abs(temp.errfinal)));
    if ploton
        figure
        if ~ismember(solvents{i},solvents_1)
            [~,p] = ismember(test_ion,mol_list);
            plot(temp.refE(allCases(p)),temp.calcE(allCases(p)),'bo','markers',12)
            hold on
        else
            plot(800,800,'bo','markers',12)
            hold on
        end
        plot(temp.refE(allCases(n)),temp.calcE(allCases(n)),'rs','markers',12)
        %plot(temp.refE,temp.calcE,'bo','markers',12)
        set(gca,'FontSize',15)
        minax = round(min(min(temp.refE(allCases(n)),temp.calcE(allCases(n)))));
        maxax = round(max(max(temp.refE(allCases(n)),temp.calcE(allCases(n)))));
        axis([minax-2 maxax+2 minax-2 maxax+2]);
        foo = refline(1,0);
        set(foo,'Linewidth',2,'color','k');
        xlabel(['\Delta G_{expt}^{solv, ',solvents{i},'}'])
        ylabel(['\Delta G_{calc}^{solv, ',solvents{i},'}'])
        legend('Predictions','Training Set','Location','southeast')
      filename = sprintf('Output/Figures/DeltaG-%s.PDF',solvents{i});
      export_fig(filename,'-painters','-transparent');

    end
    %[~,m1] = ismember(sorted_errors(end),abs(temp.errfinal));
    %[~,m2] = ismember(sorted_errors(end-1),abs(temp.errfinal));
    %max_err(i) = struct('Solvent', solvents{i},'First_Max',...
    %    m1, 'Second_Max',m2);
    clear temp 
end


% solute_errors = transpose(solute_errors);

% Plot a histogram of errors
% if ploton
%     figure 
%     nbins = 25;
%     for i = 1:length(results)
%         subplot(3,3,i);
%         histogram(results(i).errfinal,nbins)
%         title(['Errors for ',results(i).Solvent])
%         xlabel('Error')
%         ylabel('Number of Occurances') 
%     end
%     filename = sprintf('Output/Figures/HistogramOfErrors.pdf');
%     print(gcf, '-dpdf', filename); 
% %     export_fig(filename,'-painters','-transparent');
% end
% 
% % Create a structure conataining the errors associated with each solute
% for i = 1:length(common_solutes)
%     solute_struct(i) = struct('Solute_Name',common_solutes(i),'RMS',rms(solute_errors(i,:)),...
%         'Mean_Abs_error',mean(abs(solute_errors(i,:))));
% end
% 
% for i = 1:length(solvents)
%     for j = i:length(solvents)
%         rmsErr = calcRMSTrans(solvents{i},solvents{j});
%         for k = 1:length(rmsErr)
%             rmsdGTransErrorArray(i,j,k) = rmsErr(k);
%         end
%     end
%     
%     for j = i+1:length(solvents)
%         [~,meanAbsError] = calcRMSTrans(solvents{i},solvents{j});
%         for k = 1:length(meanAbsError)
%             rmsdGTransErrorArray(j,i,k) = meanAbsError(k);
%         end
%     end
% end

% writeDat('Output/Tables/SolventErrors.tex',results,solvents);
% writeDat('Output/Tables/SoluteErrors.tex',solute_struct,solvents);
% writeDat('Output/Tables/TransferRMSErrors.tex',rmsdGTransErrorArray,solvents);


Outliers = readErr(max_err);
save('Solvent Outliers','Outliers');
