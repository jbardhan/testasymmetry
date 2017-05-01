close all; clear all
addpath('export_fig/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Set your toggles and define the solute and solvent lists %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploton = 1;

solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'}; 
        
common_solutes = {'ethanol','butanone','toluene','n_octane','nitromethane',...
            '14_dioxane','phenol','acetic_acid','methanol','propanoic_acid',...
            'propan_1_ol','ethyl_acetate','aniline','ethylamine','n_butyl_acetate',...
            'methyl_acetate','hexan_1_ol','pentan_1_ol','n_propyl_acetate','butan_1_ol',...
            'heptan_1_ol','methyl_pentanoate','p_cresol','methyl_propanoate',...
            'o_cresol','propanone','pyridine'};
        
testset  = {'acetic_acid', 'ethanol', 'methanol', 'p_cresol',...
        'propanoic_acid', 'toluene', 'ethylamine', 'n_octane', 'pyridine',...
        'nitromethane', 'heptan_1_ol', 'n_butyl_acetate'};


% This routine will plot all of the calulated vs. experimental results if 
% plot is set to 1.  Futhermore, it creates structures that contain all of
% the results for calculated solvation free energies as a function of each 
% solvent and again as a function of solute
for i = 1:length(solvents)
           
    fid = fopen(['mnsol/',lower(solvents{i}),'.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};

    [~,m] = ismember(testset,mol_list);
    [~,n] = ismember(common_solutes,mol_list);
    
    allCases = 1:length(mol_list);
    allCases(m) = 0;
    temp = load(['Run',solvents{i}]);
    solute_errors(i,:) = temp.errfinal(n);
    sorted_errors = sort(abs(temp.errfinal));
    results(i) = struct('Solvent', solvents{i},'Num_Solutes',round(length(mol_list),3),'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE),...
                        'Mean_Abs_error',mean(abs(temp.errfinal)),'RMS_Training',...
                        rms(temp.errfinal(m)),'Mean_Abs_error_Training',...
                        mean(abs(temp.errfinal(m))),'RMS_Cons',...
                        rms(temp.errfinal(n)),'Mean_Abs_error_Cons',...
                        mean(abs(temp.errfinal(n))));
    if ploton
        figure
        plot(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)),'bo','markers',12)
        set(gca,'FontSize',15)
        hold on
        plot(temp.refE(m),temp.calcE(m),'rs','markers',12)
        plot([min(temp.refE)-2 max(temp.refE)+2] , [min(temp.calcE)-2 max(temp.calcE)+2],...
                    'k-','LineWidth',2)
        axis([min(temp.refE)-2 max(temp.refE)+2 min(temp.calcE)-2 max(temp.calcE)+2]);
        xlabel(['\Delta G_{expt}^{solv, ',solvents{i},'}'])
        ylabel(['\Delta G_{calc}^{solv, ',solvents{i},'}'])
        legend('Predictions','Training Set','Experiment','Location','southeast')
      filename = sprintf('Output/Figures/DeltaG-%s.PDF',solvents{i});
      export_fig(filename,'-painters','-transparent');

    end
    [~,m1] = ismember(sorted_errors(end),abs(temp.errfinal));
    [~,m2] = ismember(sorted_errors(end-1),abs(temp.errfinal));
    max_err(i) = struct('Solvent', solvents{i},'First_Max',...
        m1, 'Second_Max',m2);
    clear temp 
end


solute_errors = transpose(solute_errors);

% Plot a histogram of errors
if ploton
    figure 
    nbins = 25;
    for i = 1:length(results)
        subplot(3,3,i);
        histogram(results(i).errfinal,nbins)
        title(['Errors for ',results(i).Solvent])
        xlabel('Error')
        ylabel('Number of Occurances') 
    end
    filename = sprintf('Output/Figures/HistogramOfErrors.pdf');
    print(gcf, '-dpdf', filename); 
%     export_fig(filename,'-painters','-transparent');
end

% Create a structure conataining the errors associated with each solute
for i = 1:length(common_solutes)
    solute_struct(i) = struct('Solute_Name',common_solutes(i),'RMS',rms(solute_errors(i,:)),...
        'Mean_Abs_error',mean(abs(solute_errors(i,:))));
end

for i = 1:length(solvents)
    for j = i:length(solvents)
        rmsErr = calcRMSTrans(solvents{i},solvents{j});
        for k = 1:length(rmsErr)
            rmsdGTransErrorArray(i,j,k) = rmsErr(k);
        end
    end
    
    for j = i+1:length(solvents)
        [~,meanAbsError] = calcRMSTrans(solvents{i},solvents{j});
        for k = 1:length(meanAbsError)
            rmsdGTransErrorArray(j,i,k) = meanAbsError(k);
        end
    end
end

writeDat('Output/Tables/SolventErrors.tex',results,solvents);
writeDat('Output/Tables/SoluteErrors.tex',solute_struct,solvents);
writeDat('Output/Tables/TransferRMSErrors.tex',rmsdGTransErrorArray,solvents);


Outliers = readErr(max_err);
save('Solvent Outliers','Outliers');
