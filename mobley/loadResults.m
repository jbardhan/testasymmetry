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
solvents = solvents;
        
common_solutes = {'ethanol','butanone','toluene','n_octane','nitromethane',...
            '14_dioxane','phenol','acetic_acid','methanol','propanoic_acid',...
            'propan_1_ol','ethyl_acetate','aniline','ethylamine','n_butyl_acetate',...
            'methyl_acetate','hexan_1_ol','pentan_1_ol','n_propyl_acetate','butan_1_ol',...
            'heptan_1_ol','methyl_pentanoate','p_cresol','methyl_propanoate',...
            'o_cresol','propanone','pyridine'};
        

% This routine will plot all of the calulated vs. experimental results if 
% plot is set to 1.  Futhermore, it creates structures that contain all of
% the results for calculated solvation free energies as a function of each 
% solvent and again as a function of solute
k = 1;
for i = 1:length(solvents)
           
    fid = fopen(['mnsol/',lower(solvents{i}),'.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};

    [~,n] = ismember(common_solutes,mol_list);
    
    temp = load(['Run',solvents{i}]);
    for j = 1:length(temp.calcE(1,:))
        testset = temp.testSets(:,j);
        [~,m] = ismember(testset,mol_list);
        allCases = 1:length(mol_list);
        allCases(m) = 0;
        solute_errors(j,:,i) = temp.errfinal(n,j);
        sorted_errors = sort(abs(temp.errfinal(:,j)));
        results(k) = struct('Solvent', solvents{i},'Num_Solutes',round(length(mol_list),3),'CalcE',temp.calcE(:,j),...
                            'errfinal',temp.errfinal(:,j),'es',temp.es(:,j),'np',...
                            temp.np(:,j),'refE',temp.refE(:,j),'RMS',rms(temp.calcE(:,j)-temp.refE(:,j)),...
                            'Mean_Abs_error',mean(abs(temp.errfinal(:,j))),'RMS_Training',...
                            rms(temp.errfinal(m,j)),'Mean_Abs_error_Training',...
                            mean(abs(temp.errfinal(m,j))),'RMS_Cons',...
                            rms(temp.errfinal(n,j)),'Mean_Abs_error_Cons',...
                            mean(abs(temp.errfinal(n,j))),'Test_Set',...
                            {testset});
        
        if ploton
            figure
            refDat = temp.refE(:,j);
            calcDat = temp.calcE(:,j);
            plot(refDat(allCases(allCases~=0)),calcDat(allCases(allCases~=0)),'bo','markers',12)
            set(gca,'FontSize',15)
            hold on
            plot(temp.refE(m),temp.calcE(m),'rs','markers',12)
            minax = round(min(min(refDat,calcDat)));
            maxax = round(max(max(refDat,calcDat)));
            axis([minax-2 maxax+2 minax-2 maxax+2]);
            foo = refline(1,0);
            set(foo,'Linewidth',2,'color','k');
            xlabel(['\Delta G_{expt}^{solv, ',solvents{i},'}'])
            ylabel(['\Delta G_{calc}^{solv, ',solvents{i},'}'])
            legend('Predictions','Training Set','Location','southeast')
          filename = sprintf('Output/Figures/DeltaG-%s%s.PDF',solvents{i},string(j));
          export_fig(filename,'-painters','-transparent');
        end
        [~,m1] = ismember(sorted_errors(end),abs(temp.errfinal(:,j)));
        [~,m2] = ismember(sorted_errors(end-1),abs(temp.errfinal(:,j)));
        max_err(k) = struct('Solvent', solvents{i},'First_Max',...
            m1, 'Second_Max',m2);
        k = k + 1;
    end
end


solute_errors = permute(solute_errors,[2 1 3]);

% Plot a histogram of errors
if ploton 
    nbins = 25;
    for i = 1:length(solvents)
        figure
        m = ismember({results.Solvent},solvents{i});
        [~,n] = ismember(solvents{i},{results.Solvent});
        for j = 0:length(m(m~=0))-1
            subplot(3,4,j+1);
            histogram(results(n(1)+j).errfinal,nbins)
            title(['Errors for ',results(n(1)+j).Solvent])
            xlabel('Error')
            ylabel('Number of Occurances') 
        end
        filename = sprintf('Output/Figures/HistogramOfErrors%s.PDF',solvents{i});
        print(gcf, '-dpdf', filename); 
%     export_fig(filename,'-painters','-transparent');
    end
    
end

% Create a structure conataining the errors associated with each solute
k = 1;
for i = 1:length(common_solutes)
    for j = 1:length(solvents)
        solute_struct(k) = struct('Solvnet',solvents{j},'Solute_Name',common_solutes(i),'RMS',rms(solute_errors(i,:,j)),...
            'Mean_Abs_error',mean(abs(solute_errors(i,:,j))));
        k = k + 1;
    end
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

for i = 1:length(rmsdGTransErrorArray(:,1,1))
    for j = 1:length(rmsdGTransErrorArray(1,:,1)) 
        stdevdGTransOverTestSets(i,j) = std(rmsdGTransErrorArray(i,j,:));
    end
end

writeDat('Output/Tables/SolventErrors.tex',results,solvents);
writeDat('Output/Tables/SoluteErrors.tex',solute_struct,solvents);
writeDat('Output/Tables/TransferRMSErrors.tex',rmsdGTransErrorArray,solvents);
writeDat('Output/Tables/sdevRMSdGtrans.tex',stdevdGTransOverTestSets,solvents);

Outliers = readErr(max_err);
save('Solvent Outliers','Outliers');
