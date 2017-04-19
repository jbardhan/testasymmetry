close all; clear all
ploton = 0;
solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'}; 
        
testset  = {'acetic_acid', 'ethanol', 'methanol', 'p_cresol',...
        'propanoic_acid', 'toluene', 'ethylamine', 'n_octane', 'pyridine',...
        'nitromethane', 'heptan_1_ol', 'n_butyl_acetate'};
addpath('export_fig/')

         
for i = 1:length(solvents)
           
    fid = fopen(['mnsol/',lower(solvents{i}),'.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};

    [~,m] = ismember(testset,mol_list);
    temp = load(['Run',solvents{i}]);
    sorted_errors = sort(abs(temp.errfinal));
    results(i) = struct('Solvent', solvents{i},'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE));
    if ploton
        figure
        plot(temp.refE,temp.calcE,'bx','markers',12)
        set(gca,'FontSize',15)
        hold on
        plot(temp.refE(m),temp.calcE(m),'rx','markers',12)
        plot([min(temp.refE) max(temp.refE)] , [min(temp.calcE) max(temp.calcE)],...
                'k-')
        xlabel('Experimental Solvation Free Energies')
        ylabel('Calculated Solvation Free Energies')
        title(['Calculated vs. Experiment Solvation Free Energies for ',...
                solvents{i}])
      filename = sprintf('Output/DeltaG-%s.PDF',solvents{i});
      export_fig(filename,'-painters','-transparent');

    end
    [~,m1] = ismember(sorted_errors(end),abs(temp.errfinal));
    [~,m2] = ismember(sorted_errors(end-1),abs(temp.errfinal));
    max_err(i) = struct('Solvent', solvents{i},'First_Max',...
        m1, 'Second_Max',m2);
    clear temp 
end

Outliers = readErr(max_err);
save('Solvent Outliers','Outliers');
