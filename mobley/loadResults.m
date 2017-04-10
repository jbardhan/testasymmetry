close all; clear all
solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'};
        
for i = 1:length(solvents)
    temp = load(['Run',solvents{i}]);
    results(i) = struct('Solvent', solvents{i},'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE));
    figure
    plot(temp.calcE,temp.refE,'bx')
    hold on
    plot([min(temp.calcE) max(temp.calcE)] , [min(temp.refE) max(temp.refE)],...
            'k-')
    xlabel('Calculated Solvation Free Energies')
    ylabel('Experimental Solvation Free Energies')
    title(['Experiment vs. Calculated Solvation Free Energies for ',...
            solvents{i}])
    clear temp 
end

