close all; clear all
ploton = 0;
solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'};
        
for i = 1:length(solvents)
    temp = load(['Run',solvents{i}]);
    sorted_errors = sort(abs(temp.errfinal));
    results(i) = struct('Solvent', solvents{i},'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE));
    if ploton
        figure
        plot(temp.calcE,temp.refE,'bx')
        hold on
        plot([min(temp.calcE) max(temp.calcE)] , [min(temp.refE) max(temp.refE)],...
                'k-')
        xlabel('Calculated Solvation Free Energies')
        ylabel('Experimental Solvation Free Energies')
        title(['Experiment vs. Calculated Solvation Free Energies for ',...
                solvents{i}])
    end
    max_err(i) = struct('Solvent', solvents{i},'Fist_Max',...
        sorted_errors(end), 'Second_Max',sorted_errors(end-1));
    clear temp 
end

save('largestErrors','max_err');

