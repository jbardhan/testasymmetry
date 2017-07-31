close all; clear all;
solvents = {'Water','Methanol','Ethanol','Propanol','Butanol','Octanol'};
clc
for i = 1:length(solvents)
    ans = load(['Opt',solvents{i}]);
    rmsE=rms(ans.calc-ans.ref);
    formatSpec = 'testset RMS error for %s is %5.3f\n';
    fprintf(formatSpec,solvents{i},rmsE)
end

    fprintf('\n')

for i = 1:length(solvents)
    ans = load(['Run',solvents{i}]);
    rmsE=rms(ans.calcE-ans.refE);
    formatSpec = 'All neutral compounds RMS error for %s is %5.3f\n';
    fprintf(formatSpec,solvents{i},rmsE)
end
    
    fprintf('\n')

for i = 1:length(solvents)
    ans = load(['Run',solvents{i}]);
    meanabsE=meanabs(ans.calcE-ans.refE);
    formatSpec = 'All neutral compounds mean absolute error for %s is %5.3f\n';
    fprintf(formatSpec,solvents{i},meanabsE)
end
