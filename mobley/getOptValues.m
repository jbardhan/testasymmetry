clear all; clear global
Home = getenv('HOME');
global LoadData
addpath(sprintf('%s/repos/testasymmetry',Home));
chdir(sprintf('%s/repos/testasymmetry/mobley',Home));
solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'}; 
solvents = strcat('Opt',solvents);        
dirInfo = dir(pwd);
filesInDir = {dirInfo.name};

for i = 1:length(filesInDir)
    if dirInfo(i).isdir == 0
        curFile = strsplit(filesInDir{i},'_');
        if ismember(curFile{1},solvents)
           loadFile = strjoin(curFile,'_');
           solvent = strsplit(curFile{1},'Opt');
           addLoadData(solvent{2},loadFile);
        end
    end
end
addLoadData('Water','OptWater_25-Apr-2017 12:50:20.mat');
save('OptFiles','LoadData');
clear global; clear all