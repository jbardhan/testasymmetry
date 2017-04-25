Home = getenv('HOME');
addpath(sprintf('%s/repos/testasymmetry',Home));
chdir(sprintf('%s/repos/testasymmetry/mobley',Home));
solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'}; 
solvents = strcat('Opt',solvents);        
dirInfo = dir(pwd);
filesInDir = {dirInfo.name};

for i = 1:length(filesInDir)
    curFile = strsplit(filesInDir{i}
   if  
end