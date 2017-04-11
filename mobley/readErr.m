function [Solutes] = readErr(struct)

solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'};

for i = 1:length(solvents)
    fid = fopen(['mnsol/',lower(solvents{i}),'.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};
    [~,n] = ismember(solvents{i},{struct.Solvent});
    Solutes(i) = struct('Solvent', solvents{i}, 'Solute_One',...
        string(mol_list{struct(n).First_Max}),'Solute_Two',...
        string(mol_list{struct(n).Second_Max}));
end