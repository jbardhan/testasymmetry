function [Solutes] = readErr(structured)

solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'};

for i = 1:length(solvents)
    fid = fopen(['mnsol/',lower(solvents{i}),'.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};
    [~,n] = ismember(solvents{i},{structured.Solvent});
    Sol_One = mol_list{structured(n).First_Max};
    Sol_Two = mol_list{structured(n).Second_Max};
    Solutes(i) = struct('Solvent', solvents{i},...
        'First_Solute',Sol_One,'Second_Solute',Sol_Two);
end