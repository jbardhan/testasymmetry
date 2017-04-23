function [dG_list,error_code] = determinePaluchSolvationFreeEnergy(filename,mol_list, ...
						  temp)
dG_list = [];
error_code = 0;
[mol,t1,g1,t2,g2,t3,g3,ds,junk]=textread('reference-data/mobley_paluch_octanol.csv','%s%f%f%f%f%f%f%f%f','delimiter',',');

T = [t1 t2 t3];
G = [g1 g2 g3];

for i=1:length(mol_list)
  found_solute = 0;
  for j=1:length(mol)
    if strcmp(mol_list{i},mol{j})
      found_solute = j;
    end
  end
  
  if found_solute ==0 
    error_code = 1;
    return;
  end
  
  polycoeffs = polyfit(T(found_solute,:),G(found_solute,:),2);
  newDeltaG = polyval(polycoeffs, temp);
  dG_list(i) = newDeltaG;
end
