function epsilon = epsilon_octanol(temp)

% data from 
% A. Ghanadzadeh Gilani and M. Ansari
% J. Chem. Thermodynamics, v. 66, pp. 161-169 (2013)
% table 3

temp_list = [298.2 303.2 308.2 313.2 318.2];
eps_list  = [9.81 9.38 8.94 8.52 8.11];

polycoeffs = polyfit(temp_list,eps_list,2);
if temp < min(temp_list)
  fprintf('Warning in epsilon_octanol: requested temp %f < T_min %f!\n',temp,min(temp_list));
end
if temp > max(temp_list)
  fprintf('Warning in epsilon_octanol: requested temp %f > T_max %f!\n',temp,max(temp_list));
end
epsilon = polyval(polycoeffs,temp);
