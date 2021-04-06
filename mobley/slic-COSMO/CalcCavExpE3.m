function cavE = CalcCavExpE3(solute_vdw_a,solute_vdw_v,total_atom_vol,cav_coeff1,cav_coeff2)
sphericity = pi^(1/3)*((6*solute_vdw_v)^(2/3))/solute_vdw_a;
shape_factor = (1/sphericity)^cav_coeff1;
packing_ratio = ((total_atom_vol-solute_vdw_v)/total_atom_vol)^cav_coeff2;

cavE = (2*shape_factor-1)*packing_ratio*(4-3*packing_ratio)/((1-packing_ratio)^2)-...
       (2*shape_factor-2)*log((1-0.5*packing_ratio)/(1-packing_ratio)^3);
