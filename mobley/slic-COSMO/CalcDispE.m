function dispE = CalcDispE(solute_data,solute_atom_types,...
    solvent_data,solvent_atom_types,...
    disp_coeffs,q_s,Temp)

all_atom_types = {'br','c-sp','c-sp2','c-sp3','cl',...
                  'f','h','i','n-sp','n-sp3',...
                  'n-sp2','o-sp2','o-sp3-h',...
                  'o-sp2-n','o-sp3','p','s'};
solute_areas = solute_data(7:end);
unique_atom_types_solute = unique(string(solute_atom_types));
unique_atom_types_solute = unique_atom_types_solute(2:end);

vol_mix = solvent_data(4);
solvent_areas = solvent_data(7:end);
unique_atom_types_solvent = unique(string(solvent_atom_types));
unique_atom_types_solvent = unique_atom_types_solvent(2:end);

%dispSolventSolute = zeros(1,length(all_atom_types));
solute_types = cell2table(solute_atom_types);
solvent_types = cell2table(solvent_atom_types);

dispE=0;


for i=1:length(unique_atom_types_solvent)
    
    solvent_atom_type = unique_atom_types_solvent{i};
    [~,id_coeff_solvent] = ismember(solvent_atom_type,all_atom_types);
    coeff_solvent = disp_coeffs(id_coeff_solvent);
    [m, ~] = ismember(solvent_types.Variables,{solvent_atom_type});
    areas_solvent = solvent_areas(m);
    m_tau_solvent = sum(areas_solvent.^q_s);
    dispSouteSolvent=0;
    for j=1:length(unique_atom_types_solute)
        
        solute_atom_type = unique_atom_types_solute{j};
        [~,id_coeff_solute] = ismember(solute_atom_type,all_atom_types);
        coeff_solute = disp_coeffs(id_coeff_solute);
        [m, ~] = ismember(solute_types.Variables,{solute_atom_type});
        areas_solute = solute_areas(m);
        m_tau_solute = sum(areas_solute.^q_s);
        dispSouteSolvent = dispSouteSolvent + m_tau_solute*sqrt(coeff_solute);
        
    end
    
    dispE = dispE + m_tau_solvent*sqrt(coeff_solvent)*dispSouteSolvent;
    
end
kB = 0.001987;
dispE = kB*dispE/vol_mix; %in kcal/mol
