for i=1:length(E_omega_vec)
  E_omega = E_omega_vec(i);
  for j=1:length(distance)
	 protonatedCurData = centerAtoms(protonatedPqrData);
	 protonatedCurData = translateAtoms(protonatedCurData, [0 0 ...
						  distance(j)]);
	 deprotonatedCurData = centerAtoms(deprotonatedPqrData);
	 deprotonatedCurData = translateAtoms(deprotonatedCurData, [0 0 ...
						  distance(j)]);
	 
	 E_prot_Coul(j,i) = getCoulomb(protonatedCurData, E_omega);
	 E_deprot_Coul(j,i) = getCoulomb(deprotonatedCurData, E_omega);
  end
end


Enlprot_solv = Enlprot - Erefprot;
Enldeprot_solv = Enldeprot - Erefdeprot;
Elprot_solv = Elprot - Erefprot;
Eldeprot_solv = Eldeprot - Erefdeprot;
dE_vac = Erefprot - Erefdeprot;
dE_Coul = E_prot_Coul - E_deprot_Coul;

dE_prot_total_nl = -Enldeprot_solv + (dE_vac + dE_Coul) + ...
	 Enlprot_solv;
dE_prot_total_l = -Eldeprot_solv + (dE_vac + dE_Coul) + ...
	 Elprot_solv;

Enlres_prot_solv = estimateDetails.prot_nl-estimateDetails.prot_ref;
Elres_prot_solv  = estimateDetails.prot_l-estimateDetails.prot_ref;

Enlres_deprot_solv = estimateDetails.deprot_nl-estimateDetails.deprot_ref;
Elres_deprot_solv  = estimateDetails.deprot_l-estimateDetails.deprot_ref;

dE_res_vac = estimateDetails.prot_ref  - estimateDetails.deprot_ref;
dE_res_Coul = estimateDetails.dE_Coul;
dE_res_prot_total_nl = -Enlres_deprot_solv + (dE_res_vac+dE_res_Coul) ...
	 + Enlres_prot_solv;
dE_res_prot_total_l  = -Elres_deprot_solv + (dE_res_vac+dE_res_Coul) ...
	 + Elres_prot_solv;