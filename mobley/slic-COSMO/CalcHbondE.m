function hbondE = CalcHbondE(solute_hbond_data,solvent_hbond_data,hbond_coeffs)

solute_acceptor_grp = solute_hbond_data(1:end-5);
solvent_acceptor_grp = solvent_hbond_data(1:end-5);
solute_donor_grp = solute_hbond_data(end-4:end-2);
solvent_donor_grp = solvent_hbond_data(end-4:end-2);
solute_acceptor_ids = find(solute_acceptor_grp);
solvent_acceptor_ids = find(solvent_acceptor_grp);
solute_donor_ids = find(solute_donor_grp);
solvent_donor_ids = find(solvent_donor_grp);
solute_acceptor_coeffs = hbond_coeffs(solute_acceptor_ids);
solvent_acceptor_coeffs = hbond_coeffs(solvent_acceptor_ids);
hbond_coeffs_donor=hbond_coeffs(end-2:end);
solute_donor_coeffs = hbond_coeffs_donor(solute_donor_ids);
solvent_donor_coeffs = hbond_coeffs_donor(solvent_donor_ids);
solute_acceptor_array = solute_acceptor_grp(solute_acceptor_ids);
solvent_acceptor_array = solvent_acceptor_grp(solvent_acceptor_ids);
solute_donor_array = solute_donor_grp(solute_donor_ids);
solvent_donor_array = solvent_donor_grp(solvent_donor_ids);

hbondE = -sum(solvent_donor_array.*solvent_donor_coeffs)*sum(solute_acceptor_array.*solute_acceptor_coeffs) + ...
         -sum(solute_donor_array.*solute_donor_coeffs)*sum(solvent_acceptor_array.*solvent_acceptor_coeffs) ;

