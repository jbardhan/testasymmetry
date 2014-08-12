doComputation = 1;
if doComputation
global quadrule_order quadrule_x quadrule_w;
quadrule_order = 10;
[quadrule_x,quadrule_w] = setupquad(quadrule_order);

loadconstants;
global E_omega E_sigma kappa;
E_inf   = 1.8;  % see note on AH thesis page 57 (65 of PDF)
E_vac   = 1;
E_sigma = 80;

lambda_vec  = [2 4 6 8 10 12 15 20 24];
E_omega_vec = [1 2 4 6 8 10 12 16 20];

prot_radius    = 13.4; % all following Gong, Hocky, and Freed (PNAS 2008)
ion_radius     = 1.4;
ion_separation = 2*ion_radius;  % simplest possible model
distance = [0 (prot_radius - ion_radius)];

protFilename='../geometry/protsphere/protsphere_0.25.srf';
[meshBase, rootDir] = readsrf(protFilename);
meshData = readmesh(meshBase,1);
[centroids,normals,areas] = genmeshcolloc(meshData);
[A] = collocation_mesh(meshData,centroids,normals,areas);
singleIonPQR = '../geometry/protsphere/singleIon.pqr';
singleIonSRF = '../geometry/protsphere/singleIon.srf';
pqrData = readpqr(singleIonPQR);

for lambda_index=1:length(lambda_vec)
  lambda = lambda_vec(lambda_index);
  bigLambda = lambda * sqrt(E_inf/E_sigma);
  kappa = 1/bigLambda;
  [A_Y] = collocation_mesh_Yukawa(meshData,centroids,normals,areas,kappa);

  for E_omega_index=1:length(E_omega_vec)
	 E_omega = E_omega_vec(E_omega_index);

	 % cost of burying a single ion again
	 E_nl_part0_singleIon(lambda_index,E_omega_index) = solvEnergyNonlocal(singleIonPQR,singleIonSRF,E_0,E_inf,E_omega,E_sigma,kappa);
	 E_nl_part0_singleIon_vac(lambda_index,E_omega_index) = solvEnergyLocal(singleIonPQR,singleIonSRF,E_0,E_inf,E_omega,E_vac);
	 E_nl_part0_singleIon_solv(lambda_index,E_omega_index) = E_nl_part0_singleIon(lambda_index,E_omega_index) ...
		  - E_nl_part0_singleIon_vac(lambda_index,E_omega_index);

	 for i=1:length(distance)
		pqrData.xyz(1,3) = distance(i);
		[B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrData);
		
		[Anl,Bnl,Cnl]=generate_nonlocal_matrices(A,B,C,A_Y,areas,E_0, ...
															  E_inf,E_omega,E_sigma);
		xnl = gmres(Anl,Bnl*pqrData.q,[],1e-6,100);
		Enl_twosweep(i,lambda_index,E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
			 (pqrData.q' * Cnl * xnl); 
		
		[Al,Bl,Cl]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega,E_sigma);
		xl = gmres(Al,Bl*pqrData.q,[],1e-6,100);
		El_twosweep(i,lambda_index,E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
			 (pqrData.q' * Cl * xl); 
		
		[Aref,Bref,Cref]=generate_local_matrices(A,B,C,areas,...
															  E_0,E_inf,E_omega, ...
															  E_vac);
		xref = gmres(Aref,Bref*pqrData.q,[],1e-6,100);
		Eref_twosweep(i,lambda_index,E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
			 (pqrData.q' * Cref * xref);
		
	 end % distances
  end % E_omega loop
end % lambda loop

save table1
end % docomputation
