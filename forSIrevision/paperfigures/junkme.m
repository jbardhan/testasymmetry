doComputation = 1;  % leave this alone for now =)

% these files are what you 
protFilename='../geometry/protsphere/protsphere_0.15.srf';
prot_radius    = 13.4; % all following Gong, Hocky, and Freed (PNAS 2008)
%  for protsphere2, change the above to prot_radius=20.0;
deprotonatedPQR = '../geometry/pka_shift_1/glu_capped_charmm22.pqr';
protonatedPQR   = '../geometry/pka_shift_1/glu_capped_prot_charmm22.pqr';
gluSRF = '../geometry/pka_shift_1/glu_2.srf';
% for the lysine, use pka_shift_2 pqr files and srf
% end easily changed parameters 

if doComputation
global quadrule_order quadrule_x quadrule_w;
quadrule_order = 10;
[quadrule_x,quadrule_w] = setupquad(quadrule_order);

loadconstants;
global E_omega E_sigma kappa;
E_inf   = 1.8;  % see note on AH thesis page 57 (65 of PDF)
E_vac   = 1;
E_sigma = 80;
lambda_vec  = [15];
%bigLambda = lambda * sqrt(E_inf/E_sigma);
%kappa = 1/bigLambda;

E_omega_vec =[2];



[meshBase, rootDir] = readsrf(protFilename);
meshData = readmesh(meshBase,1);
[centroids,normals,areas] = genmeshcolloc(meshData);
[A] = collocation_mesh(meshData,centroids,normals,areas);

protonatedPqrData = readpqr(protonatedPQR);
deprotonatedPqrData = readpqr(deprotonatedPQR);
residueRadius = getPointSetRadius(protonatedPqrData);


distance = [0 (prot_radius - (residueRadius+1))];

[pKa_nl, pKa_l, estimateDetails]= estimate_pKa(protonatedPQR, ...
															  deprotonatedPQR, ...
															  gluSRF,E_omega_vec, ...
															  E_sigma, lambda_vec,E_vac);
for lambda_index=1:length(lambda_vec)
  lambda = lambda_vec(lambda_index);
  bigLambda = lambda * sqrt(E_inf/E_sigma);
  kappa = 1/bigLambda;
  [A_Y] = collocation_mesh_Yukawa(meshData,centroids,normals,areas, ...
											 kappa);
  
  for E_omega_index=1:length(E_omega_vec)
  E_omega = E_omega_vec(E_omega_index);

  for i=1:length(distance)
	 protonatedCurData = centerAtoms(protonatedPqrData);
	 protonatedCurData = translateAtoms(protonatedCurData, [0 0 ...
						  distance(i)]);
	 deprotonatedCurData = centerAtoms(deprotonatedPqrData);
	 deprotonatedCurData = translateAtoms(deprotonatedCurData, [0 0 ...
						  distance(i)]);
	 
	 E_prot_Coul(i,E_omega_index) = getCoulomb(protonatedCurData, E_omega);
	 E_deprot_Coul(i,E_omega_index) = getCoulomb(deprotonatedCurData, E_omega);

	 [Bprot,Cprot] = chargeCollocation_mesh(meshData,centroids,normals,areas,protonatedCurData);
	 [Bdeprot,Cdeprot] = chargeCollocation_mesh(meshData,centroids,normals,areas,deprotonatedCurData);
	 
	 [Anlprot,Bnlprot,Cnlprot]=generate_nonlocal_matrices(A,Bprot,Cprot,A_Y, ...
																  areas,E_0, ...
																  E_inf,E_omega,E_sigma);
	 xnlprot = gmres(Anlprot,Bnlprot*protonatedCurData.q,[],1e-6,300);
	 Enlprot(i,lambda_index, E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
		  (protonatedCurData.q' * Cnlprot * xnlprot); 
	 
	 [Anldeprot,Bnldeprot,Cnldeprot]=generate_nonlocal_matrices(A,Bdeprot,Cdeprot,A_Y, ...
																  areas,E_0, ...
																  E_inf,E_omega,E_sigma);
	 xnldeprot = gmres(Anldeprot,Bnldeprot*deprotonatedCurData.q,[],1e-6,300);
	 Enldeprot(i,lambda_index,E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
		  (deprotonatedCurData.q' * Cnldeprot * xnldeprot); 
	 
	 [Alprot,Blprot,Clprot]=generate_local_matrices(A,Bprot,Cprot,areas,E_0,E_inf,E_omega,E_sigma);
	 xlprot = gmres(Alprot,Blprot*protonatedCurData.q,[],1e-6,300);
	 Elprot(i,E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
		  (protonatedCurData.q' * Clprot * xlprot); 
	 
	 [Aldeprot,Bldeprot,Cldeprot]=generate_local_matrices(A,Bdeprot, ...
																  Cdeprot,areas,E_0,E_inf,E_omega,E_sigma);
	 
	 xldeprot = gmres(Aldeprot,Bldeprot*deprotonatedCurData.q,[],1e-6,300);
	 Eldeprot(i,E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
		  (deprotonatedCurData.q' * Cldeprot * xldeprot); 
	 
	 [Arefprot,Brefprot,Crefprot]=generate_local_matrices(A,Bprot,Cprot,areas,...
															E_0,E_inf,E_omega, ...
															E_vac);
	 xrefprot = gmres(Arefprot,Brefprot*protonatedCurData.q,[],1e-6,300);
	 Erefprot(i,E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
		  (protonatedCurData.q' * Crefprot * xrefprot);
	 
	 [Arefdeprot,Brefdeprot,Crefdeprot]=generate_local_matrices(A,Bdeprot,Cdeprot,areas,...
															E_0,E_inf,E_omega, ...
															E_vac);
	 xrefdeprot = gmres(Arefdeprot,Brefdeprot*deprotonatedCurData.q,[],1e-6,300);
	 Erefdeprot(i,E_omega_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
		  (deprotonatedCurData.q' * Crefdeprot * xrefdeprot);
  end % index i over distances

end % E_omega loop

end % lambda loop
save fig5

end % docomputation
return
