doComputation = 0;
if doComputation
global quadrule_order quadrule_x quadrule_w;
quadrule_order = 10;
[quadrule_x,quadrule_w] = setupquad(quadrule_order);

loadconstants;
global E_omega E_sigma kappa;
E_inf   = 1.8;  % see note on AH thesis page 57 (65 of PDF)
E_omega = 4;
E_omega_highdiel = 15;
E_vac   = 1;
E_sigma = 80;
lambda = 15;  % 23 A is referenced by AH
bigLambda = lambda * sqrt(E_inf/E_sigma);
kappa = 1/bigLambda;

prot_radius    = 13.4; % all following Gong, Hocky, and Freed (PNAS 2008)
ion_radius     = 1.4;
ion_separation = 2*ion_radius;  % simplest possible model
distance = [6 prot_radius-ion_radius];

protFilename='../geometry/protsphere/protsphere_0.25.srf';
[meshBase, rootDir] = readsrf(protFilename);
meshData = readmesh(meshBase,1);
[centroids,normals,areas] = genmeshcolloc(meshData);
[A] = collocation_mesh(meshData,centroids,normals,areas);
[A_Y] = collocation_mesh_Yukawa(meshData,centroids,normals,areas,kappa);

% salt bridge perpendicular to radial direction
ionPairPQR = '../geometry/protsphere/ionPair.pqr';
ionPairSRF = '../geometry/protsphere/ionPair10.srf';
pqrData = readpqr(ionPairPQR);
end
%E_nl_part1_ionpair = solvEnergyNonlocal(ionPairPQR,ionPairSRF,E_0,E_inf,E_omega,E_sigma,kappa);
for i=1:length(distance)
  pqrData.xyz(1,3) = distance(i);
  pqrData.xyz(2,3) = distance(i);

  gridSpacing = 0.5;
  [pqrDataTemp,Ygrid,Zgrid] = addPqrGrid(prot_radius, pqrData, gridSpacing);

  [B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrDataTemp);

  [Anl,Bnl,Cnl]=generate_nonlocal_matrices(A,B,C,A_Y,areas,E_0, ...
														 E_inf,E_omega,E_sigma);
  xnl = gmres(Anl,Bnl*pqrDataTemp.q,[],1e-6,100);
  phinl_pair(:,i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (Cnl * xnl)/joulesPerCalorie; 
  Enl_pair(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrDataTemp.q' * Cnl * xnl);       % kJ/mol
  
  [Al,Bl,Cl]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega,E_sigma);
  xl = gmres(Al,Bl*pqrDataTemp.q,[],1e-6,100);
  phil_pair(:,i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (Cl * xl)/joulesPerCalorie; 
  El_pair(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrDataTemp.q' * Cl * xl);         % in kJ/mol

  [Aref,Bref,Cref]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega,E_vac);
  xref = gmres(Aref,Bref*pqrDataTemp.q,[],1e-6,100);
  Eref_pair(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrDataTemp.q' * Cref * xref);     % in kJ/mol

  [Al,Bl,Cl]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega_highdiel,E_sigma);
  xl = gmres(Al,Bl*pqrDataTemp.q,[],1e-6,100);
  phil_pair_highdiel(:,i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (Cl * xl)/joulesPerCalorie; 
  El_pair_highdiel(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrDataTemp.q' * Cl * xl);         % in kJ/mol

  [Aref,Bref,Cref]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega_highdiel,E_vac);
  xref = gmres(Aref,Bref*pqrDataTemp.q,[],1e-6,100);
  Eref_pair_highdiel(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrDataTemp.q' * Cref * xref);     % in kJ/mol

  fprintf('Done %i positions...\n',i);

end


