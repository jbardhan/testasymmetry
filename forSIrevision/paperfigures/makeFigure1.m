doComputation = 1;
if doComputation
global quadrule_order quadrule_x quadrule_w;
quadrule_order = 10;
[quadrule_x,quadrule_w] = setupquad(quadrule_order);

loadconstants;
global E_omega E_sigma kappa;
E_inf   = 1.8;  % see note on AH thesis page 57 (65 of PDF)
E_omega = 4;
E_vac   = 1;
E_sigma = 80;
lambda = 15;  
bigLambda = lambda * sqrt(E_inf/E_sigma);
kappa = 1/bigLambda;

prot_radius    = 13.4; % following Gong, Hocky, and Freed (PNAS 2008)
ion_radius     = 1.4;
ion_separation = 2*ion_radius;  % simple 
distance = 0;%:0.1:(prot_radius - ion_radius);

% burying a single ion at varying depths: what is the solvation
% energy of that ion alone
singleIonPQR = '../geometry/protsphere/singleIon.pqr';
singleIonSRF = '../geometry/protsphere/singleIon.srf';
pqrData = readpqr(singleIonPQR);
E_nl_part0_singleIon = solvEnergyNonlocal(singleIonPQR,singleIonSRF,E_0,E_inf,E_omega,E_sigma,kappa);
E_nl_part0_singleIon_vac = solvEnergyLocal(singleIonPQR,singleIonSRF,E_0,E_inf,E_omega,E_vac);
E_nl_part0_singleIon_solv = E_nl_part0_singleIon - E_nl_part0_singleIon_vac;

% now handle the protein-sized sphere
protFilename='../geometry/protsphere/protsphere_0.25.srf';
protFilename2='../geometry/protsphere/protsphere_0.5.srf';  % used later
[meshBase, rootDir] = readsrf(protFilename);
meshData = readmesh(meshBase,1);
[centroids,normals,areas] = genmeshcolloc(meshData);
[A] = collocation_mesh(meshData,centroids,normals,areas);
[A_Y] = collocation_mesh_Yukawa(meshData,centroids,normals,areas,kappa);

for i=1:length(distance)
  pqrData.xyz(1,3) = distance(i);
  [B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrData);

  [Anl,Bnl,Cnl]=generate_nonlocal_matrices(A,B,C,A_Y,areas,E_0, ...
														 E_inf,E_omega,E_sigma);
  xnl = gmres(Anl,Bnl*pqrData.q,[],1e-6,100);
  Enl_part0(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cnl * xnl); 

  [Al,Bl,Cl]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega,E_sigma);
  xl = gmres(Al,Bl*pqrData.q,[],1e-6,100);
  El_part0(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cl * xl); 

  [Aref,Bref,Cref]=generate_local_matrices(A,B,C,areas,...
														 E_0,E_inf,E_omega, ...
														 E_vac);
  xref = gmres(Aref,Bref*pqrData.q,[],1e-6,100);
  Eref_part0(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cref * xref);

  if mod(i,10)==0
	 fprintf('Done with %d positions.\n',i);
  end
end


% 4X slower calculations that use twice as many boundary elements
[meshBase2, rootDir2] = readsrf(protFilename2);
meshData2 = readmesh(meshBase2,1);
[centroids2,normals2,areas2] = genmeshcolloc(meshData2);
[A2] = collocation_mesh(meshData2,centroids2,normals2,areas2);
[A_Y2] = collocation_mesh_Yukawa(meshData2,centroids2,normals2,areas2,kappa);
distance2 = 0:0.5:(prot_radius-ion_radius);  % five times fewer 
for i=1:length(distance2)
  pqrData.xyz(1,3) = distance2(i);
  [B2,C2] = chargeCollocation_mesh(meshData2,centroids2,normals2,areas2,pqrData);

  [Anl2,Bnl2,Cnl2]=generate_nonlocal_matrices(A2,B2,C2,A_Y2,areas2,E_0, ...
														 E_inf,E_omega,E_sigma);
  xnl2 = gmres(Anl2,Bnl2*pqrData.q,[],1e-6,100);
  Enl_part0_highres(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cnl2 * xnl2); 

  [Al2,Bl2,Cl2]=generate_local_matrices(A2,B2,C2,areas2,E_0,E_inf,E_omega,E_sigma);
  xl2 = gmres(Al2,Bl2*pqrData.q,[],1e-6,100);
  El_part0_highres(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cl2 * xl2); 

  [Aref2,Bref2,Cref2]=generate_local_matrices(A2,B2,C2,areas2,...
														 E_0,E_inf,E_omega, ...
														 E_vac);
  xref2 = gmres(Aref2,Bref2*pqrData.q,[],1e-6,100);
  Eref_part0_highres(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cref2 * xref2);

  if mod(i,10)==0
	 fprintf('Done with %d positions.\n',i);
  end
end

%save fig2
end
