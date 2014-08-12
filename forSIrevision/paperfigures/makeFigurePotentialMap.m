doComputation = 1;
if 1
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
distance = [ 6 prot_radius-ion_radius];

% burying a single ion at varying depths: what is the solvation
% energy of that ion alone
singleIonPQR = '../geometry/protsphere/singleIon.pqr';
singleIonSRF = '../geometry/protsphere/singleIon.srf';
pqrData = readpqr(singleIonPQR);
%E_nl_part0_singleIon = solvEnergyNonlocal(singleIonPQR,singleIonSRF,E_0,E_inf,E_omega,E_sigma,kappa);
%E_nl_part0_singleIon_vac = solvEnergyLocal(singleIonPQR,singleIonSRF,E_0,E_inf,E_omega,E_vac);
%E_nl_part0_singleIon_solv = E_nl_part0_singleIon - E_nl_part0_singleIon_vac;

% now handle the protein-sized sphere
protFilename='../geometry/protsphere/protsphere_0.25.srf';
protFilename2='../geometry/protsphere/protsphere_0.5.srf';  % used later
[meshBase, rootDir] = readsrf(protFilename);
meshData = readmesh(meshBase,1);
[centroids,normals,areas] = genmeshcolloc(meshData);
[A] = collocation_mesh(meshData,centroids,normals,areas);
[A_Y] = collocation_mesh_Yukawa(meshData,centroids,normals,areas,kappa);
end
for i=1:length(distance)
  pqrData.xyz(1,3) = distance(i);

  gridSpacing = 0.5;
  [pqrDataTemp,Ygrid,Zgrid] = addPqrGrid(prot_radius, pqrData, gridSpacing);

  [B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrDataTemp);

  [Anl,Bnl,Cnl]=generate_nonlocal_matrices(A,B,C,A_Y,areas,E_0, ...
														 E_inf,E_omega,E_sigma);
  xnl = gmres(Anl,Bnl*pqrDataTemp.q,[],1e-6,100);
  phinl_part0(:,i) =  (Na/1000) * (q^2/E_0) * 1e10 *...
	 (Cnl * xnl)/joulesPerCalorie; 
  Enl_part0(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrDataTemp.q' * Cnl * xnl); 

  [Al,Bl,Cl]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega,E_sigma);
  xl = gmres(Al,Bl*pqrDataTemp.q,[],1e-6,100);
  phil_part0(:,i) =  (Na/1000) * (q^2/E_0) * 1e10 *...
	 (Cl * xl)/joulesPerCalorie; 
  El_part0(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrDataTemp.q' * Cl * xl); 

  [Aref,Bref,Cref]=generate_local_matrices(A,B,C,areas,...
														 E_0,E_inf,E_omega, ...
														 E_vac);
  xref = gmres(Aref,Bref*pqrDataTemp.q,[],1e-6,100);
  Eref_part0(i) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrDataTemp.q' * Cref * xref);

  if mod(i,10)==0
	 fprintf('Done with %d positions.\n',i);
  end
end

save ionPotentials