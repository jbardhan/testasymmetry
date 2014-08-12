doComputation = 1;
if doComputation
global quadrule_order quadrule_x quadrule_w;
quadrule_order = 10;
[quadrule_x,quadrule_w] = setupquad(quadrule_order);

loadconstants;
global E_omega E_sigma kappa;
E_inf   = 1.8;  % see note on AH thesis page 57 (65 of PDF)
E_omega_vec = [1:0.5:20];
E_omega = 1;
E_vac   = 1;
E_sigma = 80;
lambda = 15;  
bigLambda = lambda * sqrt(E_inf/E_sigma);
kappa = 1/bigLambda;

prot_radius    = 13.4; % following Gong, Hocky, and Freed (PNAS 2008)
ion_radius     = 1.4;
ion_separation = 2*ion_radius;  % simple 
distance = [0 (prot_radius - ion_radius)];

% now handle the protein-sized sphere
singleIonPQR = '../geometry/protsphere/singleIon.pqr';
singleIonSRF = '../geometry/protsphere/singleIon.srf';
pqrData = readpqr(singleIonPQR);

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
for E_index = 1:length(E_omega_vec)
  E_omega = E_omega_vec(E_index);
  [Anl,Bnl,Cnl]=generate_nonlocal_matrices(A,B,C,A_Y,areas,E_0, ...
														 E_inf,E_omega,E_sigma);
  xnl = gmres(Anl,Bnl*pqrData.q,[],1e-6,100);
  Enl_sweepEps(i,E_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cnl * xnl); 

  [Al,Bl,Cl]=generate_local_matrices(A,B,C,areas,E_0,E_inf,E_omega,E_sigma);
  xl = gmres(Al,Bl*pqrData.q,[],1e-6,100);
  El_sweepEps(i,E_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cl * xl); 

  [Aref,Bref,Cref]=generate_local_matrices(A,B,C,areas,...
														 E_0,E_inf,E_omega, ...
														 E_vac);
  xref = gmres(Aref,Bref*pqrData.q,[],1e-6,100);
  Eref_sweepEps(i,E_index) =  0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * Cref * xref);

end
  if mod(i,10)==0
	 fprintf('Done with %d positions.\n',i);
  end

end

save fig2inset
end

