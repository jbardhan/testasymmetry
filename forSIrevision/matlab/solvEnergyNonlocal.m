function [E_transfer,A,B,C,A_Y,A_nl,B_nl,C_nl] = solvEnergyNonlocal(pqrFile,srfFile,E_0,E_inf, ...
													  E_omega,E_sigma,kappa)
% E_transfer is the energy of transferring from a low-dielectric
% local solvent to the high-dielectric non-local solvent
q = 1.60217646e-19;
Na = 6.0221415e23;

printOn = 1;
if printOn
  fprintf('Reading data...');
end

[meshBase,rootDir] = readsrf(srfFile);
meshData = readmesh(meshBase,1);
[centroids,normals,areas] = genmeshcolloc(meshData);

pqrData = readpqr(pqrFile);

if printOn
  fprintf(' done with %d panels.\n',length(areas));
  fprintf('Generating discretized Laplace operators...');
end

% create Laplace operators
[A] = collocation_mesh(meshData,centroids,normals,areas);
[B,C] = chargeCollocation_mesh(meshData,centroids,normals,areas,pqrData);

if printOn
  fprintf(' done.\n');
  fprintf('Generating discretized Yukawa operators...');
end

% actually compute Yukawa operator
[A_Y] = collocation_mesh_Yukawa(meshData,centroids,normals,areas,kappa);

if printOn
  fprintf(' done.\n');
  fprintf('Forming BEM matrix and solving...');
end

[A_nl, B_nl, C_nl] = generate_nonlocal_matrices(A, B, C, A_Y, areas,...
																E_0, E_inf, E_omega, ...
																E_sigma);

E_transfer = 0.5 * (Na/1000) * (q^2/E_0) * 1e10 *...
	 (pqrData.q' * C_nl * (A_nl \ (B_nl * pqrData.q)));  % comes out
                                                        % in kJ/mol
if printOn
  fprintf(' done.\n');
  fprintf('Returning from solvEnergyNonlocal.\n');
end
