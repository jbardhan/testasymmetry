function [A,B,C] = generate_local_matrices(surfaceToSurfaceOperators,...
														 chargeToSurfaceOperators,...
														 surfaceToChargeOperators,...
														 areas,...
														 E_0,E_inf,E_in,...
														 E_out)
np = length(areas);
I = diag(ones(np,1));

A =  0.5 * (1 + E_in/E_out) * I;
A = A + (E_in/E_out - 1) * surfaceToSurfaceOperators.K;


Bmat1 = (surfaceToSurfaceOperators.K - 0.5 * I) * ...
		 chargeToSurfaceOperators.DirichletTrace;
Bmat2 = -1 * (E_in / E_out) * surfaceToSurfaceOperators.V * ...
		 chargeToSurfaceOperators.NeumannTrace;
B = (Bmat1 + Bmat2)/E_in;


DirichletToNeumannMap = surfaceToSurfaceOperators.V \ (0.5*I + surfaceToSurfaceOperators.K);
Cmat1 = surfaceToChargeOperators.slpToCharges * DirichletToNeumannMap;
Cmat2 = -1 * surfaceToChargeOperators.dlpToCharges;

C = Cmat1 + Cmat2;

%keyboard