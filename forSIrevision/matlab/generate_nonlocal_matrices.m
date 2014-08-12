function [A,B,C] = generate_nonlocal_matrices(surfaceToSurfaceOperators,...
														 chargeToSurfaceOperators,...
														 surfaceToChargeOperators,...
														 surfaceToSurfaceOperators_Y,...
														 areas,...
														 E_0,E_inf,E_in,...
														 E_out)
np = length(areas);
I = diag(ones(np,1));
sigma = 0.5 * I;

V  = surfaceToSurfaceOperators.V;
K  = surfaceToSurfaceOperators.K;
VY = surfaceToSurfaceOperators_Y.V;
KY = surfaceToSurfaceOperators_Y.K;
%keyboard
%adjacent to Eq. 5.192, p121 in Hildebrandt thesis
Bmat1 = -(I - sigma - KY + (E_in/E_out)*(KY-K) ) ...
		  * chargeToSurfaceOperators.DirichletTrace/E_in;
Bmat2 = -((E_in/E_inf)*VY - (E_in/E_out)*(VY-V) ) ...
		  * chargeToSurfaceOperators.NeumannTrace/E_in;
B = [Bmat1+Bmat2; 
	  0 * chargeToSurfaceOperators.DirichletTrace; 
	  0 * chargeToSurfaceOperators.DirichletTrace];


% A
A11 = I - sigma - KY;
A21 = sigma + K;
A31 = 0 * I;

A12 = (E_in/E_inf) * VY	- (E_in/E_out) * (VY-V);
A22 = -V;
A32 = (E_in/E_inf)*V;

A13 = (E_inf/E_out) * (KY-K);
A23 = 0 * I;
A33 = I - sigma - K;

A = [A11 A12 A13; A21 A22 A23; A31 A32 A33;];

% i believe that the same expression as the local electrostatic
% model holds to compute the reaction field inside
%     not needed: DirichletToNeumannMap = V \ (0.5*I + K);
Cmat1 = surfaceToChargeOperators.slpToCharges;
Cmat2 = -1 * surfaceToChargeOperators.dlpToCharges;
C = [Cmat2 Cmat1 0*surfaceToChargeOperators.slpToCharges];
