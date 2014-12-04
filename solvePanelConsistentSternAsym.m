function [phiReac, phiBndy, dphiDnBndy] = ...
    solvePanelConsistentSternAsym(dielSurfData, sternSurfData, pqr, bem, ...
			     epsIn, epsOut, kappa, convFactor, ...
				  asymParams, asymBem)

picardIterLimit = 10;
maxGMRESIter = min(100, size(bem.A,1));
numDielPanels = length(dielSurfData.areas);
numSternPanels = length(sternSurfData.areas);
x = zeros(size(bem.A,1),1); 
rhs = bem.B * pqr.q;

% initial guesses for Cauchy data.  We could use Coulomb field 
phiBndy = x(1:numDielPanels);
dphiDnBndy = x(numDielPanels+1:2*numDielPanels);

for picardIter = 1:picardIterLimit
  [curA, curP] = calcPanelYoonSternAsymPert(asymParams, asymBem, dielSurfData, sternSurfData, bem, ...
				      pqr, epsIn, epsOut, kappa, phiBndy, dphiDnBndy);
  [L,U]=lu(curP);
  [x, flag, relres, iter, resvec] = gmres(curA, rhs, [], 1e-5, ...
					  maxGMRESIter, L,U);

  phiBndy = x(1:numDielPanels);
  dphiDnBndy = x(numDielPanels+1:2*numDielPanels);

end

phiReac = convFactor * bem.C * x;