function [phiReac, phiBndy, dphiDnBndy] = ...
    solveConsistentSternAsym(dielSurfData, sternSurfData, pqr, bem, ...
			     epsIn, epsOut, kappa, convFactor, asymParams)

picardIterLimit = 10;
maxGMRESIter = min(100, size(bem.A,1));
numDielPanels = length(dielSurfData.weights);
numSternPanels = length(sternSurfData.weights);
x = zeros(size(bem.A,1),1); 
rhs = bem.B * pqr.q;

% initial guesses for Cauchy data.  We could use Coulomb field 
phiBndy = x(1:numDielPanels);
dphiDnBndy = x(numDielPanels+1:2*numDielPanels);

for picardIter = 1:picardIterLimit
  [curA, curP] = calcYoonSternAsymPert(asymParams, dielSurfData, sternSurfData, bem, ...
				      pqr, epsIn, epsOut, kappa, phiBndy, dphiDnBndy);

  [x, flag, relres, iter, resvec] = gmres(curA, rhs, [], 1e-5, ...
					  maxGMRESIter, curP);

  phiBndy = x(1:numDielPanels);
  dphiDnBndy = x(numDielPanels+1:2*numDielPanels);

end

phiReac = convFactor * bem.C * x;