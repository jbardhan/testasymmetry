function [phiReac, phiBndy, dphiDnBndy] = ...
    solveConsistentYoonNoSternAsym(surfdata, bem, epsIn, epsOut, ...
				convFactor, pqr, asymParams)

picardIterLimit = 10;
maxGMRESIter = min(100, size(bem.A,1));
numPanels = length(surfdata.weights);
x = zeros(size(bem.A,1),1); 
rhs = bem.B * pqr.q;

% initial guesses for Cauchy data.  We could use Coulomb field 
phiBndy = x(1:numPanels);
dphiDnBndy = x(numPanels+1:end);

for picardIter = 1:picardIterLimit
  [curA, curP] = calcYoonNoSternAsymPert(asymParams, surfdata, bem, ...
				      pqr, epsIn, epsOut, phiBndy, dphiDnBndy);

  [x, flag, relres, iter, resvec] = gmres(curA, rhs, [], 1e-5, ...
					  maxGMRESIter, curP);

  phiBndy = x(1:numPanels);
  dphiDnBndy = x(numPanels+1:end);

end

phiReac = convFactor * bem.C * x;