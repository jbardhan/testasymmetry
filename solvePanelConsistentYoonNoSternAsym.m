function [phiReac, phiBndy, dphiDnBndy] = ...
    solvePanelConsistentYoonNoSternAsym(surfdata, bem, epsIn, epsOut, ...
					convFactor, pqr, asymParams, ...
					asymBem)

picardIterLimit = 5;
maxGMRESIter = min(100, size(bem.A,1));
numPanels = length(surfdata.areas);
x = zeros(size(bem.A,1),1);
rhs = bem.B * pqr.q;

phiBndy = x(1:numPanels);
dphiDnBndy = x(numPanels+1:end);

for picardIter = 1:picardIterLimit
  [curA, curP] = calcPanelYoonNoSternAsymPert(asymParams, asymBem, surfdata, ...
					      bem, pqr, epsIn, epsOut, ...
					      phiBndy, dphiDnBndy);
  [L,U]=lu(curP);
  [x, flag, relres, iter, resvec] = gmres(curA, rhs, [], 1e-5, ...
					  maxGMRESIter, L,U);

  phiBndy = x(1:numPanels);
  dphiDnBndy = x(numPanels+1:end);
end

phiReac = convFactor * bem.C * x;
