function [phiReac, sigma] = solveAsymmetric(surfdata, surfsurfop, ...
					    chargesurfop, bem, epsIn, ...
					    epsOut, convFactor, pqr, asymParams)

% this is the simpler, not self-consistent version (self-consistency
% means adding bem.A*sigma to asymPert and looping)
maxIter = 100;
sigma = zeros(size(bem.A,1),1);
rhs = bem.B * pqr.q;
asymPert = calcAsymmetricPert(asymParams, surfdata, surfsurfop, ...
			      chargesurfop, pqr, sigma);
[sigma,flag,relres,iter,resvec] = gmres(bem.A+asymPert, rhs, [], 1e-5, maxIter);
phiReac = convFactor * bem.C * sigma;

if iter >= maxIter
  keyboard;
end
