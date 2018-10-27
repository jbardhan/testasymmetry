function [phiReac, sigma] = solveAsymmetric2(surf1data, surf12data, surfsurf1op, surfsurf2op, ...
					    chargesurfop1, chargesurfop2, bem, epsIn1, epsIn2, ...
					    epsOut, convFactor, pqr1,pqr2, asymParams)

% this is the simpler, not self-consistent version (self-consistency
% means adding bem.A*sigma to asymPert and looping)
maxIter = 100;
sigma = zeros(size(bem.A,1),1);
rhs = bem.B * [pqr1.q;pqr2.q;];
asymPert = calcAsymmetricPert2(asymParams, surf1data, surf2data, surfsurfop1, surfsurfop2, ...
			      chargesurfop1,chargesurfop2, pqr1, pqr2, sigma1, sigma2);
[sigma,flag,relres,iter,resvec] = gmres(bem.A+asymPert, rhs, [], 1e-5, maxIter);
phiReac = convFactor * bem.C * sigma;

if iter >= maxIter
  keyboard;
end
