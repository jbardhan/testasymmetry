function [phiReac, sigma] = solveConsistentAsymmetric(surfdata, ...
						  surfsurfop, ...
						  chargesurfop, bem, ...
						  epsIn, epsOut, ...
						  convFactor, pqr, asymParams)

picardIterLimit = 10;
maxGMRESIter = min(100,size(bem.A,1)); % to avoid warning from gmres.
sigma = zeros(size(bem.A,1),1);
rhs = bem.B * pqr.q;

for picardIter = 1:picardIterLimit
  asymPert = calcAsymmetricPert(asymParams, surfdata, surfsurfop, ...
				chargesurfop, pqr, sigma);
  [sigma,flag,relres,iter,resvec] = gmres(bem.A+asymPert, rhs, [], 1e-5, maxGMRESIter);
  phiReac = convFactor * bem.C * sigma;
  
  if iter >= maxGMRESIter
    keyboard;
  end
end
