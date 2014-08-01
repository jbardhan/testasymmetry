function deltaA = calcAsymmetricPert(asymParams, surfdata, surfsurfop, ...
				     chargesurfop, pqr, sigma)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);

Efield = -chargesurfop.dphidnCoul * pqr.q - surfsurfop.K'*sigma;% - 0.5*sigma;
deltaA = diag(surfdata.weights.*(alpha*(tanh(beta*Efield-EfieldOffset)) ...
				 +deltaOffset));
