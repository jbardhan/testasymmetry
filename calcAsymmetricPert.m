function deltaA = calcAsymmetricPert(asymParams, surfdata, surfsurfop, ...
				     chargesurfop, pqr)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = asymParams.deltaOffset;

Efield = -chargesurfop.dphidnCoul * pqr.q;
deltaA = diag(surfdata.weights.*(alpha*(tanh(beta*Efield-EfieldOffset)) ...
				 +deltaOffset));
