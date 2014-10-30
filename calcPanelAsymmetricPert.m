function deltaA = calcPanelAsymmetricPert(asymParams, asymBem, surfdata, surfsurfop, ...
				     chargesurfop, pqr, sigma)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);

Efield = -asymBem.dphidnCoul_NLBC * pqr.q - asymBem.Kp_NLBC*sigma;% - 0.5*sigma;
h = (alpha*(tanh(beta*Efield-EfieldOffset)) +deltaOffset);

deltaA = diag(asymBem.h_quadrature*h);
