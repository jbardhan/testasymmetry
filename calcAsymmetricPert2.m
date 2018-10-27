function deltaA = calcAsymmetricPert2(asymParams, surf1data, surf2data, surfsurfop1, surfsurfop2, ...
				     chargesurfop1,chargesurfop2, pqr1, pqr2, sigma1, sigma2)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);

Efield1 = -chargesurfop1.dphidnCoul1 * pqr1.q - surfsurfop1.K'*sigma1;% - 0.5*sigma;
Efield2 = -chargesurfop2.dphidnCoul2 * pqr2.q - surfsurfop2.K'*sigma2;% - 0.5*sigma;
h1 = (alpha*(tanh(beta*Efield1-EfieldOffset)) +deltaOffset);
h2 = (alpha*(tanh(beta*Efield2-EfieldOffset)) +deltaOffset);
deltaA = diag([diag(surf1data.weights1.*h1);diag(surf2data.weights2.*h2)]);
