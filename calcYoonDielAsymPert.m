function [A, P] = calcYoonDielAsymPert(asymParams, surfdata, bem, ...
				       pqr, epsIn, epsOut, phiBndy, dphiDnBndy)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);

Efield = -bem.dielChargeOp.dphidnCoul * pqr.q - bem.dielDielOp.K'*dphiDnBndy ...
	 + bem.dielDielOp.W*phiBndy;

h = (alpha*(tanh(beta*Efield-EfieldOffset)) +deltaOffset);
f = (epsIn/(epsOut-epsIn)) - h;

A = [bem.A11 bem.A12; bem.A21 (f/(1+f))*bem.A22_base];
P = [diag(diag(bem.A11)) diag(diag(bem.A12));
     diag(diag(bem.A21)) diag(diag(bem.A22_base).*(f./(1+f)))];
