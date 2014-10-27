function [A, P] = calcYoonNoSternAsymPert(asymParams, surfdata, bem, ...
				       pqr, epsIn, epsOut, phiBndy, dphiDnBndy)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);

sigmaEff = bem.dielDielOp.V\(bem.dielDielOp.V * dphiDnBndy - bem.dielDielOp.K*phiBndy);
Efield = -bem.dielChargeOp.dphidnCoul * pqr.q - bem.dielDielOp.K'*sigmaEff;

h = (alpha*(tanh(beta*Efield-EfieldOffset)) +deltaOffset);
f = (epsIn/(epsOut-epsIn)) - h;

A = [bem.A11 bem.A12; bem.A21 bem.A22_base*diag(f./(1+f))];
P = [diag(diag(bem.A11)) diag(diag(bem.A12));
     diag(diag(bem.A21)) diag(diag(bem.A22_base).*(f./(1+f)))];
