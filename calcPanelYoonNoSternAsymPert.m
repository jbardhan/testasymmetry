function [A, P] = calcPanelYoonNoSternAsymPert(asymParams, asymBem, ...
					       surfdata, bem, pqr, ...
					       epsIn, epsOut, phiBndy, ...
					       dphiDnBndy)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);

I = eye(length(phiBndy));
sigmaEff = bem.dielDielOp.V\(bem.dielDielOp.V * dphiDnBndy - (-0.5*I+bem.dielDielOp.K)*phiBndy);

% calculating h is then the same as calcPanelAsymmetricPert
Efield = -asymBem.dphidnCoul_NLBC * pqr.q - asymBem.Kp_NLBC*sigmaEff;
h = (alpha*(tanh(beta*Efield-EfieldOffset)) +deltaOffset);

% last, need to convert to f(En) and reform our matrices
f = (epsIn/(epsOut-epsIn)) - h;

A = [bem.A11 bem.A12; bem.A21 bem.A22_base*diag(f./(1+f))];
P = [diag(diag(bem.A11)) diag(diag(bem.A12));
     diag(diag(bem.A21)) diag(diag(bem.A22_base).*(f./(1+f)))];
