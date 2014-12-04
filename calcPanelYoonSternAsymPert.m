function [A, P] = calcPanelYoonSternAsymPert(asymParams, asymBem, dielSurfData, ...
					sternSurfData, bem, pqr, ...
					epsIn, epsOut, kappa, phiBndy, ...
					dphiDnBndy)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);
I = eye(length(dielSurfData.areas));

sigmaEff = bem.dielDielOp.V\(bem.dielDielOp.V * dphiDnBndy - (-0.5*I+bem.dielDielOp.K)*phiBndy);
Efield = -asymBem.dphidnCoul_NLBC * pqr.q - asymBem.Kp_NLBC*sigmaEff; 
	 
h = (alpha*(tanh(beta*Efield-EfieldOffset)) +deltaOffset);
f = (epsIn/(epsOut-epsIn)) - h;

A = [bem.A11 bem.A12 bem.A13 bem.A14;
     bem.A21 bem.A22_base*diag(f./(1+f)) bem.A23 bem.A24;
     bem.A31 bem.A32_base*diag(f./(1+f)) bem.A33 bem.A34;
     bem.A41 bem.A42 bem.A43 bem.A44_base*(epsOut/epsOut)];

P = sparse([diag(diag(bem.A11)) diag(diag(bem.A12)) 0*bem.A13 0*bem.A14;
     diag(diag(bem.A21)) diag(diag(bem.A22_base))*diag(f./(1+f)) ...
     0*bem.A23 0*bem.A24;
     0*bem.A31 0*bem.A32_base diag(diag(bem.A33)) diag(diag(bem.A34));
     0*bem.A41 0*bem.A42 diag(diag(bem.A43)) diag(diag(bem.A44_base))*(epsOut/epsOut)]);
