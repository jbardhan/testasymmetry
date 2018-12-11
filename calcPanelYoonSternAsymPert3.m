function [A, P] = calcPanelYoonSternAsymPert3(asymParams, asymBem1, asymBem2, ...
					diel1SurfData,diel2SurfData,stern3SurfData, ...
                    bem, pqr1, pqr2, epsIn1, epsIn2, ...
					epsOut, kappa, ...
                    phiDiel1Bndy, dphiDnDiel1Bndy, ...
                    phiDiel2Bndy, dphiDnDiel2Bndy, ...
                    phiStern3Bndy, dphiDnStern3Bndy)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
%deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);
deltaOffset  = asymParams.mu;
I1 = eye(length(diel1SurfData.areas));
I2 = eye(length(diel2SurfData.areas));

sigmaEff1 = bem.diel1Diel1Op.V\(bem.diel1Diel1Op.V * dphiDnDiel1Bndy - (-0.5*I1+bem.diel1Diel1Op.K)*phiDiel1Bndy);
Efield1 = -asymBem1.dphidnCoul_NLBC * pqr1.q - asymBem1.Kp_NLBC*sigmaEff1; 

sigmaEff2 = bem.diel2Diel2Op.V\(bem.diel2Diel2Op.V * dphiDnDiel2Bndy - (-0.5*I2+bem.diel2Diel2Op.K)*phiDiel2Bndy);
Efield2 = -asymBem2.dphidnCoul_NLBC * pqr2.q - asymBem2.Kp_NLBC*sigmaEff2; 

 
 
h1 = (alpha*(tanh(beta*Efield1-EfieldOffset)) +deltaOffset);
f1 = (epsIn1/(epsOut-epsIn1)) - h1;
 
h2 = (alpha*(tanh(beta*Efield2-EfieldOffset)) +deltaOffset);
f2 = (epsIn2/(epsOut-epsIn2)) - h2;

r1 = -sum(sum(pqr1.q) + sum(pqr2.q))/epsOut;
r11 = sum(stern3SurfData.areas' .* dphiDnStern3Bndy);


if (r1==r11) || (abs(r1)<1e-2)
  r1Overr11 = 1;
else
  r1Overr11 = r1/r11; %epsOut/epsOut;
end




A = [bem.A11 bem.A12 bem.A13 bem.A14 bem.A15 bem.A16;
     bem.A21 bem.A22_base*diag(f1./(1+f1)) bem.A23 bem.A24_base*diag(f2./(1+f2)) bem.A25 bem.A26;
     bem.A31 bem.A32 bem.A33 bem.A34 bem.A35 bem.A36;
     bem.A41 bem.A42_base*diag(f1./(1+f1)) bem.A43 bem.A44_base*diag(f2./(1+f2)) bem.A45 bem.A46;
     bem.A51 bem.A52_base*diag(f1./(1+f1)) bem.A53 bem.A54_base*diag(f2./(1+f2)) bem.A55 bem.A56;
     bem.A61 bem.A62 bem.A63 bem.A64 bem.A65 bem.A66_base*diag(r1Overr11)];
     
if nargout > 1
    P = sparse([diag(diag(bem.A11)) diag(diag(bem.A12)) 0*bem.A13 0*bem.A14 0*bem.A15 0*bem.A16;
         diag(diag(bem.A21)) diag(diag(bem.A22_base))*diag(f1./(1+f1)) 0*bem.A23 0*bem.A24_base 0*bem.A25 0*bem.A26;
         0*bem.A31 0*bem.A32 diag(diag(bem.A33)) diag(diag(bem.A34)) 0*bem.A35 0*bem.A36;
         0*bem.A41 0*bem.A42_base diag(diag(bem.A43)) diag(diag(bem.A44_base))*diag(f2./(1+f2)) 0*bem.A45 0*bem.A46;
         0*bem.A51 0*bem.A52_base 0*bem.A53 0*bem.A54_base diag(diag(bem.A55)) diag(diag(bem.A56));
         0*bem.A61 0*bem.A62 0*bem.A63 0*bem.A64 diag(diag(bem.A65)) diag(diag(bem.A66_base))*diag(r1Overr11)]);
end     
     