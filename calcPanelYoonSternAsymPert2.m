function [A, P] = calcPanelYoonSternAsymPert2(asymParams, asymBem1, asymBem2, ...
					diel1SurfData,stern1SurfData, ...
                    diel2SurfData,stern2SurfData, ...
                    bem, pqr1, pqr2, epsIn1, epsIn2, ...
					epsOut, kappa, ...
                    phiDiel1Bndy, dphiDnDiel1Bndy, ...
                    phiStern1Bndy, dphiDnStern1Bndy, ...
                    phiDiel2Bndy, dphiDnDiel2Bndy, ...
                    phiStern2Bndy, dphiDnStern2Bndy)

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

r1 = -sum(pqr1.q)/epsOut;
r11 = sum(stern1SurfData.areas' .* dphiDnStern1Bndy);

r2 = -sum(pqr2.q)/epsOut;
r22 = sum(stern2SurfData.areas' .* dphiDnStern2Bndy);


if (r1==r11) || (abs(r1)<1e-2)
  r1Overr11 = 1;
else
  r1Overr11 = r1/r11; %epsOut/epsOut;
end

if (r2==r22) || (abs(r2)<1e-2)
  r2Overr22 = 1;
else
  r2Overr22 = r2/r22; %epsOut/epsOut;
end




A = [bem.A11 bem.A12 bem.A13 bem.A14 bem.A15 bem.A16 bem.A17 bem.A18;
     bem.A21 bem.A22_base*diag(f1./(1+f1)) bem.A23 bem.A24 bem.A25 bem.A26 bem.A27 bem.A28;
     bem.A31 bem.A32_base*diag(f1./(1+f1)) bem.A33 bem.A34 bem.A35 bem.A36 bem.A37 bem.A38;
     bem.A41 bem.A42 bem.A43 bem.A44_base*(r1Overr11) bem.A45 bem.A46 bem.A47 bem.A48_base*(r2Overr22);
     bem.A51 bem.A52 bem.A53 bem.A54 bem.A55 bem.A56 bem.A57 bem.A58;
     bem.A61 bem.A62 bem.A63 bem.A64 bem.A65 bem.A66_base*diag(f2./(1+f2)) bem.A67 bem.A68;
     bem.A71 bem.A72 bem.A73 bem.A74 bem.A75 bem.A76_base*diag(f2./(1+f2)) bem.A77 bem.A78;
     bem.A81 bem.A82 bem.A83 bem.A84_base*(r1Overr11) bem.A85 bem.A86 bem.A87 bem.A88_base*(r2Overr22)];
     
if nargout > 1
    P = sparse([diag(diag(bem.A11)) diag(diag(bem.A12)) 0*bem.A13 0*bem.A14 0*bem.A15 0*bem.A16 0*bem.A17 0*bem.A18;
         diag(diag(bem.A21)) diag(diag(bem.A22_base))*diag(f1./(1+f1)) 0*bem.A23 0*bem.A24 0*bem.A25 0*bem.A26 0*bem.A27 0*bem.A28;
         0*bem.A31 0*bem.A32_base diag(diag(bem.A33)) diag(diag(bem.A34)) 0*bem.A35 0*bem.A36 0*bem.A37 0*bem.A38;
         0*bem.A41 0*bem.A42 diag(diag(bem.A43)) diag(diag(bem.A44_base))*diag(r1Overr11) 0*bem.A45 0*bem.A46 0*bem.A47 0*bem.A48_base;
         0*bem.A51 0*bem.A52 0*bem.A53 0*bem.A54 diag(diag(bem.A55)) diag(diag(bem.A56)) 0*bem.A57 0*bem.A58;
         0*bem.A61 0*bem.A62 0*bem.A63 0*bem.A64 diag(diag(bem.A65)) diag(diag(bem.A66_base))*diag(f2./(1+f2)) 0*bem.A67 0*bem.A68;
         0*bem.A71 0*bem.A72 0*bem.A73 0*bem.A74 0*bem.A75 0*bem.A76_base diag(diag(bem.A77)) diag(diag(bem.A78));
         0*bem.A81 0*bem.A82 0*bem.A83 0*bem.A84_base 0*bem.A85 0*bem.A86 diag(diag(bem.A87)) diag(diag(bem.A88_base))*diag(r2Overr22)]);
end     
     