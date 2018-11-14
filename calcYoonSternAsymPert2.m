function [A, P] = calcYoonSternAsymPert2(asymParams, ...
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
%deltaOffset  = asymParams.mu;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);

Efield1 = -bem.diel1ChargeOp.dphidnCoul * pqr1.q - bem.diel1Diel1Op.K'*dphiDnDiel1Bndy ...
	 + bem.diel1Diel1Op.W*phiDiel1Bndy;

Efield2 = -bem.diel2ChargeOp.dphidnCoul * pqr2.q - bem.diel1Diel2Op.K'*dphiDnDiel2Bndy ...
	 + bem.diel2Diel2Op.W*phiDiel2Bndy; 
 
h1 = (alpha*(tanh(beta*Efield1-EfieldOffset)) +deltaOffset);
f1 = (epsIn1/(epsOut-epsIn1)) - h1;
 
h2 = (alpha*(tanh(beta*Efield2-EfieldOffset)) +deltaOffset);
f2 = (epsIn2/(epsOut-epsIn2)) - h2;

r1 = -sum(pqr1.q)/epsOut;
r11 = sum(stern1SurfData.weights .* dphiDnStern1Bndy);
r12 = sum(stern1SurfData.weights .* dphiDnStern2Bndy);

r2 = -sum(pqr2.q)/epsOut;
r21 = sum(stern2SurfData.weights .* dphiDnStern1Bndy);
r22 = sum(stern2SurfData.weights .* dphiDnStern2Bndy);

r1Overr11 = r1/r11; %epsS1/epsOut;
r1Overr12 = r1/r12; %epsS2/epsOut;
r2Overr21 = r2/r21; %epsS1/epsOut;
r2Overr22 = r2/r22; %epsS2/epsOut;

A = [bem.A11 bem.A12 bem.A13 bem.A14 bem.A15 bem.A16 bem.A17 bem.A18;
     bem.A21 bem.A22_base*diag(f1./(1+f1)) bem.A23 bem.A24 bem.A25 bem.A26 bem.A27 bem.A28;
     bem.A31 bem.A32_base*diag(f1./(1+f1)) bem.A33 bem.A34 bem.A35 bem.A36 bem.A37 bem.A38;
     bem.A41 bem.A42 bem.A43 bem.A44_base*(r1Overr11) bem.A45 bem.A46 bem.A47 bem.A48_base*(r1Overr12);
     bem.A51 bem.A52 bem.A53 bem.A54 bem.A55 bem.A56 bem.A57 bem.A58;
     bem.A61 bem.A62 bem.A63 bem.A64 bem.A65 bem.A66_base*diag(f2./(1+f2)) bem.A67 bem.A68;
     bem.A71 bem.A72 bem.A73 bem.A74 bem.A75 bem.A76_base*diag(f2./(1+f2)) bem.A77 bem.A78;
     bem.A81 bem.A82 bem.A83 bem.A84_base*(r2Overr21) bem.A85 bem.A86 bem.A87 bem.A88_base*(r2Overr22)];
     

P = sparse([diag(diag(bem.A11)) diag(diag(bem.A12)) 0*bem.A13 0*bem.A14 0*bem.A15 0*bem.A16 0*bem.A17 0*bem.A18;
     diag(diag(bem.A21)) diag(diag(bem.A22_base))*diag(f1./(1+f1)) 0*bem.A23 0*bem.A24 0*bem.A25 0*bem.A26 0*bem.A27 0*bem.A28;
     0*bem.A31 0*bem.A32_base diag(diag(bem.A33)) diag(diag(bem.A34)) 0*bem.A35 0*bem.A36 0*bem.A37 0*bem.A38;
     0*bem.A41 0*bem.A42 diag(diag(bem.A43)) diag(diag(bem.A44_base))*diag(r1Overr11) 0*bem.A45 0*bem.A46 0*bem.A47 0*bem.A48_base;
     0*bem.A51 0*bem.A52 0*bem.A53 0*bem.A54 diag(diag(bem.A55)) diag(diag(bem.A56)) 0*bem.A57 0*bem.A58;
     0*bem.A61 0*bem.A62 0*bem.A63 0*bem.A64 diag(diag(bem.A65)) diag(diag(bem.A66_base))*diag(f2./(1+f2)) 0*bem.A67 0*bem.A68;
     0*bem.A71 0*bem.A72 0*bem.A73 0*bem.A74 0*bem.A75 0*bem.A76_base diag(diag(bem.A77)) diag(diag(bem.A78));
     0*bem.A81 0*bem.A82 0*bem.A83 0*bem.A84_base 0*bem.A85 0*bem.A86 diag(diag(bem.A87)) diag(diag(bem.A88_base))*diag(r2Overr22)]);
     
     
     
     
     
     
     
     
     
     