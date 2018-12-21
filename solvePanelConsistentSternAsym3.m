function [phiReac1, phiReac2, ...
    phiDiel1Bndy, dphiDnDiel1Bndy, phiDiel2Bndy, dphiDnDiel2Bndy, ...
    x, curA] = ...
    solvePanelConsistentSternAsym3(diel1SurfData, diel2SurfData,stern3SurfData,...
    pqr1, pqr2, bem, epsIn1, epsIn2, epsOut, kappa, convFactor, ...
				  asymParams, asymBem1, asymBem2)

picardIterLimit = 30;
maxGMRESIter = min(100, size(bem.A,1));
numDiel1Panels = length(diel1SurfData.areas);
numDiel2Panels = length(diel2SurfData.areas);
numStern3Panels = length(stern3SurfData.areas);
x = zeros(size(bem.A,1),1); 
rhs = bem.B1 * pqr1.q +  bem.B2 * pqr2.q;

% initial guesses for Cauchy data.  We could use Coulomb field 
phiDiel1Bndy = x(1:numDiel1Panels);
dphiDnDiel1Bndy = x(numDiel1Panels+1:2*numDiel1Panels);
phiDiel2Bndy = x(2*numDiel1Panels+1:2*numDiel1Panels+numDiel2Panels);
dphiDnDiel2Bndy = x(2*numDiel1Panels+numDiel2Panels+1:2*(numDiel1Panels+numDiel2Panels));
phiStern3Bndy = x(2*(numDiel1Panels+numDiel2Panels)+1:2*(numDiel1Panels+numDiel2Panels)+numStern3Panels);
dphiDnStern3Bndy = ones(numStern3Panels,1)*(-sum(sum(pqr1.q) + sum(pqr2.q))/ epsOut)/sum(stern3SurfData.areas);

for picardIter = 1:picardIterLimit
  if picardIter ==1
    [curA, curP] = calcPanelYoonSternAsymPert3(asymParams, asymBem1,asymBem2, ...
                          diel1SurfData,diel2SurfData,stern3SurfData,bem, ...
					      pqr1, pqr2, epsIn1, epsIn2, epsOut, kappa,...
					      phiDiel1Bndy,dphiDnDiel1Bndy,...
                          phiDiel2Bndy,dphiDnDiel2Bndy,...
					      phiStern3Bndy,dphiDnStern3Bndy);
    [L,U]=lu(curP);
  else
    curA = calcPanelYoonSternAsymPert3(asymParams, asymBem1, asymBem2, ...
                          diel1SurfData,diel2SurfData,  stern3SurfData, bem, ...
					      pqr1, pqr2, epsIn1, epsIn2, epsOut, kappa, ...
					      phiDiel1Bndy, dphiDnDiel1Bndy, ...
					      phiDiel2Bndy, dphiDnDiel2Bndy, ...
					      phiStern3Bndy, dphiDnStern3Bndy);
  end
  [x, flag, relres, iter, resvec] = gmres(curA, rhs, [], 1e-5, ...
					  maxGMRESIter, L,U,x);
  fprintf('took %d GMRES iters at picard iteration %d\n', max(iter), ...
	  picardIter)
  %  keyboard
  phiDiel1Bndy = x(1:numDiel1Panels);
  dphiDnDiel1Bndy = x(numDiel1Panels+1:2*numDiel1Panels);
  phiDiel2Bndy = x(2*numDiel1Panels+1:2*numDiel1Panels+numDiel2Panels);
  dphiDnDiel2Bndy = x(2*numDiel1Panels+numDiel2Panels+1:2*(numDiel1Panels+numDiel2Panels));
  phiStern3Bndy = x(2*(numDiel1Panels+numDiel2Panels)+1:2*(numDiel1Panels+numDiel2Panels)+numStern3Panels);
  dphiDnStern3Bndy = x(2*(numDiel1Panels+numDiel2Panels)+numStern3Panels+1:2*(numDiel1Panels+numDiel2Panels+numStern3Panels));
end

phiReac1 = convFactor * bem.C1 * x;
phiReac2 = convFactor * bem.C2 * x;


