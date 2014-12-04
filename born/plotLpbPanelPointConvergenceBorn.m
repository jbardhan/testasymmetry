printOn = 0;

addpath('../');
addpath('../../pointbem');
addpath('../../panelbem');
addpath('../../nonlocal');
loadConstants
global quadrule_order quadrule_x quadrule_w;
quadrule_order = 25;
[quadrule_x, quadrule_w]=setupquad(quadrule_order);

asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);


origin = [0 0 0];
q_list = [-1 1];
R_list = {'1'};%, '1.375', '1.75', '2.125', '2.5'}; 
kappa = 0.125; % kappa^-1 = 8 A
density_list = 2:2:10;
epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;

maxPicardIterations = 30;
for i=1:length(R_list)
  R= str2double(R_list{i});
  for j=1:length(q_list)
    q = q_list(j);
    E_picard_lpb(i,j) = bornPicardNoStern(R, q, epsIn, epsOut, kappa, ...
					  asymParams, conv_factor, ...
					  maxPicardIterations);
  end
end


for k=1:length(density_list)
  density = density_list(k);
  for i=1:length(R_list)
    Rstr = R_list{i};
    R =  str2double(Rstr);
    srfFile = sprintf('born_%sA_%d.srf',Rstr,density);
    numPoints(i,k) = floor(4 * pi * R^2 * density);
    srfPoint   = makeSphereSurface(origin, R, numPoints(i,k));
    
    srfPanel=loadSrfIntoPanels(srfFile);
    numPanels(i,k) = length(srfPanel.areas);
    for j=1:length(q_list)
      q = q_list(j);
      pqrData = struct('xyz',[0 0 0],'q',q,'R',0);
      bemPoint = makeBemYoonLPBMatrices(srfPoint, pqrData, epsIn, epsOut,kappa);
      [phiReacPoint, phiBndyPoint, dphiDnBndyPoint] = ...
	  solveConsistentYoonNoSternAsym(srfPoint, bemPoint, epsIn, ...
					 epsOut, conv_factor, pqrData, ...
					 asymParams);
      L_point(i,j,k) = 0.5 * pqrData.q' * phiReacPoint;

      bemPanelPcm = makePanelBemEcfQualMatrices(srfPanel, pqrData, epsIn, epsOut);
      bemPanel = makePanelBemYoonLPBMatrices(srfPanel, pqrData, ...
					     epsIn,epsOut,kappa);
      asymBemPanel = makePanelAsymEcfCollocMatrices(srfPanel, bemPanelPcm,pqrData);
      [phiReacPanel, phiBndyPanel,dphiDnBndyPanel] = ...
	  solvePanelConsistentYoonNoSternAsym(srfPanel, bemPanel, ...
					      epsIn, epsOut, conv_factor, ...
					      pqrData, asymParams, ...
					      asymBemPanel);
      L_panel(i,j,k) = 0.5 * pqrData.q' * phiReacPanel;
    end
  end
end
