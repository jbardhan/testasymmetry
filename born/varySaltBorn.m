addpath('../../pointbem');
addpath('../../panelbem');
addpath('..');
loadConstants

origin = [0 0 0];
R_list = 0.92 * linspace(1,2.5,10);% [1.41075 2.27];
% sodium RminOver2 1.41075 from newest Roux paper
% chloride RminOver2 2.27 

sternLayerThickness = 5.0;
chargeLocation = [0 0 0];
q_list = [1 -1];

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

surftype = {'1','2'};
kappa_list = [0.000];% 0.02 0.125 0.200 0.500 1.0];

k = 1; 
densityDiel = 4;
densityStern = 1;

for i=1:length(R_list)
  R = R_list(i);
  Rstern = R + sternLayerThickness;

  numDielPoints(i) = ceil(4 * pi *  densityDiel * R^2);
  numSternPoints(i) = ceil(4 * pi * densityStern * Rstern^2);
  numTotalPoints(i) = numDielPoints(i) + numSternPoints(i);
  
  dielSurfData = makeSphereSurface(origin, R, numDielPoints(i));
  sternSurfData = makeSphereSurface(origin, Rstern, numSternPoints(i));

  for j=1:length(q_list)
    q = q_list(j);
    pqrData = struct('xyz',chargeLocation,'q',q,'R',R);

    for k=1:length(kappa_list)
      kappa = kappa_list(k);


      bemEcfAsym      = makeBemEcfQualMatrices(dielSurfData, pqrData,  epsIn, epsOut);
      bemStern        = makeBemSternMatrices(dielSurfData, sternSurfData, pqrData, ...
				   epsIn, epsOut, kappa);
      phiReac = solveConsistentSternAsym(dielSurfData, sternSurfData, ...
					 pqrData, bemStern, epsIn, epsOut,...
					 kappa, conv_factor, ...
					 asymParams);
      dG_stern(i,j,k) = 0.5 * pqrData.q' * phiReac;

    end
  end
end