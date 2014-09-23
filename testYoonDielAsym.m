addpath('../pointbem');
loadConstants
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

origin = [0 0 0];
R = 2.0;
chargeLocation = [0 0 0];
q_list = [1 -1];
epsIn  =  2;
epsOut = 80;
conv_factor = 332.112;

densities = 0.5:1:8;

for i=1:length(densities)
  density = densities(i);
  numPoints(i) = ceil(4 * pi * density * R^2);
  surfdata   = makeSphereSurface(origin, R, numPoints(i));
  
  for j=1:length(q_list)
    q = q_list(j);
    pqr = struct('xyz',chargeLocation,'q',q,'R',R);
    
    bemEcf = makeBemEcfQualMatrices(surfdata, pqr,  epsIn, epsOut);
    bemYoonDiel = makeBemYoonDielMatrices(surfdata, pqr,  epsIn, epsOut);
    [phiReacEcf, sigma] = solveConsistentAsymmetric(surfdata, bemEcf, ...
						    epsIn, epsOut, ...
						    conv_factor, pqr, asymParams);
    [phiReacYoonDiel, phiBndy,dPhiBndy] = solveConsistentYoonNoSternAsym(surfdata, ...
						  bemYoonDiel, epsIn, ...
						  epsOut, conv_factor, ...
						  pqr, asymParams);
    E_ecf(i,j) = 0.5 * q'*phiReacEcf;
    E_yoondiel(i,j) = 0.5 * q'*phiReacYoonDiel;
    fprintf('numPoints = %d, ECF = %f, YL = %f\n',numPoints(i), E_ecf(i), ...
	    E_yoondiel(i));
  end
end

if length(densities) < 5
  fprintf('Error: not enough densities to do Richardson extrapolation\n');
  return
end

expectedOrder = -0.5; % with respect to the number of points
index1 = length(densities) - 4; 
index2 = length(densities);

for j=1:length(q_list)
E_ecf_richardson_extrap(j) = RichardsonExtrapolation(index1, index2, ...
						  numPoints, E_ecf(:,j), ...
						  expectedOrder);
E_yoondiel_richardson_extrap(j) = RichardsonExtrapolation(index1, index2,...
						  numPoints, E_yoondiel(:,j),...
						  expectedOrder);
end

