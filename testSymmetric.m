addpath('../pointbem');
loadConstants

origin = [0 0 0];
R      = 6.0;
epsIn  =  4;
epsOut = 80;
numCharges = -1;
conv_factor = 332.112;

% set up sphere and boundary integral operators
density = 1.0;
numPoints = ceil(4 * pi * density * R^2)
surfdata   = makeSphereSurface(origin, R, numPoints);

% set up interior grid of charges
h       = 1.0;
hFine   = h; % for now
pqrReacGrid  = makeSphereChargeDistribution(R, numCharges, h); 

% now we can set up the full BEM system (charges to surface and back)
bemReacGridECF = makeBemEcfQualMatrices(surfdata, pqrReacGrid,  epsIn, epsOut);
bemReacGridYoonDiel = makeBemYoonDielMatrices(surfdata, pqrReacGrid, ...
					      epsIn, epsOut);

% repeat the above three steps for the sources on the line segment from the center to the
% surface along the Z axis
pqrSourceChargeLine = makeSourceChargeLine(R, hFine, h);
bemEcfSourceChargeLine = makeBemEcfQualMatrices(surfdata, pqrSourceChargeLine, ...
					     epsIn, epsOut);
bemYoonDielSourceChargeLine = makeBemYoonDielMatrices(surfdata, ...
						  pqrSourceChargeLine, ...
						  epsIn, epsOut);


