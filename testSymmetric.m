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
bemReacGrid = makeBemEcfQualMatrices(surfdata, pqrReacGrid,  epsIn, epsOut);

% solve BEM system for individual self energies, and get Born radii
E_ReacGrid_preconvfactor = getSelfEnergies(pqrReacGrid,bemReacGrid);
R_eff_ReacGrid = getEffectiveRadii(E_ReacGrid_preconvfactor,epsIn,epsOut);

% repeat the above three steps for the sources on the line segment from the center to the
% surface along the Z axis
pqrSourceChargeLine = makeSourceChargeLine(R, hFine, h);
bemSourceChargeLine = makeBemEcfQualMatrices(surfdata, pqrSourceChargeLine, ...
					     epsIn, epsOut);
E_SourceChargeLine_preconvfactor = getSelfEnergies(pqrSourceChargeLine,bemSourceChargeLine);
R_eff_SourceChargeLine = ...
    getEffectiveRadii(E_SourceChargeLine_preconvfactor,epsIn,epsOut);

% now compute the reaction potentials at the grid locations for
% each single charge along the line
L_total = zeros(length(pqrReacGrid.q),length(pqrSourceChargeLine.q));
for i=1:length(pqrSourceChargeLine.q)
  for j=1:length(pqrReacGrid.q)
  L_total(j,i) = ...
      getStandardGBReactionPotential(pqrSourceChargeLine.xyz(i,:), ...
				     R_eff_SourceChargeLine(i), 1, ...
				     pqrReacGrid.xyz(j,:), R_eff_ReacGrid(j),epsIn,epsOut);
  end
end
