addpath('../pointbem');
loadConstants

origin = [0 0 0];
R      = 6.0;
epsIn  =  4;
epsOut = 80;
numCharges = -1;
conv_factor = 332.112;
waterModel = struct('tau',1,'R_OH',0.58,'rho_w',1.4,'R_s',0.52);

% set up sphere and boundary integral operators
density = 1.0;
numPoints = ceil(4 * pi * density * R^2)
surfdata   = makeSphereSurface(origin, R, numPoints);

% set up interior grid of charges
h      = 0.25;
hMin   = 1.0; % for now
[pqrReacGrid,Ygrid,Zgrid]  = addPqrGridSpherePlane(R-hMin, struct('q',[],'xyz',[],'R',[]),h); 

% now we can set up the full BEM system (charges to surface and back)
bemReacGrid = makeBemEcfQualMatrices(surfdata, pqrReacGrid, epsIn, epsOut);

% solve BEM system for individual self energies, and get Born radii
E_ReacGrid_preconvfactor = getSelfEnergies(pqrReacGrid,bemReacGrid);
R_eff_ReacGrid = getEffectiveRadii(E_ReacGrid_preconvfactor,epsIn,epsOut);

% repeat the above three steps for the sources on the line segment from the center to the
% surface along the Z axis
pqrSourceChargeLine = struct('xyz', [0 0.5 R-hMin; 0 -0.5 R-hMin], 'q', [1; ...
		    -1], 'R', [0; 0]);
bemSourceChargeLine = makeBemEcfQualMatrices(surfdata, pqrSourceChargeLine, ...
					     epsIn, epsOut);
E_SourceChargeLine_preconvfactor = getSelfEnergies(pqrSourceChargeLine,bemSourceChargeLine);
R_eff_SourceChargeLine = ...
    getEffectiveRadii(E_SourceChargeLine_preconvfactor,epsIn,epsOut);

% now compute the reaction potentials at the grid locations for
% each single charge along the line
L_total = zeros(length(pqrReacGrid.q),length(pqrSourceChargeLine.q));
for i=1:length(pqrSourceChargeLine.q)
  [pqrNew, R_eff_New] = combineSourcesSpecial(pqrReacGrid,R_eff_ReacGrid,...
					      pqrSourceChargeLine, ...
					      R_eff_SourceChargeLine, ...
					      i);
  R_scaled_New = getScaledEffectiveRadii(pqrNew, R_eff_New, waterModel);

  for j=1:length(pqrReacGrid.q)
  L_asym(j,i) = ...
      getAsymmetricGBReactionPotential(pqrSourceChargeLine.xyz(i,:),...
				       R_eff_SourceChargeLine(i),...
				       R_scaled_New(length(R_eff_ReacGrid)+1),...
				       pqrSourceChargeLine.q(i),...
				       pqrReacGrid.xyz(j,:),...
				       R_eff_ReacGrid(j),...
				       R_scaled_New(j),...
				       epsIn,epsOut);
  end
end


for i=1:length(pqrSourceChargeLine.q)
  for j=1:length(pqrReacGrid.q)
  L_sym(j,i) = getStandardGBReactionPotential(pqrSourceChargeLine.xyz(i,:), ...
					      R_eff_SourceChargeLine(i), ...
					      pqrSourceChargeLine.q(i), ...
					      pqrReacGrid.xyz(j,:), ...
					      R_eff_ReacGrid(j),epsIn,epsOut);
  end
end

phiReac_asym = conv_factor * (L_asym * [1; 1]); % the 1;1 vector is
                                                % "add both contributions"
phiReac_sym  = conv_factor * (L_sym * [1; 1]);