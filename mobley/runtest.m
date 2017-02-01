% Path information
Home = getenv('HOME');
addpath(sprintf('%s/repos/pointbem',Home));
addpath(sprintf('%s/repos/panelbem',Home));
addpath(sprintf('%s/repos/testasymmetry',Home));

loadConstants

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
staticpotential = 10.7; 
symParams = struct('alpha',0.0, 'beta', -60.0,'EfieldOffset',-0.5);
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

surftype = {'4'};
kappa_list = 0.0;

for k=1:length(surftype)
  pqrData = readpqr('test.pqr')

  srfSternFile   = sprintf('test_%s.srf', surftype{k});
  srfSternData = loadSternSrfIntoPanels(srfSternFile);

  bemPcm = makePanelBemEcfQualMatrices(srfSternData.dielBndy(1), pqrData, epsIn, epsOut);
  asymBemPcm = makePanelAsymEcfCollocMatrices(srfSternData.dielBndy(1), bemPcm, pqrData);
  fprintf('Finished setting up %d/%d asymBemPcm\n',...
	  k,length(surftype));

  for l=1:length(kappa_list)
    kappa = kappa_list(l)
    
    bemYoonStern = makePanelBemSternMatrices(srfSternData, ...
					     pqrData, epsIn, ...
					     epsOut, kappa);
    fprintf('Done setting up bemYoonStern.\n');
    % salt with Stern, symmetric
    [phiReacSym, phiBndy, dphiDnBndy] = ...
	solvePanelConsistentSternAsym(srfSternData.dielBndy(1), ...
				      srfSternData.sternBndy(1), ...
				      pqrData, bemYoonStern, ...
				      epsIn, epsOut, kappa, ...
				      conv_factor, symParams, ...
				      asymBemPcm);
    dG_sym(k,l) = 0.5 * pqrData.q' * phiReacSym;
    withstatic_sym(k,l) = dG_sym(k,l) + staticpotential*sum(pqrData.q);
    
    % salt with Stern, NLBC 
    % again, can re-use asymBemYL
    [phiReac, phiBndy, dphiDnBndy] = ...
	solvePanelConsistentSternAsym(srfSternData.dielBndy(1), ...
				      srfSternData.sternBndy(1), ...
				      pqrData, bemYoonStern, ...
				      epsIn, epsOut, kappa, ...
				      conv_factor, asymParams, ...
				      asymBemPcm);
    dG_asym(k,l) = 0.5 * pqrData.q' * phiReac;
    withstatic_asym(k,l) = dG_asym(k,l) + staticpotential*sum(pqrData.q);
  end

end

fprintf('nlbc out %8.3f\n', withstatic_asym(1));
quit