addpath('../pointbem');
addpath('../panelbem');
loadConstants

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

res = {'arg', 'asp', 'cys', 'glu', 'his', 'lys', 'tyr'};
oneletter = {'r','d','c','e','h','k','y'};

pstate = {'prot'};%,'deprot'};
surftype = {'2'}%,'2','4'};
kappa_list = [0.00];% 0.500];

for k=1:length(surftype)
    
  for i=1:1 % ASP ONLY
    for j=1:length(pstate)
      % note that for crgFile we use j rather than index into pstate(j)!!!
      pdbFile = sprintf('saltresidues/%s/%s.pdb',res{i},res{i});
      crgFile = sprintf('saltresidues/%s/j%s%d.crg',res{i},oneletter{i},j);
      pqrData = loadPdbAndCrg(pdbFile,crgFile);

      srfNoSternFile = sprintf('saltresidues/%s/%s_scaledcharmm_%s.srf',res{i},res{i},surftype{k});
      srfData = loadSrfIntoPanels(srfNoSternFile);

      srfSternFile   = sprintf('saltresidues/%s/%s_scaledcharmm_stern_%s.srf',res{i},res{i},surftype{k});
      srfSternData = loadSternSrfIntoPanels(srfSternFile);

      % part 1: no-salt reference calculation: PCM/NLBC formulation
      bemPcm = makePanelBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
      asymBemPcm = makePanelAsymEcfCollocMatrices(srfData, bemPcm, pqrData);
      [phiReac, sigma] = solvePanelConsistentAsymmetric(srfData, bemPcm, ...
							epsIn, epsOut, ...
							conv_factor, ...
							pqrData, ...
							asymParams,asymBemPcm);
      dG_nosalt_pcm(i,j,k) = 0.5 * pqrData.q' * phiReac;

      % part 2: no-salt reference calculation: YL/NLBC formulation
      bemYoonDiel = makePanelBemYoonDielMatrices(srfData,pqrData,epsIn,epsOut);
      asymBemYL   = asymBemPcm; % we're reusing asymBemPcm
      [phiReacYL, phiBndy, dphiDnBndy] = ...
	  solvePanelConsistentYoonNoSternAsym(srfData, bemYoonDiel, ...
					      epsIn, epsOut, ...
					      conv_factor, pqrData, ...
					      asymParams, asymBemYL);
      dG_nosalt_yl(i,j,k) = 0.5 * pqrData.q' * phiReacYL;

      for l=1:length(kappa_list)
	kappa = kappa_list(l)

	% part 4: salt with Stern, YL/NLBC formulation
	bemYoonStern = makePanelBemSternMatrices(srfSternData, ...
						 pqrData, epsIn, ...
						 epsOut, kappa);
	% again, can re-use asymBemYL
	[phiReacStern, phiBndy, dphiDnBndy] = ...
	    solvePanelConsistentSternAsym(srfSternData.dielBndy(1), ...
					  srfSternData.sternBndy(1), ...
					  pqrData, bemYoonStern, ...
					  epsIn, epsOut, kappa, ...
					  conv_factor, asymParams, ...
					  asymBemYL);
	dG_stern(i,j,k,l) = 0.5 * pqrData.q' * phiReacStern;
	
      end


      fprintf('Done w/ %s%s.srf and j%s%d.crg.\n',res{i},surftype{k},oneletter{i},j);
    end
  end
end

qBasic = pqrData.q;
for i=1:length(pqrData.q)
  pqrData.q = 0*pqrData.q; pqrData.q(i) = qBasic(i);
  [phiReacStern, phiBndy, dphiDnBndy] = ...
      solvePanelConsistentSternAsym(srfSternData.dielBndy(1), ...
				    srfSternData.sternBndy(1), ...
				    pqrData, bemYoonStern, ...
				    epsIn, epsOut, kappa, ...
				    conv_factor, asymParams, ...
				    asymBemYL);
  selfEnergyRoux(i) = 0.5 * pqrData.q' * phiReacStern
end
