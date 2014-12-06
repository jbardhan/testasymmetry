addpath('../pointbem');
addpath('../panelbem');
loadConstants

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
staticpotential = 10.;
symParams = struct('alpha',0.0, 'beta', -60.0,'EfieldOffset',-0.5);
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

res = {'arg', 'asp', 'cys', 'glu', 'his', 'lys', 'tyr'};
oneletter = {'r','d','c','e','h','k','y'};

pstate = {'prot','deprot'};
surftype = {'1','2','4'};
kappa_list = [0.00 0.125 0.5 1.0];

for k=1:length(surftype)
    
  for i=1:7
    for j=1:length(pstate)
      % note that for crgFile we use j rather than index into pstate(j)!!!
      pdbFile = sprintf('saltresidues/%s/%s.pdb',res{i},res{i});
      crgFile = sprintf('saltresidues/%s/j%s%d.crg',res{i},oneletter{i},j);
      pqrData = loadPdbAndCrg(pdbFile,crgFile);

      srfSternFile   = sprintf('saltresidues/%s/%s_scaledcharmm_stern_%s.srf',res{i},res{i},surftype{k});
      srfSternData = loadSternSrfIntoPanels(srfSternFile);

      bemPcm = makePanelBemEcfQualMatrices(srfSternData.dielBndy(1), pqrData, epsIn, epsOut);
      asymBemPcm = makePanelAsymEcfCollocMatrices(srfSternData.dielBndy(1), bemPcm, pqrData);
      fprintf('Finished setting up %d/%d, %d/%d, %d/%d asymBemPcm\n',...
	      k,length(surftype),i,7,j,length(pstate));
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
	dG_sym(i,j,k,l) = 0.5 * pqrData.q' * phiReacSym;
	withstatic_sym(i,j,k,l) = dG_sym(i,j,k,l) + staticpotential*sum(pqrData.q);

	% salt with Stern, NLBC 
	% again, can re-use asymBemYL
	[phiReac, phiBndy, dphiDnBndy] = ...
	    solvePanelConsistentSternAsym(srfSternData.dielBndy(1), ...
					  srfSternData.sternBndy(1), ...
					  pqrData, bemYoonStern, ...
					  epsIn, epsOut, kappa, ...
					  conv_factor, asymParams, ...
					  asymBemPcm);
	dG_asym(i,j,k,l) = 0.5 * pqrData.q' * phiReac;
	withstatic_asym(i,j,k,l) = dG_asym(i,j,k,l) + staticpotential*sum(pqrData.q);
      end


      fprintf('Done w/ %s%s.srf and j%s%d.crg.\n',res{i},surftype{k},oneletter{i},j);
    end
  end
  clear bemYoonStern asymBemPcm bemPcm;
  save
end

