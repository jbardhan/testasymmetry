addpath('../pointbem');
addpath('../panelbem');
loadConstants

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

res = {'arg', 'asp', 'cys', 'glu', 'his', 'lys', 'tyr'};
oneletter = {'r','d','c','e','h','k','y'};

pstate = {'prot','deprot'};
surftype = {'_scaledcharmm_1','_scaledcharmm_2','_scaledcharmm_3', ...
	    '_scaledcharmm_4','_scaledcharmm_5','_scaledcharmm_6'};%'_scaledcharmm_7',...
%	    '_scaledcharmm_8'};

for k=1:length(surftype)
    
  for i=1:length(res)
    for j=1:length(pstate)
      srfFile = sprintf('saltresidues/%s/%s%s.srf',res{i},res{i},surftype{k});
      pdbFile = sprintf('saltresidues/%s/%s.pdb',res{i},res{i});
      crgFile = sprintf('saltresidues/%s/j%s%d.crg',res{i},oneletter{i},j);
      % note that we use j rather than index into pstate(j)!!!

      pqrData = loadPdbAndCrg(pdbFile,crgFile);
      srfData = loadSrfIntoPanels(srfFile);
      
      bemPcm = makePanelBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
      asymBemPcm = makePanelAsymEcfCollocMatrices(srfData, bemPcm, pqrData);
      [phiReac, sigma] = solvePanelConsistentAsymmetric(srfData, bemPcm, ...
							epsIn, epsOut, ...
							conv_factor, ...
							pqrData, ...
							asymParams,asymBemPcm);
      dG_pcm(i,j,k) = 0.5 * pqrData.q' * phiReac;

      bemYoonDiel = makePanelBemYoonDielMatrices(srfData,pqrData,epsIn,epsOut);
      asymBemYL   = asymBemPcm; % we're reusing asymBemPcm
      [phiReacYL, phiBndy, dphiDnBndy] = ...
	  solvePanelConsistentYoonNoSternAsym(srfData, bemYoonDiel, ...
					      epsIn, epsOut, ...
					      conv_factor, pqrData, ...
					      asymParams, asymBemYL);

      dG_yl(i,j,k) = 0.5 * pqrData.q' * phiReacYL;
      fprintf('Done w/ %s%s.srf and j%s%d.crg.\n',res{i},surftype{k},oneletter{i},j);
    end
  end
end

density_list = 1:6;
for i=1:7;%length(R_list)
  for j=1:2%length(q_list)
    E_extrap_pcm(i,j) = RichardsonExtrapolation(length(density_list)-1,length(density_list),...
						density_list, squeeze(dG_pcm(i,j,:)),-1);
    E_extrap_yl(i,j) = RichardsonExtrapolation(length(density_list)-1,length(density_list),...
						density_list, squeeze(dG_yl(i,j,:)),-1);
  end
end
