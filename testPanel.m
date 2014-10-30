addpath('../pointbem');
loadConstants

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

res = {'arg'};
pstate = {'prot','deprot'};
oneletter = {'r','d','c','e','h','k','y'};
surftype = {'_scaledcharmm_1'};

for k=1:length(surftype)
    
  for i=1:length(res)
    for j=1:length(pstate)
      srfFile = sprintf('saltresidues/%s/%s%s.srf',res{i},res{i},surftype{k})
      pdbFile = sprintf('saltresidues/%s/%s.pdb',res{i},res{i});
      crgFile = sprintf('saltresidues/%s/j%s%d.crg',res{i},oneletter{i},j)
      % note that we use j rather than index into pstate(j)!!!

      pqrData = loadPdbAndCrg(pdbFile,crgFile);
      srfData = loadSrfIntoPanels(srfFile);
      
      bem = makePanelBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
      asymBem = makePanelAsymQuadMatrices(srfData, bem, pqrData,1);
      [phiReac, sigma] = solvePanelConsistentAsymmetric(srfData, bem, ...
						   epsIn, epsOut, ...
						   conv_factor, ...
						   pqrData, asymParams,asymBem);
      dG_pcm(i,j) = 0.5 * pqrData.q' * phiReac;
    end
  end
end
