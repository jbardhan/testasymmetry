addpath('../pointbem');
loadConstants

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
asymParams = struct('alpha',0.0, 'beta', -60.0,'EfieldOffset',-0.5);

res = {'arg'};%, 'asp', 'cys', 'glu', 'his', 'lys', 'tyr'};
oneletter = {'r','d','c','e','h','k','y'};
pstate = {'prot','deprot'};  %prot always 1, deprot always 2-N
surftype = {'_scaledcharmm_1','_scaledcharmm_2'};%,'_scaledcharmm_3','_scaledcharmm_4','_scaledcharmm_5','_scaledcharmm_6'};%, '_scaledcharmm_7','_scaledcharmm_8'};  % 


for k=1:length(surftype)
    
  for i=1:length(res)
    for j=1:length(pstate)
      srfFile = sprintf('saltresidues/%s/%s%s.srf',res{i},res{i},surftype{k})
      pdbFile = sprintf('saltresidues/%s/%s.pdb',res{i},res{i});
      crgFile = sprintf('saltresidues/%s/j%s%d.crg',res{i},oneletter{i},j)
      % note that we use j rather than index into pstate(j)!!!

      pqrData = loadPdbAndCrg(pdbFile,crgFile);
      srfData = loadSrfIntoSurfacePoints(srfFile);
      
      bem = makeBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
      [phiReac, sigma] = solveConsistentAsymmetric(srfData, bem, ...
						   epsIn, epsOut, ...
						   conv_factor, ...
						   pqrData, asymParams);
      bemYoonDiel = makeBemYoonDielMatrices(srfData, pqrData, epsIn, ...
					    epsOut);
      [phiReacYoonDiel, phiBndy,dphiBndy] = ...
	  solveConsistentYoonNoSternAsym(srfData, bemYoonDiel, ...
					 epsIn, epsOut, conv_factor, ...
					 pqrData, asymParams);

      dG_pcm(i,j) = 0.5 * pqrData.q' * phiReac;
      dG_yl(i,j) = 0.5 * pqrData.q' * phiReacYoonDiel;
      
    end
  end
  asym_pcm(:,:,k) = dG_pcm;
  asym_yl(:,:,k) = dG_yl;
end

