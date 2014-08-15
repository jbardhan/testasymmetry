addpath('../pointbem');
loadConstants

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;

modeltype = {'sym', 'asym'};
res = {'arg', 'asp', 'cys', 'glu', 'his', 'lys', 'tyr'};
oneletter = {'r','d','c','e','h','k','y'};
pstate = {'prot','deprot'};  %prot always 1, deprot always 2-N
surftype = {'_roux_6','_scaledcharmm_6'};  % s/roux/parse/ for that set
% 1st is for std symmetric model, standard Roux radii
% 2nd is for NLBC: use uniform 0.92 scalefact like Mobley examples


for k=1:1%length(modeltype)
  if k == 1
    asymParams  = struct('alpha', 0.0, 'beta', 0.0, 'EfieldOffset', ...
			0.0);
  else 
    asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);
  end

  for i=1:length(res)
    for j=1:length(pstate)
      srfFile = sprintf('residues/%s/%s%s.srf',res{i},res{i},surftype{k})
      pdbFile = sprintf('residues/%s/%s.pdb',res{i},res{i});
      crgFile = sprintf('residues/%s/j%s%d.crg',res{i},oneletter{i},j)
      % note that we use j rather than index into pstate(j)!!!

      pqrData = loadPdbAndCrg(pdbFile,crgFile);
      srfData = loadSrfIntoSurfacePoints(srfFile);
      surfsurfop = makeSurfaceToSurfaceOperators(srfData);
      chargesurfop = makeSurfaceToChargeOperators(srfData, pqrData);
      
      bem = makeBemMatrices(srfData, pqrData, surfsurfop, chargesurfop, ...
			    epsIn, epsOut);
      
      [phiReac, sigma] = solveConsistentAsymmetric(srfData, surfsurfop, chargesurfop, ...
						   bem, epsIn, epsOut, conv_factor, ...
						   pqrData, asymParams);
      dG(i,j) = 0.5 * pqrData.q' * phiReac;
    end
  end
  if k == 1
    sym = dG;
  else
    asym = dG;
  end
  clear dG;
end

