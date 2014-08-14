addpath('../pointbem');
loadConstants

% radii.siz should say 1.75 for alpha=0.5, beta=-60, Eoff=-0.5
epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
asymParams = struct('alpha',0.0, 'beta', 0.0,'EfieldOffset',- ...
		    0.0);

res = 'asp'; 
crgname = 'jd2';
surftype = '_roux_6';
srfFile = sprintf('residues/%s/%s%s.srf',res,res,surftype)
pdbFile = sprintf('residues/%s/%s.pdb',res,res);
crgFile = sprintf('residues/%s/%s.crg',res,crgname);


pqrData = loadPdbAndCrg(pdbFile,crgFile);
srfData = loadSrfIntoSurfacePoints(srfFile);
surfsurfop = makeSurfaceToSurfaceOperators(srfData);
chargesurfop = makeSurfaceToChargeOperators(srfData, pqrData);

bem = makeBemMatrices(srfData, pqrData, surfsurfop, chargesurfop, ...
		      epsIn, epsOut);

[phiReac, sigma] = solveConsistentAsymmetric(srfData, surfsurfop, chargesurfop, ...
					     bem, epsIn, epsOut, conv_factor, ...
					     pqrData, asymParams);
dG = 0.5 * pqrData.q' * phiReac
