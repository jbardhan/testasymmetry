addpath('../pointbem');
loadConstants

asymParams = struct('alpha', 0.6, 'beta', -20.0,'EfieldOffset',0,'deltaOffset',0);

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
waterModel = struct('tau',1,'R_OH',0.58,'rho_w',1.4,'R_s',0.52); 

vertexDensity = 6;
var1_numSides = 3:8;
var2_posOrNeg = {'N','P'};
var3_typeOfChargeDist = {'dipole','distrib','oppose'};

testQ = var3_typeOfChargeDist{1};

for i=1:length(var1_numSides)
  for j=1:length(var2_posOrNeg)
    srfFile = sprintf('./bracelet/n%d/n%d_%d.srf',var1_numSides(i), ...
		      var1_numSides(i), vertexDensity);
    pdbFile = sprintf('./bracelet/n%d/N%d.pdb',var1_numSides(i), var1_numSides(i));
    crgFile = sprintf('./bracelet/n%d/%s%d_%s.crg', var1_numSides(i), ...
		      var2_posOrNeg{j}, var1_numSides(i), testQ)

    pqrData = loadPdbAndCrg(pdbFile,crgFile);
    srfData = loadSrfIntoSurfacePoints(srfFile);

    bem = makeBemEcfQualMatrices(srfData, pqrData, ...
			  epsIn, epsOut);
    
    [phiReac, sigma] = solveAsymmetric(srfData, bem, epsIn, epsOut, ...
				       conv_factor, pqrData, asymParams);
    
    deltaG(i,j) = 0.5 * pqrData.q' * phiReac;

  end
end
