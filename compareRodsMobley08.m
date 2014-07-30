addpath('../pointbem');
loadConstants

asymParams = struct('alpha', 0.6, 'beta', -20.0,'EfieldOffset',0,'deltaOffset',0);

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;

vertexDensity = 6;

Prods = [
1	-18.301472
2	-59.386276
3	-87.306433
4	-14.972744
5	-17.830256
6	-58.915060
7	-87.595533
8	-105.389039];

Nrods = [
1	-21.591817
2	-61.917939
3	-90.850761
4	-14.718761
5	-21.374584
6	-61.700706
7	-90.380362
8	-108.680201];

Arods = [
1	3.327128
2	2.654256
3	3.500477
4	-0.248935
5	3.425414
6	2.666508
7	3.407631
8	3.298420];

var1_length = [5 5 5 5 6 6 6 6];
var2_type   = [12 13 14 23 12 13 14 15 ];
var3_posOrNeg = {'N','P'};

for i=1:length(var1_length) % simultaneously looping over var2_type
  for j=1:length(var3_posOrNeg)
    srfFile = sprintf('./rod/l%d_q%d/L%d_%d.srf',var1_length(i), ...
		      var2_type(i), var1_length(i), vertexDensity);
    pdbFile = sprintf('./rod/l%d_q%d/L%d.pdb',var1_length(i), ...
		      var2_type(i), var1_length(i));
    crgFile = sprintf('./rod/l%d_q%d/%s%d.crg', var1_length(i), ...
		      var2_type(i), var3_posOrNeg{j},var1_length(i));

    pqrData = loadPdbAndCrg(pdbFile,crgFile);
    srfData = loadSrfIntoSurfacePoints(srfFile);

    surfsurfop = makeSurfaceToSurfaceOperators(srfData);
    chargesurfop = makeSurfaceToChargeOperators(srfData, pqrData);

    bem = makeBemMatrices(srfData, pqrData, surfsurfop, chargesurfop, ...
			  epsIn, epsOut);
    
    [phiReac, sigma] = solveAsymmetric(srfData, surfsurfop, chargesurfop, ...
				       bem, epsIn, epsOut, conv_factor, ...
				       pqrData, asymParams);
    
    rodBEM(i,j) = 0.5 * pqrData.q' * phiReac;

  end
end

ProdBEM = rodBEM(:,2);
NrodBEM = rodBEM(:,1);
ArodBEM = ProdBEM-NrodBEM;

