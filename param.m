addpath('../pointbem');
loadConstants

asymParams = struct('alpha', 0.5, 'beta', -60.0,'EfieldOffset',-0.5);

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;

vertexDensity = 6;
var1_numSides = [3 8]; %3:8;
var2_posOrNeg = {'N','P'};
var3_typeOfChargeDist = {'dipole','distrib','oppose'};

Pdistrib = [
3	-14.448498
4	-14.322179
5	-15.635084
6	-17.287700
7	-19.025060
8	-20.592930];

Ndistrib = [
3	-19.957700
4	-21.017095
5	-23.685932
6	-26.948716
7	-29.872518
8	-32.880338];

Adistrib = [
3	5.149787
4	6.641088
5	8.062742
6	9.519200
7	10.766596
8	12.153367];

testQ = var3_typeOfChargeDist{2}; % distrib

for i=1:length(var1_numSides)
  for j=1:length(var2_posOrNeg)
    srfFile = sprintf('./bracelet/n%d/n%d_%d.srf',var1_numSides(i), ...
		      var1_numSides(i), vertexDensity);
    pdbFile = sprintf('./bracelet/n%d/N%d.pdb',var1_numSides(i), var1_numSides(i));
    crgFile = sprintf('./bracelet/n%d/%s%d_%s.crg', var1_numSides(i), ...
		      var2_posOrNeg{j}, var1_numSides(i), testQ)

    pqrData = loadPdbAndCrg(pdbFile,crgFile);
    srfData = loadSrfIntoSurfacePoints(srfFile);

    bem = makeBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
    
    [phiReac, sigma] = solveConsistentAsymmetric(srfData, bem, epsIn, ...
						 epsOut, conv_factor, ...
						 pqrData, asymParams);
    
    distribBEM(i,j) = 0.5 * pqrData.q' * phiReac;

  end
end

PdistribBEM = distribBEM(:,2);
NdistribBEM = distribBEM(:,1);
AdistribBEM = PdistribBEM - NdistribBEM;

Poppose = [
3	-14.637305
4	-8.646444
5	-36.000172
6	-29.255987
7	-53.689315
8	-46.192020];

Noppose = [
3	-19.911872
4	-15.051752
5	-42.781820
6	-40.087913
7	-59.811830
8	-57.306093];

Aoppose = [
3	5.185615
4	6.090487
5	6.612529
6	10.440835
7	5.916473
8	10.858469];

testQ = var3_typeOfChargeDist{3}; % oppose


for i=1:length(var1_numSides)
  for j=1:length(var2_posOrNeg)
    srfFile = sprintf('./bracelet/n%d/n%d_%d.srf',var1_numSides(i), ...
		      var1_numSides(i), vertexDensity);
    pdbFile = sprintf('./bracelet/n%d/N%d.pdb',var1_numSides(i), var1_numSides(i));
    crgFile = sprintf('./bracelet/n%d/%s%d_%s.crg', var1_numSides(i), ...
		      var2_posOrNeg{j}, var1_numSides(i), testQ)

    pqrData = loadPdbAndCrg(pdbFile,crgFile);
    srfData = loadSrfIntoSurfacePoints(srfFile);

    bem = makeBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
    
    [phiReac, sigma] = solveAsymmetric(srfData, bem, epsIn, epsOut, ...
				       conv_factor, pqrData, asymParams);
    
    opposeBEM(i,j) = 0.5 * pqrData.q' * phiReac;

  end
end

PopposeBEM = opposeBEM(:,2);
NopposeBEM = opposeBEM(:,1);
AopposeBEM = PopposeBEM - NopposeBEM;

Pdipole = [
3	-14.590809
4	-11.454522
5	-10.193053
6	-8.675998
7	-7.670310
8	-6.579393];

Ndipole = [
3	-19.534023
4	-18.357976
5	-17.778427
6	-16.687316
7	-15.766953
8	-14.420353];

Adipole = [
3	4.878557
4	6.708016
5	7.451335
6	7.563325
7	7.850711
8	7.472363];


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

    bem = makeBemEcfQualMatrices(srfData, pqrData, epsIn, epsOut);
    
    [phiReac, sigma] = solveAsymmetric(srfData, bem, epsIn, epsOut, ...
				       conv_factor, pqrData, asymParams);
    
    dipoleBEM(i,j) = 0.5 * pqrData.q' * phiReac;

  end
end

PdipoleBEM = dipoleBEM(:,2);
NdipoleBEM = dipoleBEM(:,1);
AdipoleBEM = PdipoleBEM - NdipoleBEM;
