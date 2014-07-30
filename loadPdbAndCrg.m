function pqr = loadPdbAndCrg(pdbFile, crgFile)
NmaxAtoms = 10000;
[crd, crg, atype, resname ] = readPdb(pdbFile, NmaxAtoms);
crg = readCrg(crgFile);
pqr = struct('xyz', crd, 'q', crg, 'R', 0 * crg);
