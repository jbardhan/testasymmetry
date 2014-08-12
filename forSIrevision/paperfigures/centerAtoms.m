function pqrNew = centerAtoms(pqr)

meanpos = mean(pqr.xyz,1);
pqrNew = pqr;
pqrNew.xyz = pqrNew.xyz - ones(length(pqr.q),1) * meanpos;
