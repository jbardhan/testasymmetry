function pqrNew = translateAtoms(pqr,vec)
pqrNew = pqr;
pqrNew.xyz = pqrNew.xyz+ones(length(pqr.q),1)*vec;