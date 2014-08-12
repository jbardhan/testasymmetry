function radius = getPointSetRadius(pqr)

pqrNew = centerAtoms(pqr);
radius = 0;
for i=1:length(pqr.q)
  radius = max(radius,norm(pqrNew.xyz(i,:)));
end
