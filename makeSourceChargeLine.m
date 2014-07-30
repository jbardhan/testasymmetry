function pqrData = makeSourceChargeLine(radius, hLine, hToSurface)

Z = (0:hLine:radius-hToSurface)';
if length(Z) < 1
  fprintf('Error: no points for R=%f, hLine=%f, hToSurface=%f\n',radius,hLine,hToSurface);
  pqrData = 0;
  return;
end
q = 0 * Z;
xyz = [0*Z 0*Z Z];
pqrData = struct('q',q,'xyz',xyz,'R',0*q);