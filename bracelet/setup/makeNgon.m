function pqr = makeNgon(d, N)

xyz = [ 0 0 0]; q = [0]; R = [0];
deltaAngle = 2*pi/N;
currentAngle = 0;
for i=2:N
  nextXyz = xyz(end,:) + d * [cos(currentAngle) sin(currentAngle) 0];
  currentAngle = currentAngle+deltaAngle;
  nextQ =  0;
  nextR =  0;
  xyz = [xyz; nextXyz];
  q   = [q; nextQ];
  R   = [R; nextR];
end

pqr = struct('xyz',xyz,'q',q,'R',R);
