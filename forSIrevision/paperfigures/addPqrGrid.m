function [pqrDataNew,Y,Z] = addPqrGrid(radius, pqrData, spacing)

thSweep = (0:8:360)*pi/180;
radSweep = 0:spacing:radius;
[t,r]=meshgrid(thSweep, radSweep);
[Y,Z]=pol2cart(t,r);
yPoints = reshape(Y, length(thSweep)*length(radSweep),1);
zPoints = reshape(Z, length(thSweep)*length(radSweep),1);
numPoints = length(yPoints);
qNew = zeros(numPoints, 1);
xyzNew = [zeros(numPoints,1) yPoints zPoints];
pqrDataNew = struct('q', [pqrData.q; qNew], 'xyz', [pqrData.xyz; xyzNew]);
