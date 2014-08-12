function vfstruct = readmesh(base,noHeader)

vertfile = sprintf('%s.vert', base);
facefile = sprintf('%s.face', base);

if nargin < 2 
  [Xv,Yv,Zv,nx,ny,nz,j1,j2,j3] = textread(vertfile,...
												 '%f %f %f %f %f %f %d %d %d',...
												 'headerlines',3);

  [v1,v2,v3,j1,j2] = textread(facefile,...
										'%d %d %d %d %d',...
										'headerlines',3);
else
  [Xv,Yv,Zv,nx,ny,nz,j1,j2,j3] = textread(vertfile,...
												 '%f %f %f %f %f %f %d %d %d');

  [v1,v2,v3,j1,j2] = textread(facefile,...
										'%d %d %d %d %d');
end

v = [Xv Yv Zv];
f = [v1 v3 v2];
for i=1:size(f, 1)
  X(:,i) = v(f(i,1:3),1);
  Y(:,i) = v(f(i,1:3),2);
  Z(:,i) = v(f(i,1:3),3);
end

vfstruct = struct('vert', v, 'face', f, 'X', X, 'Y', Y, 'Z', Z);
