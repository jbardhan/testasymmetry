function [V,K] = colloc_Laplace(meshData,centroids,normals,areas)

np = size(meshData.face,1);
for i=1:np
  numverts = 3;
  panel = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];

  [area,collocpt,Z,fss,fds,fess] = calcp(panel,centroids,normals);
  V(:,i) = fss'/4/pi;
  K(:,i) = fds'/4/pi;
end
