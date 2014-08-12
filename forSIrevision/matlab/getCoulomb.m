function E = getCoulomb(pqrData,E_omega)
global Na q joulesPerCalorie E_0
nq = length(pqrData.q);

Cmat = zeros(nq,nq);
for i=1:nq
  for j=1:i-1
	 Cmat(i,j) = 1/norm(pqrData.xyz(i,:)-pqrData.xyz(j,:))/4/pi/E_omega/E_0;
  end
end
Cmat = Cmat+Cmat';
E = 1/2 * (Na/1000) * q^2 * 1e10 ...
	 * pqrData.q' * Cmat * pqrData.q;
