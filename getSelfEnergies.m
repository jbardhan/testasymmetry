function E_self = getSelfEnergies(pqr, bem)
E_self = zeros(length(pqr.q),1);
iterativeSolver = 1;
P = diag(diag(bem.A));

for i = 1:length(pqr.q)
  rhs = bem.B(:,i); % assume q_j = 0 except for q_i = 1
  if iterativeSolver
    sigma = gmres(bem.A,rhs,[],1e-5,100,P);
  else
    sigma = bem.A \ rhs;
  end
  E_self(i) = 0.5 * bem.C(i,:)*sigma;
end
