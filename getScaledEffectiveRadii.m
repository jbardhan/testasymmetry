function R_list = getScaledEffectiveRadii(pqr, R_eff, water)
scale_factors = 0 * R_eff;

% this can be vectorized
for i=1:length(R_eff)
  weightedSumQ = 0;
  for j=1:length(R_eff)
    r_ij = norm(pqr.xyz(i,:) - pqr.xyz(j,:));
    if i==j
      weightedSumQ = weightedSumQ + pqr.q(j);
    else
      weightedSumQ = weightedSumQ + pqr.q(j)*exp(-(water.tau*r_ij^2)/(R_eff(i)*R_eff(j)));
    end
  end
  
  scale_factors(i) = 1 + sign(weightedSumQ) * (water.R_OH / (R_eff(i) - water.R_s + water.rho_w));
end
R_list = R_eff .* scale_factors;
