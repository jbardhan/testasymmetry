function R_list = getEffectiveRadii(E_list, epsIn, epsOut)

R_list = -0.5 * (1/epsIn - 1/epsOut) ./ E_list;
