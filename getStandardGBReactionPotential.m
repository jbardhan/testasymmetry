function phi = getStandardGBReactionPotential(xyzSrc, RSrc, qSrc, ...
						      xyzDest,RDest,epsIn,epsOut)

r_ij = norm(xyzSrc-xyzDest);
f_GB = sqrt(r_ij^2 + RSrc*RDest*exp(-r_ij^2/(4*RSrc*RDest)));
phi = -0.5 * (1/epsIn - 1/epsOut) * qSrc / f_GB;