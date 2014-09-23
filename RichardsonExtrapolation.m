function finalVal = RichardsonExtrapolation(i, j, nVec, Aestimated, ...
					    k0)

n  = nVec(i);
tn = nVec(j);
t = tn/n;

An = Aestimated(i);
Atn = Aestimated(j);

finalVal = (t^(-k0)*Atn - An)/(t^(-k0)-1);