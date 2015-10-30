function Eneg = CalcEnergy(Params,R,q,epsOut,epsIn)
%              Calculates the energy of solvation for a specific trial.

tildeEps = (epsOut - epsIn)/epsOut;

alpha = Params.alpha;
beta  = Params.beta;
gamma = Params.EfieldOffset;
mu    = -alpha * tanh(-gamma);

h = alpha.*tanh(beta.*[q./(4.*pi.*R.^2)] - gamma) + mu ;

Eneg = (1 + tildeEps.*(h));