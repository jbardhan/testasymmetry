asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

origin = [0 0 0];
R = 1;
q = -1;

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
kcal_to_kJ = 4.18;
density = 4.0;

epsTilde = (epsOut-epsIn)/epsOut;

qdGdn = q * (-1/R^2);
ECoul = -qdGdn;
alpha = asymParams.alpha;
beta  = asymParams.beta;
gamma = asymParams.EfieldOffset;
mu =  -alpha * tanh(-gamma);

rhs = epsTilde * qdGdn;
sigmaBorn = epsTilde * qdGdn;

sigma = linspace(-0.2*sigmaBorn, 5*sigmaBorn, 100);
E_sigma = -0.5 * sigma;
h = (alpha*tanh(beta*E_sigma+beta*ECoul-gamma)+mu);
lhs = (1-epsTilde.*h).*sigma;