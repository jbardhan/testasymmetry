function E_yl = bornPicardNoStern2(R, q, epsIn, epsOut, kappa, asymParams, ...
				  conv_factor, numPicardIterations)           
              
tildeEps = (epsOut - epsIn)/epsOut;

% shorthand aliases.
alpha = asymParams.alpha;
beta  = asymParams.beta;
gamma = asymParams.EfieldOffset;
mu    = -alpha * tanh(-gamma);

% initialize solutions
sigma = zeros(1,1, numPicardIterations);
a     = 0 * sigma; % surface varphi
b     = 0 * sigma; % surface dvarphidn
i = 1:32;
j = 1;
dphidnCoul = -(q./4/pi./R.^2);
sigma(i,j,1) = 0;

a(i,j,1) = (epsIn/epsOut) .* q./4/pi./R;
b(i,j,1) = -q./4/pi./R.^2;
b_yl = [q./4/pi./R; zeros(length(dphidnCoul))];

[lambdaV,lambdaK] = computeYukawaEigenvalues(R, kappa, 10);
lambdaVY0 = lambdaV(1);
lambdaKY0 = lambdaK(1);

for k=1:numPicardIterations-1
  tanharg(i,j,k) = beta*(-dphidnCoul-(-1/2)*sigma(i,j,k))-gamma;
  sigma(i,j,k+1) = (tildeEps * dphidnCoul)/(1-tildeEps*(alpha* ...
						  tanh(tanharg(i,j,k))+mu));
  
  fakesigma = b(i,j,k)+a(i,j,k)/R;
  newarg(i,j,k) = beta*(-dphidnCoul -(-1/2)*fakesigma) - gamma;
  fk = (epsIn/(epsOut-epsIn)) - (alpha * tanh(newarg(i,j,k)) + mu);
  A = [0 -R; 0.5-lambdaKY0 (fk*lambdaVY0)/(1+fk)];
  xk = A\b_yl;
  a(i,j,k+1) = xk(1);
  b(i,j,k+1) = xk(2);
end
E_pcm = .5 * conv_factor *4*pi* sigma(i,j,end) * R * q / epsIn;

E_yl = .5 * conv_factor *4*pi * (a(i,j,end) + R * b(i,j,end))*q / epsIn;

