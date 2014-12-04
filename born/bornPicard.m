addpath('../../pointbem');
loadConstants
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

origin = [0 0 0];
if 0
  R_list = linspace(1,2.5,5);
else
  R_list = rscale * [sodiumRminOver2 chlorideRminOver2 potassiumRminOver2 ...
		    rubidiumRminOver2 magnesiumRminOver2 cesiumRminOver2 ...
		    calciumRminOver2 bariumRminOver2 zincRminOver2 cadmiumRminOver2];
end

q_list = linspace(-1,1,2);

epsIn  =  1;
epsOut = 80;
tildeEps = (epsOut - epsIn)/epsOut;
conv_factor = 332.112;

numPicardIterations = 20;

% shorthand aliases.
alpha = asymParams.alpha;
beta  = asymParams.beta;
gamma = asymParams.EfieldOffset;
mu    = -alpha * tanh(-gamma);
% initialize solutions
sigma = zeros(length(R_list), length(q_list), numPicardIterations);
a     = 0 * sigma; % surface varphi
b     = 0 * sigma; % surface dvarphidn
for i=1:length(R_list)
  R = R_list(i);
  for j=1:length(q_list)
    cur_q = q_list(j);
    dphidnCoul = -(cur_q/4/pi/R^2);
    sigma(i,j,1) = 0;

    a(i,j,1) = (epsIn/epsOut) * cur_q/4/pi/R;
    b(i,j,1) = -cur_q/4/pi/R^2;
    b_yl = [cur_q/4/pi/R; 0];
    for k=1:numPicardIterations-1
      tanharg(i,j,k) = beta*(-dphidnCoul-(-1/2)*sigma(i,j,k))-gamma;
      sigma(i,j,k+1) = (tildeEps * dphidnCoul)/(1-tildeEps*(alpha* ...
							tanh(tanharg(i,j,k))+mu));

      fakesigma = b(i,j,k)+a(i,j,k)/R;
      newarg(i,j,k) = beta*(-dphidnCoul -(-1/2)*fakesigma) - gamma;
      fk = (epsIn/(epsOut-epsIn)) - (alpha * tanh(newarg(i,j,k)) + mu);
      A = [0 -R; 1 (fk*R)/(1+fk)];
      xk = A\b_yl;
      a(i,j,k+1) = xk(1);
      b(i,j,k+1) = xk(2);
    end
    E(i,j) = .5 * conv_factor *4*pi* sigma(i,j,end) * R * cur_q / epsIn;

    E_yl(i,j) = .5 * conv_factor *4*pi * (a(i,j,end) + R * b(i,j,end))*cur_q ...
	/ epsIn;
  end
end

