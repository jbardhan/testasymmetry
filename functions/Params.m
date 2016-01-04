addpath('..')
addpath('../born/')
addpath('../../pointbem/')


Start = Parameters_of_temp([.5 -60 -.5]);
SIZE  = size(Start);
x0 = [0.5 -60 -0.5];
lb = [0 -500 -100];
ub = [200 0 0];

options = optimoptions('lsqnonlin');
options = optimoptions(options,'Display', 'off');

for i = 1:SIZE(1)
   [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@Parameters_of_temp,x0,lb,ub,options);
end