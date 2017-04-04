function [x,resnorm,residual,exitflag,output,lambda,jacobian] = Optimizer_MSA(x0,lb,ub)
%% This is an auto generated MATLAB file from Optimization Tool.  
% x0 is usually [0.5 -60 -0.5]. lb and ub are usually [0 -500 -100], [100 0 0] respectively.  


%% Start with the default options
options = optimoptions('lsqnonlin');
%% Modify options setting
options = optimoptions(options,'Display', 'off');
[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
lsqnonlin(@ObjectiveFunction_MSA,x0,lb,ub,options);
