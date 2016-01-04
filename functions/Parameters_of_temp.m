function Parameters = Parameters_of_temp(x)

addpath('..')
addpath('../born/')
addpath('../../pointbem/')


t = 5:10;
for i = 1:length(t)
Parameters(i,:) = ObjectiveFunction_MSA(x,t(i));
end
end
