function D = Differenceof(y,f,x)
%           Computeds the difference between data points and a fit line for
%           a polynomial of degree 2.
for i = 1:length(x)
D(i) = [(y(i) - f(i)).^2];
end