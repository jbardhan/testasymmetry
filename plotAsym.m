x = linspace(-1,1,1000);
ap = asymParams;
figure;
plot(x,ap.alpha*tanh(ap.beta*(x-ap.EfieldOffset))+ (-ap.alpha* ...
						  tanh(ap.beta*(0-ap.EfieldOffset))),'linewidth',2);
hold on;
plot(0,0,'ks','linewidth',2);
clear ap

