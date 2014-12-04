printOn = 0;
addpath('..');
addpath('../../pointbem');
loadConstants
asymParamsMobley = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;

% from Bardhan12_Jungwirth_Makowski
rscale = 0.92;
sodiumRminOver2 = 1.41075; % new Roux toppar 1.36375;  % standard charmm
sodiumPlus = -93.4 ;
sodiumMinus = -175.7 ;

chlorideRminOver2 = 2.27;
chloridePlus = -57.0 ;
chlorideMinus = -95.3 ;

potassiumRminOver2 =1.76375; % new Roux toppar
potassiumPlus  = -73.4 ;
potassiumMinus = -128.3895 ;

rubidiumRminOver2 = 1.90;
rubidiumPlus = -66.78 ;
rubidiumMinus = -114.1 ;

magnesiumRminOver2 = 1.185; % new Roux toppar
magnesiumPlus = -108.6 ;
magnesiumMinus = -218.5 ;

cesiumRminOver2 = 2.1;
cesiumPlus = -60.42 ;
cesiumMinus = -101.9 ; 

calciumRminOver2  = 1.367;
calciumPlus = -88.91 ;
calciumMinus = -163.4 ; 

bariumRminOver2 = 1.89;
bariumPlus = -67.03 ; 
bariumMinus = -115.1 ; 

zincRminOver2 = 1.09;
zincPlus = -99.05 ; 
zincMinus = -191.2 ; 

cadmiumRminOver2 = 1.357;
cadmiumPlus = -89.08 ; 
cadmiumMinus = -164.3 ; 

radii = [sodiumRminOver2 chlorideRminOver2 potassiumRminOver2 ...
	 rubidiumRminOver2 magnesiumRminOver2 cesiumRminOver2 ...
	 calciumRminOver2 bariumRminOver2 zincRminOver2 ...
	 cadmiumRminOver2];
plus = [sodiumPlus chloridePlus potassiumPlus rubidiumPlus magnesiumPlus ...
	cesiumPlus calciumPlus bariumPlus zincPlus cadmiumPlus];
minus = [sodiumMinus chlorideMinus potassiumMinus rubidiumMinus magnesiumMinus ...
	cesiumMinus calciumMinus bariumMinus zincMinus cadmiumMinus];

[sortedRadii,I] = sort(radii*rscale);
sortedPlus = plus(I)/conv_factor;
sortedMinus = minus(I)/conv_factor;
Q = 1;
fsym = epsIn/(epsOut-epsIn);

sigmaPlus = 2*epsIn*sortedPlus/(+Q)./sortedRadii;
sigmaMinus = 2*epsIn*sortedMinus/(-Q)./sortedRadii;
EnPlus = (1./(sortedRadii.^2)) * (1+0.5*(1/(1+fsym)));
EnMinus = (-1./(sortedRadii.^2)) * (1+0.5*(1/(1+fsym)));

dPhidnPlus = -(+Q)./sortedRadii.^2;
dPhidnMinus = -(-Q)./sortedRadii.^2;

dPhidn = [dPhidnPlus dPhidnMinus];
sigma = [sigmaPlus sigmaMinus];
h = (1 + fsym) - dPhidn./sigma;

figure;set(gca,'fontsize',16);
plot(dPhidn, dPhidn ./ sigma,'bo','linewidth',2,'markersize',10); % don't change this
                                                  % without changing
                                                  % the last plot
                                                  % command
hold on;

plot(dPhidn, ones(size(dPhidn))*(1+fsym),'r','linewidth',2);

alpha = 0.3;
beta  = -60.0;
gamma = 0.3742;
mu = -0.1073;
plot(dPhidn, 1+fsym-(alpha*tanh(beta*(-dPhidn-0.5*sigma)-gamma)+ ...
		     mu),'ks','linewidth',2);

plot(dPhidn, dPhidn ./ sigma,'bo','linewidth',2); % makes sure this
                                                  % is on top
axis([-1, 1, 0.6, 1.8]);
legend('Explicit-solvent molecular dynamics FEP',...
       'Standard Maxwell boundary condition',...
       'Proposed nonlinear boundary condition');

xlabel('E_n^{Coul}');
ylabel('E_n^{Coul} / \sigma');
if printOn
  print -depsc2 boundary-condition-motivation.eps
  print -dpng boundary-condition-motivation.png
end
