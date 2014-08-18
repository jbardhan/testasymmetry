printOn = 1;

addpath('../pointbem');
loadConstants
% alpha = 0.3, beta = -60, EfieldOffset = 0 works well for the ions
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

origin = [0 0 0];
R_list = linspace(1,2.5,30);
q_list = [-1 1]; %linspace(-1,1,20);

epsIn  =  1;
epsOut = 80;
numCharges = -1;
conv_factor = 332.112;
kcal_to_kJ = 4.18;
waterModel = struct('tau',1,'R_OH',0.58,'rho_w',1.4,'R_s',0.52); 
staticPotential = 10.7;
density = 4.0;


for i=1:length(R_list)
  R = R_list(i);
  numPoints = ceil(4 * pi * density * R^2)
  surfdata   = makeSphereSurface(origin, R, numPoints);
  surfsurfop = makeSurfaceToSurfaceOperators(surfdata);

  for j=1:length(q_list)
    q = q_list(j);
    pqr = struct('xyz',[0 0 0],'q',q,'R',0);
    chargesurfop = makeSurfaceToChargeOperators(surfdata, pqr);
    bem = makeBemMatrices(surfdata, pqr, surfsurfop, ...
			  chargesurfop,  epsIn, epsOut);

    [phiReac, sigma] = solveConsistentAsymmetric(surfdata, surfsurfop, ...
						 chargesurfop, bem, ...
						 epsIn, epsOut, ...
						 conv_factor, pqr, asymParams);
    L(i,j) = 0.5 * q'*phiReac;
  end 
  E_born(i) = 0.5 * conv_factor * (1/epsOut - 1/epsIn)/R;
end
  
figure; set(gca,'fontsize',16);
plot(R_list, L(:,1),'r-','linewidth',2,'markersize',10); hold on;
plot(R_list, L(:,end),'b-.','linewidth',2,'markersize',10);
plot(R_list, E_born,'k--','linewidth',1.5, 'markersize',10);
legend('q = -1', 'q = +1','location','southeast');
xlabel('R_{ion} (Angstrom)');
ylabel('Charging free energy (kcal/mol)');

% from Bardhan12_Jungwirth_Makowski
rscale = 0.92;
sodiumRminOver2 = 1.41075; % new Roux toppar 1.36375;  % standard charmm
sodiumStatic = 0; % 11.16;
sodiumPlus = -93.4 - (1.0) * sodiumStatic; 
sodiumMinus = -175.7 - (-1.0) * sodiumStatic; 

chlorideRminOver2 = 2.27;
chlorideStatic = 0; % 9.15;
chloridePlus = -57.0 - (1.0) * chlorideStatic;
chlorideMinus = -95.3 - (-1.0) * chlorideStatic;

potassiumRminOver2 =1.76375; % new Roux toppar
potassiumStatic = 0; % 10.028;
potassiumPlus  = -73.4 - (1.0) * potassiumStatic;
potassiumMinus = -128.3895 - (-1.0) * potassiumStatic;

rubidiumRminOver2 = 1.90;
rubidiumStatic = 0; % 9.3987;
rubidiumPlus = -66.78 - (1.0) * rubidiumStatic;
rubidiumMinus = -114.1 - (-1.0) * rubidiumStatic;

magnesiumRminOver2 = 1.185; % new Roux toppar
magnesiumStatic = 0; %13.4133;
magnesiumPlus = -108.6 - (1.0) * magnesiumStatic;
magnesiumMinus = -218.5 - (-1.0) * magnesiumStatic;;

cesiumRminOver2 = 2.1;
cesiumStatic  = 0; %9.373;
cesiumPlus = -60.42 - (1.0) * cesiumStatic;
cesiumMinus = -101.9 - (-1.0) * cesiumStatic;

calciumRminOver2  = 1.367;
calciumStatic = 0; %10.8085;
calciumPlus = -88.91 - (1.0) * calciumStatic;
calciumMinus = -163.4 - (-1.0) * calciumStatic;

bariumRminOver2 = 1.89;
bariumStatic = 0; %9.7107;
bariumPlus = -67.03 - (1.0) * bariumStatic;
bariumMinus = -115.1 - (-1.0) * bariumStatic;

zincRminOver2 = 1.09;
zincStatic = 0; %12.3103;
zincPlus = -99.05 - (1.0) * zincStatic;
zincMinus = -191.2 - (-1.0) * zincStatic;

cadmiumRminOver2 = 1.357;
cadmiumStatic = 0; %10.7392;
cadmiumPlus = -89.08 - (1.0) * cadmiumStatic;
cadmiumMinus = -164.3 - (-1.0) * cadmiumStatic;

deltaR = 0.15;
plot(rscale*[sodiumRminOver2 sodiumRminOver2], [sodiumMinus sodiumPlus], 'ks','markersize',10,'linewidth',2);

plot(rscale*[potassiumRminOver2 potassiumRminOver2], [potassiumMinus ...
		    potassiumPlus], 'bs','markersize',10,'linewidth',2);

plot(rscale*[rubidiumRminOver2 rubidiumRminOver2], [rubidiumMinus rubidiumPlus], 'rs','markersize',10,'linewidth',2);

plot(rscale*[cesiumRminOver2 cesiumRminOver2], [cesiumMinus cesiumPlus], 'ms','markersize',10,'linewidth',2);

plot(rscale*[chlorideRminOver2 chlorideRminOver2], [chlorideMinus chloridePlus], 'k*','markersize',10,'linewidth',2);

plot(rscale*[magnesiumRminOver2 magnesiumRminOver2], [magnesiumMinus magnesiumPlus], 'bo','markersize',10,'linewidth',2);

plot(rscale*[calciumRminOver2 calciumRminOver2], [calciumMinus calciumPlus], 'ro','markersize',10,'linewidth',2);

plot(rscale*[bariumRminOver2 bariumRminOver2], [bariumMinus bariumPlus], 'ko','markersize',10,'linewidth',2);

plot(rscale*[zincRminOver2 zincRminOver2], [zincMinus zincPlus], 'md','markersize',10,'linewidth',2);

plot(rscale*[cadmiumRminOver2 cadmiumRminOver2], [cadmiumMinus cadmiumPlus], 'kd','markersize',10,'linewidth',2);

legend('NLBC, q = -1', 'NLBC, q = +1', 'Poisson', 'Na \pm 1', 'K \pm 1',...
       'Rb \pm 1', 'Cs \pm 1', 'Cl \pm 1', ...
       'Mg \pm 1', 'Ca \pm 1', 'Ba \pm 1', 'Zn \pm 1', 'Cd \pm 1', 'Location','EastOutside');

%plot([sodiumRminOver2-deltaR sodiumRminOver2+deltaR],[sodiumMinus sodiumMinus],'k-.');
%plot([sodiumRminOver2-deltaR sodiumRminOver2+deltaR],[sodiumPlus sodiumPlus],'k-.');
%plot([chlorideRminOver2-deltaR chlorideRminOver2+deltaR],[chlorideMinus chlorideMinus],'k-.');
%plot([chlorideRminOver2-deltaR chlorideRminOver2+deltaR],[chloridePlus chloridePlus],'k-.');

if printOn 
  print -depsc2 vary-ion-radius.eps
  print -dpng   vary-ion-radius.png
end
