addpath('../../pointbem');
addpath('../../panelbem');
addpath('..');
loadConstants

origin = [0 0 0];
R_list = linspace(1,2.3,20);
% sodium RminOver2 1.41075 from newest Roux paper
% chloride RminOver2 2.27 

sternLayerThickness = 2.0;
chargeLocation = [0 0 0];
q_list = [-1 1];
kappa_list = [0.000 0.001 0.02 0.125 0.200 0.500 1.0];

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);
symParams  = struct('alpha',0.0, 'beta', -60.0,'EfieldOffset',-0.5);

k = 1; 
densityDiel = 4;
densityStern = 1;
staticPotential = 10.0; % approximate static potential. weakly size
                        % dependent, see Ashbaugh00, Bardhan12
			
for i=1:length(R_list)
  R = R_list(i);
  Rstern = R + sternLayerThickness;

  numDielPoints(i) = ceil(4 * pi *  densityDiel * R^2);
  numSternPoints(i) = ceil(4 * pi * densityStern * Rstern^2);
  numTotalPoints(i) = numDielPoints(i) + numSternPoints(i);
  
  dielSurfData = makeSphereSurface(origin, R, numDielPoints(i));
  sternSurfData = makeSphereSurface(origin, Rstern, numSternPoints(i));

  for j=1:length(q_list)
    q = q_list(j);
    pqrData = struct('xyz',chargeLocation,'q',q,'R',R);

    for k=1:length(kappa_list)
      kappa = kappa_list(k);


      bemEcfAsym      = makeBemEcfQualMatrices(dielSurfData, pqrData,  epsIn, epsOut);
      bemStern        = makeBemSternMatrices(dielSurfData, sternSurfData, pqrData, ...
				   epsIn, epsOut, kappa);
      phiReacAsym = solveConsistentSternAsym(dielSurfData, sternSurfData, ...
					 pqrData, bemStern, epsIn, epsOut,...
					 kappa, conv_factor, ...
					 asymParams);
      dG_stern(i,j,k) = 0.5 * pqrData.q' * phiReacAsym;
      withstatic_stern(i,j,k) = dG_stern(i,j,k) + sum(pqrData.q)*staticPotential;
      phiReacSym = solveConsistentSternAsym(dielSurfData, sternSurfData, ...
					 pqrData, bemStern, epsIn, epsOut,...
					 kappa, conv_factor, ...
					 symParams);
      dG_Sym_stern(i,j,k) = 0.5 * pqrData.q' * phiReacSym;
      withstatic_sym(i,j,k) = dG_Sym_stern(i,j,k) + sum(pqrData.q)*staticPotential;
    end
  end
end

% from Bardhan12_Jungwirth_Makowski
rscale = 0.92;
sodiumRminOver2 = 1.41075; % new Roux toppar 1.36375;  % standard charmm
sodiumPlus = -93.4 ;
sodiumMinus = -175.7;

chlorideRminOver2 = 2.27;
chloridePlus = -57.0;
chlorideMinus = -95.3 ;

potassiumRminOver2 =1.76375; % new Roux toppar
potassiumPlus  = -73.4;
potassiumMinus = -128.3895;

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
cadmiumPlus = -89.08; 
cadmiumMinus = -164.3 ; 

radii = [sodiumRminOver2 chlorideRminOver2 potassiumRminOver2 ...
	 rubidiumRminOver2 magnesiumRminOver2 cesiumRminOver2 ...
	 calciumRminOver2 bariumRminOver2 zincRminOver2 ...
	 cadmiumRminOver2];
plus = [sodiumPlus chloridePlus potassiumPlus rubidiumPlus magnesiumPlus ...
	cesiumPlus calciumPlus bariumPlus zincPlus cadmiumPlus];
minus = [sodiumMinus chlorideMinus potassiumMinus rubidiumMinus magnesiumMinus ...
	cesiumMinus calciumMinus bariumMinus zincMinus cadmiumMinus];

stern_nosalt_anion = squeeze(withstatic_stern(:,1,1));
stern_nosalt_cation = squeeze(withstatic_stern(:,2,1));
stern_nosalt_sym = squeeze(withstatic_sym(:,1,1));

figure; 
plot(R_list, stern_nosalt_anion,'r','linewidth',2);
hold on;
set(gca,'fontsize',16);
xlabel('R_{ion} (Angstrom)');
ylabel('Charging free energy (kcal/mol)');

plot(R_list, stern_nosalt_cation,'b','linewidth',2);
plot(R_list, stern_nosalt_sym, 'k--','linewidth',2);

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

axis([1 2.3 -230 -50])
print -depsc2 withstatic-renormalized.eps
print -dpng withstatic-renormalized.png


rootI = 0.304 * kappa_list;

figure;
plot(rootI,squeeze(withstatic_sym(1,1,:))-withstatic_sym(1,1,1), ...
     'k*','markersize',12,'linewidth',2);
set(gca,'fontsize',16);
xlabel('I^{1/2}');
ylabel('\Delta \Delta G^{salt} (kcal/mol)');
hold
plot(rootI,squeeze(withstatic_stern(1,2,:))-withstatic_stern(1,2, ...
						  1),'rs','markersize',8, ...
     'linewidth',2)
plot(rootI,squeeze(withstatic_stern(1,1,:))-withstatic_stern(1,1, ...
						  1),'bo', ...
     'markersize',8,'linewidth',2)
plot(rootI,squeeze(withstatic_sym(1,1,:))-withstatic_sym(1,1,1), ...
     'ks','markersize',10,'linewidth',2);
legend('LPB, q=\pm 1','NLBC, q=-1','NLBC, q=+1','location','southwest');
print -depsc2 salt-dependent-solvation-1A.eps
print -dpng salt-dependent-solvation-1A.png


figure;
plot(rootI,squeeze(withstatic_sym(end,1,:))-withstatic_sym(end,1,1), ...
     'k*','markersize',12,'linewidth',2);
set(gca,'fontsize',16);
xlabel('I^{1/2}');
ylabel('\Delta \Delta G^{salt} (kcal/mol)');
hold
plot(rootI,squeeze(withstatic_stern(end,2,:))-withstatic_stern(end,2, ...
						  1),'rs','markersize',8, ...
     'linewidth',2)
plot(rootI,squeeze(withstatic_stern(end,1,:))-withstatic_stern(end,1, ...
						  1),'bo', ...
     'markersize',8,'linewidth',2)
plot(rootI,squeeze(withstatic_sym(end,1,:))-withstatic_sym(end,1,1), ...
     'ks','markersize',10,'linewidth',2);
%legend('LPB, q=\pm 1','NLBC, q=-1','NLBC, q=+1','location','southwest');

print -depsc2 salt-dependent-solvation-2.3A.eps
print -dpng salt-dependent-solvation-2.3A.png
