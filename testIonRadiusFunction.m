printOn = 1;

addpath('../pointbem');
loadConstants
asymParams = struct('alpha', 0.6, 'beta', -20.0,'EfieldOffset',0,'deltaOffset',0);

origin = [0 0 0];
R_list = linspace(1,4,30);
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
sodiumRminOver2 = 1.36375;
sodiumStatic = 11.16;
sodiumPlus = -93.4 - (1.0) * sodiumStatic; 
sodiumMinus = -175.7 - (-1.0) * sodiumStatic; 

chlorideRminOver2 = 2.27;
chlorideStatic = 9.15;
chloridePlus = -57.0 - (1.0) * chlorideStatic;
chlorideMinus = -95.3 - (-1.0) * chlorideStatic;

deltaR = 0.15;
plot([sodiumRminOver2 sodiumRminOver2], [sodiumMinus sodiumPlus], 'ks','markersize',8,'linewidth',2);
%plot([sodiumRminOver2-deltaR sodiumRminOver2+deltaR],[sodiumMinus sodiumMinus],'k-.');
%plot([sodiumRminOver2-deltaR sodiumRminOver2+deltaR],[sodiumPlus sodiumPlus],'k-.');

plot([chlorideRminOver2 chlorideRminOver2], [chlorideMinus chloridePlus], 'ks','markersize',8,'linewidth',2);
%plot([chlorideRminOver2-deltaR chlorideRminOver2+deltaR],[chlorideMinus chlorideMinus],'k-.');
%plot([chlorideRminOver2-deltaR chlorideRminOver2+deltaR],[chloridePlus chloridePlus],'k-.');

if printOn 
  print -depsc2 vary-ion-radius.eps
  print -dpng   vary-ion-radius.png
end
