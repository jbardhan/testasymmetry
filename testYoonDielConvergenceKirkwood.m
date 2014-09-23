addpath('../pointbem');
loadConstants
asymParams = struct('alpha',0.0, 'beta', -60.0,'EfieldOffset',-0.5);

origin = [0 0 0];
R = 5.0;
chargeLocation = [3 0 0];
q = 1;
epsIn  =  2;
epsOut = 80;

densities = 0.4:0.4:3.2;

conv_factor = 332.112;


pqr = struct('xyz',chargeLocation,'q',q,'R',R);

phi_Born = conv_factor * q * doAnalytical(R, epsIn, epsOut, pqr, ...
					     80);
E_born = 0.5 * phi_Born* q';
  
for i=1:length(densities)
  density = densities(i);
  numPoints(i) = ceil(4 * pi * density * R^2);
  
  surfdata   = makeSphereSurface(origin, R, numPoints(i));
  bemEcf = makeBemEcfQualMatrices(surfdata, pqr,  epsIn, epsOut);
  bemYoonDiel = makeBemYoonDielMatrices(surfdata, pqr,  epsIn, epsOut);
  [phiReacEcf, sigma] = solveConsistentAsymmetric(surfdata, bemEcf, ...
						  epsIn, epsOut, ...
						  conv_factor, pqr, asymParams);
  [phiReacYoonDiel, phiBndy,dPhiBndy] = solveConsistentYoonNoSternAsym(surfdata, ...
						  bemYoonDiel, epsIn, ...
						  epsOut, conv_factor, ...
						  pqr, asymParams);
  L_ecf(i) = 0.5 * q'*phiReacEcf;
  L_yoondiel(i) = 0.5 * q'*phiReacYoonDiel;
  fprintf('numPoints = %d, ECF = %f, YL = %f\n',numPoints(i), L_ecf(i), ...
	  L_yoondiel(i));
end

deltaE = [abs(E_born-L_ecf); abs(E_born-L_yoondiel)];
loglog(numPoints, abs(E_born-L_ecf),'rs--','linewidth',2,'markersize',10);
hold on;
axis([min(numPoints) max(numPoints) min(min(deltaE)) max(max(deltaE))]);
loglog(numPoints, abs(E_born-L_yoondiel),'bo-','linewidth',2,'markersize',10);
loglog([numPoints(1) numPoints(end)], 0.33*sum(deltaE(:,1))*[1 sqrt(numPoints(1)/numPoints(end))],'k-.','linewidth',2);
set(gca,'fontsize',16);
legend('ECF Qual', 'Yoon-Lenhoff w/ \kappa=0','Order N^{1/2} convergence')
xlabel('Number of Points');
ylabel('Error (kcal/mol)');
