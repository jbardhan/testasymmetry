printOn = 1;
maxDensity = 12;

addpath('../');
addpath('../../pointbem');
loadConstants
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

origin = [0 0 0];
R_list = linspace(1,2.5,5);
q_list = [-1 1];

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
density_list = 1:maxDensity;

maxPicardIterations = 20;
for i=1:length(R_list)
  R= R_list(i);
  for j=1:length(q_list)
    q = q_list(j);
    [E_picard_pcm(i,j), E_picard_yl(i,j)] = bornPicardProblem(R, q, ...
						  epsIn, epsOut, ...
						  asymParams, ...
						  conv_factor, maxPicardIterations);
  end
end

for k=1:length(density_list)
  density = density_list(k);
  for i=1:length(R_list)
    R = R_list(i);
    numPoints(i,k) = ceil(4 * pi * density * R^2);
    surfdata   = makeSphereSurface(origin, R, numPoints(i,k));
    
    for j=1:length(q_list)
      q = q_list(j);

      pqr = struct('xyz',[0 0 0],'q',q,'R',0);
      bemPcm = makeBemEcfQualMatrices(surfdata, pqr,  epsIn, epsOut);
      bemYoonDiel = makeBemYoonDielMatrices(surfdata, pqr,  epsIn, epsOut);
      [phiReacPcm, sigma] = solveConsistentAsymmetric(surfdata, bemPcm, ...
						      epsIn, epsOut, ...
						      conv_factor, pqr, asymParams);
      [phiReacYoonDiel, phiBndy,dPhiBndy] = solveConsistentYoonNoSternAsym(surfdata, bemYoonDiel, ...
						  epsIn, epsOut, ...
						  conv_factor, pqr, asymParams);
      L_pcm(i,j,k) = 0.5 * q'*phiReacPcm;
      L_yoondiel(i,j,k) = 0.5 * q'*phiReacYoonDiel;
      
    end 
  end
  fprintf('Done with density %d of %d\n',k,length(density_list));
end


for i=1:length(R_list)
  for j=1:length(q_list)
    E_extrap_pcm(i,j) = RichardsonExtrapolation(length(density_list)-1,length(density_list),...
						density_list, squeeze(L_pcm(i,j,:)),-0.5);
    E_extrap_yl(i,j) = RichardsonExtrapolation(length(density_list)-1,length(density_list),...
						density_list, squeeze(L_yoondiel(i,j,:)),-0.5);
  end
end

figure;
i = 1;  % which member of R_list
j = 1;  % which member of q_list
axis([10 200 6 70]); % specially chosen for i=1, j=1!
errPCM = abs(E_picard_pcm(i,j)-squeeze(L_pcm(i,j,:)));
loglog(numPoints(i,:),errPCM,'bs', 'linewidth',2);
hold on;
errYL = abs(E_picard_yl(i,j)-squeeze(L_yoondiel(i,j,:)));
loglog(numPoints(i,:),errYL,'ro', 'linewidth',2);
np1 = numPoints(i,1);
np2 = numPoints(i,end);
meanError2 = sqrt(errPCM(end)*errYL(end));
meanError1 = meanError2 * (np1/np2)^(-0.5);
loglog([np1 np2], [meanError1 meanError2],'k','linewidth',1.5);
ax=axis;
axis([0.9*ax(1) 1.1*ax(2) 0.9*ax(3) 1.1*ax(4)]);
set(gca,'fontsize',16);
xlabel('Number of BEM Unknowns');
ylabel('Error relative to Picard solution (kcal/mol)');

if printOn
  print -depsc2 convergence-pointbem-Born.eps
  print -dpng convergence-pointbem-Born.png
end

