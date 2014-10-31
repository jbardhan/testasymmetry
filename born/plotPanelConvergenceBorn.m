printOn = 0;

addpath('../');
addpath('../../pointbem');
addpath('../../panelbem');
loadConstants
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);


origin = [0 0 0];
q_list = [-1 1];
% R_list = linspace(1,2.5,5) becomes
R_list = {'1'};%, '1.375', '1.75', '2.125', '2.5'}; 
density_list = 1:6;
epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;


maxPicardIterations = 20;
for i=1:length(R_list)
  R= str2double(R_list{i});
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
    Rstr = R_list{i};
    R = str2double(Rstr);
    srfFile = sprintf('born_%sA_%d.srf',Rstr, density);
    srfData = loadSrfIntoPanels(srfFile);
    
    for j=1:length(q_list)
      q = q_list(j);

      pqrData = struct('xyz',[0 0 0],'q',q,'R',0);
      bemPcm = makePanelBemEcfQualMatrices(srfData, pqrData,  epsIn, epsOut);
      asymBem = makePanelAsymCollocMatrices(srfData, bemPcm, pqrData);
      [phiReacPcm, sigma] = solvePanelConsistentAsymmetric(srfData, ...
						  bemPcm, epsIn, ...
						  epsOut, conv_factor, ...
						  pqrData, asymParams,asymBem);
      L_pcm(i,j,k) = 0.5 * pqrData.q'*phiReacPcm;
      
    end 
  end
  fprintf('Done with density %d of %d\n',k,length(density_list));
end


for i=1:length(R_list)
  for j=1:length(q_list)
    E_extrap_pcm(i,j) = RichardsonExtrapolation(length(density_list)-1,length(density_list),...
						density_list, squeeze(L_pcm(i,j,:)),-0.5);
  end
end
return

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

