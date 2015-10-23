printOn = 0;

addpath('../pointbem');
loadConstants
% alpha = 0.3, beta = -50, Efield = -0.5 works a bit better. 
asymParams = struct('alpha', 0.5, 'beta', -60.0,'EfieldOffset',-0.5);
epsIn  =  1;
epsOut = 80;
kappa = 0.00;
conv_factor = 332.112;
sternLayerThickness = 2.0;
staticPotential =  43.5/4.18;  % from Bardhan12_Jungwirth_Makowski

curdir = pwd;
%%%%%%%%%%%% Needs to be replaced with actual data
clear origin
cd ~/research/static/newsphere/offsetsFine_sph4
cd ../sph4
origin
cd ../offsetsFine_sph4

xOrigin = localfit(Qdata, Edata);
off00 = struct('Qdata',Qdata,'Edata',Edata);

loadOff05
x05 = localfit(Qdata, Edata);
newLoadOff10
x10 = localfit(Qdata, Edata);
newLoadOff15
x15 = localfit(Qdata, Edata);
newLoadOff20
x20 = localfit(Qdata, Edata);
newLoadOff25
x25 = localfit(Qdata, Edata);
fixLoad30 %%%%%%% 
x30 = localfit(Qdata, Edata);
newLoadOff35
x35 = localfit(Qdata, Edata);
%%%%%%%%%%%%%%
cd(curdir);
[junk,Iall]=sort(off00.Qdata);
E_at_m1 = mean(off00.Edata(Iall(1:2))) - (-1.0)*staticPotential;
R_m1 = .5 * conv_factor * (1/epsOut - 1/epsIn) /E_at_m1;
E_at_p1 = mean(off00.Edata(Iall(200:201))) - (1.0) * staticPotential;
R_p1 = .5 * conv_factor * (1/epsOut - 1/epsIn) /E_at_p1;


origin = [0 0 0];
R      = R_m1;  % or substitute with R_m1.  results don't change
                % substantially
Rstern = R + sternLayerThickness;

% set up sphere and boundary integral operators
densityDiel = 2.0;
densityStern = 1;
numDielPoints = ceil(4 * pi * densityDiel * R^2);
numSternPoints = ceil(4 * pi * densityStern *Rstern^2);

dielSurfData   = makeSphereSurface(origin, R, numDielPoints);
sternSurfData  = makeSphereSurface(origin, Rstern, numSternPoints);
clear origin

h = 0.5;
h2 = 1.5;
lineCharges = [0 2 2.5 3 3.5];%lineCharges = 0.0:h:R-h2;
plusOne = 1.0;
minusOne = -1.0;

q_list = linspace(-1,1,40);
for i=1:length(lineCharges)
    pqr = struct('xyz',[0 0 lineCharges(i)],'q',0,'R',0);
    bemEcfAsym = makeBemEcfQualMatrices(dielSurfData, pqr,  epsIn, epsOut);
    bemStern   = makeBemSternMatrices(dielSurfData, sternSurfData, ...
				      pqr, epsIn, epsOut, kappa);
  for j=1:length(q_list)
    pqr.q = q_list(j);
    phiReacAsym = solveConsistentSternAsym(dielSurfData, sternSurfData, ...
					   pqr, bemStern, epsIn, ...
					   epsOut, kappa, conv_factor, ...
					   asymParams);
    dG_bem(i,j) = 0.5 * pqr.q' * phiReacAsym;
    withstatic(i,j) = dG_bem(i,j) + sum(pqr.q)*staticPotential;
  end
end

figure; set(gca,'fontsize',16)
[junk,Iall] = sort(off00.Qdata);
plot(off00.Qdata(Iall),off00.Edata(Iall),'ro','markersize',8,'linewidth',1.5);
hold on;

[junk,Iall] = sort(off20.Qdata);
plot(off20.Qdata(Iall),off20.Edata(Iall),'g^','markersize',8,'linewidth',1.5);

[junk,Iall] = sort(off25.Qdata);
plot(off25.Qdata(Iall),off25.Edata(Iall),'k>','markersize',8,'linewidth',1.5);

[junk,Iall] = sort(off30.Qdata);
plot(off30.Qdata(Iall),off30.Edata(Iall),'b<','markersize',8,'linewidth',1.5);

[junk,Iall] = sort(off35.Qdata);
plot(off35.Qdata(Iall),off35.Edata(Iall),'ms','markersize',8,'linewidth',1.5);

plot(q_list,withstatic(1,:),'r-','linewidth',2,'markersize',10)
plot(q_list,withstatic(2,:),'g--','linewidth',2,'markersize',10)
plot(q_list,withstatic(3,:),'k-.','linewidth',2,'markersize',10)
plot(q_list,withstatic(4,:),'b:','linewidth',2,'markersize',10)
plot(q_list,withstatic(5,:),'m-.','linewidth',2,'markersize',10)

xlabel('Charge q');
ylabel('\Delta G^{charging} (kcal/mol)');
legend('FEP 0.0','FEP 2.0','FEP 2.5','FEP 3.0','FEP 3.5',...
       'NLBC 0.0', 'NLBC 2.0', 'NLBC 2.5', 'NLBC 3.0', 'NLBC 3.5',...
       'location','south');

return


figure; set(gca,'fontsize',16)
plot(q_list,dG_bem_sym(1,:),'ro-','linewidth',2,'markersize',10)
hold
plot(q_list,dG_bem_sym(2,:),'g^-','linewidth',2,'markersize',10)
plot(q_list,dG_bem_sym(3,:),'k>-','linewidth',2,'markersize',10)
plot(q_list,dG_bem_sym(4,:),'b<-','linewidth',2,'markersize',10)
plot(q_list,dG_bem_sym(5,:),'ms-','linewidth',2,'markersize',10)
title('Affine Symmetric Model')
xlabel('Charge q');
ylabel('\Delta G^{charging} (kcal/mol)');
legend('0.0','2.0','2.5','3.0','3.5','location','southeast');

return

figure; set(gca,'fontsize',16);
plot(lineCharges,dG_cat_cha,'bo-','linewidth',2)
hold
plot(lineCharges,dG_cat_bem,'bo--','linewidth',1.5,'markersize',10)
plot(lineCharges,dG_an_cha,'rs-','linewidth',2)
plot(lineCharges,dG_an_bem,'rs--','linewidth',1.5,'markersize',10)
legend('\Delta G_{GB} q=+1', '\Delta G_{BEM} q=+1', ...
       '\Delta G_{GB} q=-1', ...
       '\Delta G_{BEM} q=-1',...
       'location','southwest');

