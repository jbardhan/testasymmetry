printOn = 1;

addpath('../pointbem');
loadConstants
% alpha = 0.3, beta = -50, Efield = -0.5 works a bit better. 
asymParams = struct('alpha', 0.5, 'beta', -60.0,'EfieldOffset',-0.5);
epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;

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
staticPotential =  43.5/4.18;  % from Bardhan12_Jungwirth_Makowski
E_at_m1 = mean(off00.Edata(Iall(1:2))) - (-1.0)*staticPotential;
R_m1 = .5 * conv_factor * (1/epsOut - 1/epsIn) /E_at_m1;
E_at_p1 = mean(off00.Edata(Iall(200:201))) - (1.0) * staticPotential;
R_p1 = .5 * conv_factor * (1/epsOut - 1/epsIn) /E_at_p1;

origin = [0 0 0];
R      = R_p1;  % or substitute with R_m1.  results don't change
                % substantially
		
numCharges = -1;
waterModel = struct('tau',1,'R_OH',0.58,'rho_w',1.4,'R_s',0.52); 
% from Mukhopadhyay

% set up sphere and boundary integral operators
density = 2.0;
numPoints = ceil(4 * pi * density * R^2)
surfdata   = makeSphereSurface(origin, R, numPoints);
surfsurfop = makeSurfaceToSurfaceOperators(surfdata);
clear origin

h = 0.5;
h2 = 1.5;
lineCharges = [0 2 2.5 3 3.5];%lineCharges = 0.0:h:R-h2;
plusOne = 1.0;
minusOne = -1.0;

q_list = linspace(-1,1,40);
for i=1:length(lineCharges)
  for j=1:length(q_list)
    q = q_list(j);
    pqr = struct('xyz',[0 0 lineCharges(i)],'q',q,'R',0);
    bem = makeBemEcfQualMatrices(surfdata, pqr,  epsIn, epsOut);
    
    E = getSelfEnergies(pqr,bem);
    R_eff = getEffectiveRadii(E,epsIn,epsOut);
    R_scaled = getScaledEffectiveRadii(pqr, R_eff, waterModel);
    
    curv_cha(i,j) =  conv_factor * ...
	getAsymmetricGBReactionPotential(pqr.xyz(1,:), R_eff, R_scaled, ...
						       q, ...
						       pqr.xyz(1,:), ...
						       R_eff, R_scaled, ...
						       epsIn, epsOut);
    dG_cha(i,j) = 0.5 * q * curv_cha(i,j);

    [phiReac, sigma] = solveConsistentAsymmetric(surfdata, bem, ...
						 epsIn, epsOut, ...
						 conv_factor, pqr, asymParams);
    dG_bem(i,j) = 0.5 * pqr.q' * phiReac + pqr.q*staticPotential;

%    Efield = -chargesurfop.dphidnCoul * pqr.q;
%    diagPert = diag(surfdata.weights.*(as(1)*(tanh(as(2)*Efield-as(3)))+as(4)));
%    posBem = bem.A + diagPert;
%    sigma = gmres(posBem, rhs, [], 1e-5, 100);
    
%    curv_bem(i,j) = conv_factor * bem.C * sigma;
%    dG_bem(i,j) = 0.5 * pqr.q * curv_bem(i,j) + pqr.q * staticPotential;
  
    rhs = bem.B*pqr.q;    
    sigma_sym = gmres(bem.A,rhs,[],1e-5,100);
    dG_bem_sym(i,j) = 0.5 * conv_factor * pqr.q' * bem.C * sigma_sym + ...
	pqr.q * staticPotential;
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

plot(q_list,dG_bem(1,:),'r-','linewidth',2,'markersize',10)
plot(q_list,dG_bem(2,:),'g--','linewidth',2,'markersize',10)
plot(q_list,dG_bem(3,:),'k-.','linewidth',2,'markersize',10)
plot(q_list,dG_bem(4,:),'b:','linewidth',2,'markersize',10)
plot(q_list,dG_bem(5,:),'m-.','linewidth',2,'markersize',10)

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

