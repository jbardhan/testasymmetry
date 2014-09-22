printOn = 0;

addpath('../pointbem');
loadConstants
asymParams = struct('alpha', 0.6, 'beta', -20.0,'EfieldOffset',0,'deltaOffset',0);

origin = [0 0 0];
R      = 5.4;
epsIn  =  1;
epsOut = 80;
numCharges = -1;
conv_factor = 332.112;
staticPotential =  43.5/4.18;  % from Bardhan12_Jungwirth_Makowski
waterModel = struct('tau',1,'R_OH',0.58,'rho_w',1.4,'R_s',0.52); 
% from Mukhopadhyay

% set up sphere and boundary integral operators
density = 2.0;
numPoints = ceil(4 * pi * density * R^2)
surfdata   = makeSphereSurface(origin, R, numPoints);

h = 0.5;
h2 = 1.5;
lineCharges = [0 2 2.5 3 3.5];%lineCharges = 0.0:h:R-h2;
plusOne = 1.0;
minusOne = -1.0;

q_list = linspace(-1,1,20);
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
plot(q_list,dG_bem(1,:),'ro-','linewidth',2,'markersize',10)
hold
plot(q_list,dG_bem(2,:),'g^-','linewidth',2,'markersize',10)
plot(q_list,dG_bem(3,:),'k>-','linewidth',2,'markersize',10)
plot(q_list,dG_bem(4,:),'b<-','linewidth',2,'markersize',10)
plot(q_list,dG_bem(5,:),'ms-','linewidth',2,'markersize',10)
xlabel('Charge q');
ylabel('\Delta G^{charging} (kcal/mol)');
legend('0.0','2.0','2.5','3.0','3.5','location','southeast');

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



