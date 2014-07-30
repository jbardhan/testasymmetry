addpath('../pointbem');
loadConstants
as = [0.6 -16 0 0];

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

    rhs = bem.B*pqr.q;
    
    Efield = -chargesurfop.dphidnCoul * pqr.q;
    diagPert = diag(surfdata.weights.*(as(1)*(tanh(as(2)*Efield-as(3)))+as(4)));
    newBem = bem.A + diagPert;
    sigma = gmres(newBem, rhs, [], 1e-5, min(100,size(bem.A,1)));
    
    L(i,j) = kcal_to_kJ* conv_factor * (0.5 * q * bem.C * sigma);
    
  end    
end
  
figure; set(gca,'fontsize',16);
plot(R_list, L(:,1),'r','linewidth',2); hold on;
plot(R_list, L(:,end),'b--','linewidth',2);
legend('q = -1', 'q = +1');
xlabel('R_{ion}');
ylabel('Charging free energy (kJ/mol)');
