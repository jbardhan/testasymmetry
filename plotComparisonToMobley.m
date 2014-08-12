printOn = 1;

% reproduce Figure 4 of Mobley, with our results added
figure; set(gca,'fontsize',16);
plot(Ndistrib(:,1),Ndistrib(:,2),'r^','linewidth',2,'markersize', 10);
hold on;
plot(Pdistrib(:,1),Pdistrib(:,2),'b^','linewidth',2,'markersize', 10);
axis([2 9 -60 0]);
plot(var1_numSides,NdistribBEM,'rs','linewidth',1.6,'markersize', ...
     8);
plot(var1_numSides,PdistribBEM,'bs','linewidth',1.6,'markersize', ...
     8);
xlabel('Number of Atoms in Ring');
ylabel('Electrostatic Solvation Free Energy (kcal/mol)');
legend('Negative (MD)','Positive (MD)',...
       'Negative (Nonlinear B.C.)', 'Positive (Nonlinear B. C.)','location','southwest');
if printOn
  print -depsc2 Mobley-fig4-a.eps
  print -dpng Mobley-fig4-a.png
end

figure; set(gca,'fontsize',16);
plot(Adistrib(:,1),Adistrib(:,2),'g^','linewidth',2,'markersize', ...
     10);
hold on;
plot(var1_numSides,AdistribBEM,'ks','linewidth',1.6,'markersize', ...
     8);
axis([2 9 0 25]);
xlabel('Number of Atoms in Ring');
ylabel('Asymmetry in Electrostatic Free Energy (kcal/mol)');
legend('MD','Nonlinear B.C.');
if printOn
  print -depsc2 Mobley-fig4-b.eps
  print -dpng Mobley-fig4-b.png
end

% Figure 5
figure; set(gca,'fontsize',16); 
plot(Noppose(:,1),Noppose(:,2),'r^','linewidth',2,'markersize', ...
     10);
hold on;
title('(a) Hydration free energies, opposing case');
plot(Poppose(:,1),Poppose(:,2),'b^','linewidth',2,'markersize', ...
     10);
axis([2 9 -70 0]);
plot(var1_numSides,NopposeBEM,'rs','linewidth',1.6,'markersize', ...
     8);
plot(var1_numSides,PopposeBEM,'bs','linewidth',1.6,'markersize', ...
     8);
xlabel('Number of Atoms in Ring');
ylabel('Electrostatic Solvation Free Energy (kcal/mol)');
legend('Negative (MD)','Positive (MD)',...
       'Negative (Nonlinear B.C.)', 'Positive (Nonlinear B. C.)','location','southwest');
if printOn
  print -depsc2 Mobley-fig5-a.eps
  print -dpng Mobley-fig5-a.png
end

figure; set(gca,'fontsize',16);
plot(Aoppose(:,1),Aoppose(:,2),'g^','linewidth',2,'markersize', ...
     10);
hold on;
title('(b) Asymmetry free energies, opposing case');
plot(var1_numSides,AopposeBEM,'ks','linewidth',1.6,'markersize', ...
     8);
axis([2 9 0 25]);
xlabel('Number of Atoms in Ring');
ylabel('Asymmetry in Electrostatic Free Energy (kcal/mol)');
legend('MD','Nonlinear B.C.');
if printOn
  print -depsc2 Mobley-fig5-b.eps
  print -dpng Mobley-fig5-b.png
end


figure; set(gca,'fontsize',16); 
plot(Ndipole(:,1),Ndipole(:,2),'r^','linewidth',2,'markersize', ...
     10);
hold on;
title('(c) Hydration free energies, fixed dipole case');
plot(Pdipole(:,1),Pdipole(:,2),'b^','linewidth',2,'markersize', ...
     10);
axis([2 9 -70 0]);
plot(var1_numSides,NdipoleBEM,'rs','linewidth',1.6,'markersize', ...
     8);
plot(var1_numSides,PdipoleBEM,'bs','linewidth',1.6,'markersize', ...
     8);
xlabel('Number of Atoms in Ring');
ylabel('Electrostatic Solvation Free Energy (kcal/mol)');
legend('Negative (MD)','Positive (MD)',...
       'Negative (Nonlinear B.C.)', 'Positive (Nonlinear B. C.)','location','southwest');
if printOn
  print -depsc2 Mobley-fig5-c.eps
  print -dpng Mobley-fig5-c.png
end

figure; set(gca,'fontsize',16);
plot(Adipole(:,1),Adipole(:,2),'g^','linewidth',2,'markersize', ...
     10);
hold on;
title('(d) Asymmetry free energies, fixed dipole case');
plot(var1_numSides,AdipoleBEM,'ks','linewidth',1.6,'markersize', ...
     8);
axis([2 9 0 25]);
xlabel('Number of Atoms in Ring');
ylabel('Asymmetry in Electrostatic Free Energy (kcal/mol)');
legend('MD','Nonlinear B.C.');
if printOn
  print -depsc2 Mobley-fig5-d.eps
  print -dpng Mobley-fig5-d.png
end
