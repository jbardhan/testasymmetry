printOn = 0;

% Figure 9
figure; set(gca,'fontsize',16);
plot(Prods(:,2),'b^','linewidth',2,'markersize',10);
hold on;
title('(a) Hydration Free Energies');
plot(Nrods(:,2),'r^','linewidth',2,'markersize',10);
plot(ProdBEM,'bs','linewidth',1.6,'markersize',8);
plot(NrodBEM,'rs','linewidth',1.6,'markersize',8);
axis([1 8 -130 -10]);
xlabel('Problem Index');
ylabel('Hydration Free Energy (kcal/mol)');
legend('Negative (MD)','Positive (MD)',...
       'Negative (Nonlinear B.C.)', 'Positive (Nonlinear B. C.)','location','southwest');
if printOn
  print -depsc2 Mobley-fig9-a.eps
  print -dpng Mobley-fig9-a.png
end

figure; set(gca,'fontsize',16);
plot(Arods(:,2),'g^','linewidth',2,'markersize',10);
hold on;
title('(b) Asymmetry Free Energies');
plot(ArodBEM,'ks','linewidth',1.6,'markersize',8);
xlabel('Problem Index');
ylabel('Asymmetry in Electrostatic Free Energy (kcal/mol)');
legend('MD','Nonlinear B.C.','location','southeast');
if printOn
  print -depsc2 Mobley-fig9-b.eps
  print -dpng Mobley-fig9-b.png
end
