printOn=1;
load saltrenormresults.mat
a=[1 0; 0 -1; 0 -1; 0 -1; 1 0; 1 0; 0 -1];
staticpotential = 10.0;

nosalt_sym = [-77.604130 -22.033053
 -13.054572 -81.110342
 -16.422474 -61.510903
 -14.121100 -83.783991
 -64.607177 -18.470316
 -71.879774 -12.394605
 -17.422928 -86.112685];

salt_sym = [-77.750984 -22.041344
 -13.057615 -81.252882
 -16.433034 -61.604308
 -14.123401 -83.927499
 -64.745914 -18.476643
 -72.020604 -12.398500
 -17.426886 -86.251827];

sym = salt_sym - nosalt_sym;

d=squeeze(dG_asym(:,:,3,1))-squeeze(dG_asym(:,:,2,1));
ff=squeeze(dG_asym(:,:,3,1))+d*2;
withstatic_asym=ff+a*staticpotential;


ksi = 2.837; % from B Brooks, see Bardhan12_Jungwirth_Makowski
correction = 0.5 * ksi * conv_factor  / 40.0; %box length
bkProt = [1 2 3 4 5 6 7]';
bkProtEnergies = [(-54.44-correction) % JR1 (prot=charged)
		  -12.57 % JD1 (prot = neutral)
		  -12.82 % JC1 (prot = neutral)
		  -12.48 % JE1 (prot = neutral) NEEDS FINAL ANSWER
		  (-49.96-correction)  %JH1 (prot=charged)
		  (-56.94-correction) % JK1 (prot=charged)
		  -18.29]; % JY1,
bkDeprot = [1 2 3 4 5 6 7]';
bkDeprotEnergies = [-29.84 % JR2, deprot = neutral
		    (-82.65-correction) % JD2, deprot = charged
		    (-84.46-correction) % JC2, deprot = charged
		    (-85.02-correction) % JE2, deprot = charged
		    -25.39 % JH2, deprot = neutral, 
		    -15.36  % JK2, deprot = neutral
		    (-104.68-correction)]; % JY2, deprot = charged,
figure;
plot([1:7 1:7], [bkProtEnergies; bkDeprotEnergies],'ks','markersize',12,'linewidth',4)
hold on
plot([1:7 1:7], [nosalt_sym(:,1); nosalt_sym(:,2)],'m>','markersize',12,'linewidth',4);
plot([1:7 1:7], [withstatic_asym(:,1); withstatic_asym(:,2)],'bo','markersize',12,'linewidth',4)

dG1 = squeeze(dG_asym(:,:,1,1))+a*staticpotential;
dG2 = squeeze(dG_asym(:,:,2,1))+a*staticpotential;
dG3 = squeeze(dG_asym(:,:,3,1))+a*staticpotential;
plot([1:7 1:7], [dG3(:,1); dG3(:,2)],'rd','markersize',10,'linewidth',2);
plot([1:7 1:7], [dG2(:,1); dG2(:,2)],'g*','markersize',10,'linewidth',2);
%plot([1:7 1:7], [dG1(:,1); dG1(:,2)],'m>','markersize',10,'linewidth',2);
plot([1:7 1:7], [bkProtEnergies; bkDeprotEnergies],'ks','markersize',12,'linewidth',4)
plot([1:7 1:7], [withstatic_asym(:,1); withstatic_asym(:,2)],'bo','markersize',12,'linewidth',4)

legend('MD FEP', 'Poisson, Roux radii', 'NLBC+Stern, extrap.','NLBC+Stern, 4 vert/A^2', ...
       'NLBC+Stern, 2 vert/A^2',...%'NLBC+Stern, 1 vert/A^2', ...
       'location','southwest');
       
set(gca,'fontsize',16);
ylabel('Charging free energy (kcal/mol)')
axis([0 8 -180 0])
set(gca,'Xtick',1:7)
set(gca,'XTickLabel',{'ARG','ASP','CYS','GLU','HIS','LYS','TYR'})

if printOn 
  print -depsc2 renormalized_residues.eps
  print -dpng renormalized_residues.png
end


figure;
salt=squeeze(dG_asym(:,:,3,2))-squeeze(dG_asym(:,:,3,1));
plot([1:7],sym(:,1)-sym(:,2),'m>','linewidth',4,'markersize',12);
hold on;
plot([1:7],salt(:,1)-salt(:,2),'bo','linewidth',4,'markersize',12);
set(gca,'fontsize',16);
set(gca,'Xtick',1:7)
set(gca,'XTickLabel',{'ARG','ASP','CYS','GLU','HIS','LYS','TYR'})
ylabel('\Delta \Delta \Delta G^{prot,salt} (kcal/mol)');
legend('Symmetric LPB','NLBC LPB');

if printOn
  print -depsc2 salt_dependent_protonation.eps
  print -dpng   salt_dependent_protonation.png
end
