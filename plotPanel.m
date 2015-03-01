printOn = 1;

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
set(gca,'fontsize',16);
plot([bkProt; bkDeprot] ,[bkProtEnergies; bkDeprotEnergies],'ks','linewidth',4,'markersize',12);
hold on
plot([bkProt; bkDeprot], [dG_pcm(:,1,6); dG_pcm(:,2,6)],'bo','linewidth',2, ...
     'markersize',14);
plot([bkProt; bkDeprot], [dG_pcm(:,1,6); dG_yl(:,2,6)],'rx','linewidth',2, ...
     'markersize',14);

ylabel('Charging free enegy (kcal/mol)')
legend('MD FEP','Panel BEM, PCM', 'Panel BEM, Yoon-Lenhoff','location','southwest');
axis([0 8 -160 0])
set(gca,'Xtick',1:7)
set(gca,'XTickLabel',{'ARG','ASP','CYS','GLU','HIS','LYS','TYR'})
if printOn 
  print -depsc2 panelbem-nosalt-residues.eps
  print -dpng panelbem-nosalt-residues.png
end
