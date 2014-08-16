printOn = 1;

% the "natural order" was alphabetical by amino acid.
% for plot clarity we may want to change that with an index set,
% but this'll be a hassle if we also want to plot the Nina et al results.
indexSet = 1:7; 


%  MD results from Nina et al. '97; the protonated histidine(5) results below are from 
% http://thallium.bsd.uchicago.edu/RouxLab/downloads/charmm/pbeq/radii_prot_na.str
rouxProt = [1 3 5 6 7]';
rouxProtEnergies = [-66.0 -15.5 -68.15 -72.9 -19.0]';
rouxDeprot = [2 4 5]';
rouxDeprotEnergies = [-91.9 -87.8 -25.7]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results below are from Roux Born radii * 1.02 using molecular surface
% the above webpage says that the OD* and OE* radii for ASP and GLU
% should be 1.40 not 1.42, as reported in the Nina97 pub
rouxProtPoisson  = [-67.0 -15.4 -68.10 -72.2 -19.1]'; 
rouxDeprotPoisson = [-91.5 -87.5 -24.4]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parseProt = [ -91.8124
  -17.1849
  -20.1325
  -17.7716
  -73.1701
  -81.1003
  -23.8325]; 

parseDeprot = [  -34.6894
  -81.2160
  -71.8941
  -84.7363
  -28.9262
  -17.3594
 -108.8292];

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

figure; set(gca,'fontsize',16);
plot([bkProt; bkDeprot] ,[bkProtEnergies; bkDeprotEnergies],'ks','linewidth',4,'markersize',12);
hold on
plot([bkProt; bkDeprot], [asym(:,1); asym(:,2)],'bo','linewidth',4,'markersize',10);

plot([bkProt; bkDeprot],[sym(:,1); sym(:,2)],'mx','linewidth',2,'markersize',12);

plot([bkProt; bkDeprot],[parseProt; parseDeprot],'g*','linewidth',2,'markersize',12);

plot([rouxProt; rouxDeprot],[rouxProtEnergies; rouxDeprotEnergies],'rs','linewidth',2,'markersize',8);

ylabel('Charging free enegy (kcal/mol)')
legend('MD, this work','NLBC', 'Poisson (Nina radii)',...
       'Poisson (PARSE radii)', 'MD, Nina et al.','location','southwest');

plot(bkProt,bkProtEnergies,'ks','linewidth',4,'markersize',12);
plot(bkDeprot,bkDeprotEnergies,'ks','linewidth',4,'markersize',12);
plot(bkProt, asym(:,1),'bo','linewidth',4,'markersize',10);
plot(bkDeprot, asym(:,2),'bo','linewidth',4,'markersize',10);

axis([0 8 -160 0])
set(gca,'Xtick',1:7)
set(gca,'XTickLabel',{'ARG','ASP','CYS','GLU','HIS','LYS','TYR'})

if printOn 
  print -depsc2 residues.eps
  print -dpng residues.png
end
