printOn = 0;

% results are Roux Born radii * 1.02 using molecular surface
% the protonated histidine(5) results below are from 
% http://thallium.bsd.uchicago.edu/RouxLab/downloads/charmm/pbeq/radii_prot_na.str

% the above webpage says that the OD* and OE* radii for ASP and GLU
% should be 1.40 not 1.42, as reported in the Nina97 pub

% the "natural order" was alphabetical by amino acid.
% for plot clarity we may want to change that with an index set,
% but this'll be a hassle if we also want to plot the Nina et al results.
indexSet = 1:7; 

rouxProt = [1 3 5 6 7];
rouxProtEnergies = [-66.0 -15.5 -68.15 -72.9 -19.0];
rouxProtPoisson  = [-67.0 -15.4 -68.10 -72.2 -19.1]; 
rouxDeprot = [2 4 5];
rouxDeprotEnergies = [-91.9 -87.8 -25.7];
rouxDeprotPoisson = [-91.5 -87.5 -24.4];

ksi = 2.837; % from B Brooks, see Bardhan12_Jungwirth_Makowski
correction = 0.5 * ksi * conv_factor  / 40.0; %box length
bkProt = [1 2 3 4 5 6 7];
bkProtEnergies = [(-54.44-correction) % JR1 (prot=charged)
		  -12.57 % JD1 (prot = neutral)
		  -12.82 % JC1 (prot = neutral)
		  0 % JE1 (prot = neutral) no answer yet
		  (-49.95-correction)  %JH1 (prot=charged)
		  -72.9 % JK1 from Roux!! (prot=neutral)
		  -17.6]; % JY1, from nsteps=2k,needs 8k redone
bkDeprot = [1 2 3 4 5 6 7];
bkDeprotEnergies = [-29.84 % JR2, deprot = neutral
		    (-82.65-correction) % JD2, deprot = charged
		    (-84.46-correction) % JC2, deprot = charged
		    (-85.02-correction) % JE2, deprot = charged
		    -20.75 % JH2, deprot = neutral, needs final answer
		    -15.36  % JK2, deprot = neutral
		    (-103.69-correction)]; % JY2, deprot = charged,
                                           % needs to be redone

figure; set(gca,'fontsize',16);
plot(bkProt,bkProtEnergies,'ks','linewidth',4,'markersize',12);
hold on;
plot(bkDeprot,bkDeprotEnergies,'ks','linewidth',4,'markersize',12);

plot(bkProt, asym(:,1),'bo','linewidth',4,'markersize',10);
plot(bkDeprot, asym(:,2),'bo','linewidth',4,'markersize',10);
% 
plot(bkProt,sym(:,1),'mx','linewidth',2,'markersize',12);
plot(bkDeprot, sym(:,2),'mx','linewidth',2,'markersize',12);

plot(rouxProt,rouxProtEnergies,'rs','linewidth',2,'markersize',8);
plot(rouxDeprot,rouxDeprotEnergies,'rs','linewidth',2,'markersize',8);

ylabel('Charging free enegy (kcal/mol)')

if printOn 
  print -depsc2 residues.eps
  print -dpng residues.png
end
