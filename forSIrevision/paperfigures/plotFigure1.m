load fig2

NL_low = plot(distance,(Enl_part0-Eref_part0)/joulesPerCalorie,'b','linewidth',1.8);
hold on;
set(gca,'FontSize',16);
NL_high = plot(distance2,(Enl_part0_highres-Eref_part0_highres)/ ...
					joulesPerCalorie,'bs','linewidth',1.8);
L_low = plot(distance,(El_part0-Eref_part0)/joulesPerCalorie,'r','linewidth',1.8);
hold on;
set(gca,'FontSize',16);
L_high = plot(distance2,(El_part0_highres-Eref_part0_highres)/ ...
				  joulesPerCalorie,'rs','linewidth',1.8);
XH = xlabel('Charge Position (Angstrom)');
YH = ylabel('Electrostatic Solvation Free Energy (kcal/mol)');
LegendH = legend('Nonlocal, N=1100','Nonlocal, N=2140',...
					  'Local, N=1100','Local, N=2140','location','SouthWest');
print -depsc2 fig2-charge-burial.eps
