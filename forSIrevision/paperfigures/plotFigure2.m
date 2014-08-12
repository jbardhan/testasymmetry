load fig3
El_pair_solv = El_pair - Eref_pair;
El_pair_highdiel_solv = El_pair_highdiel - Eref_pair_highdiel;
Enl_pair_solv = Enl_pair -Eref_pair;

El_self_solv = El_self - Eref_self;
El_self_highdiel_solv = El_self_highdiel - Eref_self_highdiel;
Enl_self_solv = Enl_self - Eref_self;

figure;
NLH = plot(distance, Enl_pair_solv/joulesPerCalorie,'b','linewidth',1.8);
set(gca,'fontsize',16);
hold on;
LH  = plot(distance, El_pair_solv/joulesPerCalorie,'r--','linewidth',1.8);
L_HD_H = plot(distance, El_pair_highdiel_solv/joulesPerCalorie,'k-.','linewidth',1.8);
XH = xlabel('Ion pair position (A)');
YH = ylabel('Electrostatic Solvation Free Energy \Delta \Delta G_{solv} (kcal/mol)');
LegH = legend('Nonlocal','Local, \epsilon_\Omega=4',...
				  'Local, \epsilon_\Omega=15','location','southeast');


El_pairwise_interaction = (El_pair_solv - 2 * El_self_solv)/2;
El_pairwise_interaction_highdiel = (El_pair_highdiel_solv - 2 * El_self_highdiel_solv)/2;
Enl_pairwise_interaction = (Enl_pair_solv - 2 * Enl_self_solv)/2;

figure;
NLH2 = plot(distance, Enl_pairwise_interaction/joulesPerCalorie, ...
				'b','linewidth',1.8);
set(gca,'fontsize',16);
hold on;
LH2 = plot(distance,El_pairwise_interaction/joulesPerCalorie,...
			  'r--','linewidth',1.8);
L_HD_H2 = plot(distance,El_pairwise_interaction_highdiel/joulesPerCalorie,...
					'k-.','linewidth',1.8);
XH2 = xlabel('Ion pair position (A)');
YH2 = ylabel('Screened Pairwise Interaction Energy (kcal/mol)');
LegH2 = legend('Nonlocal','Local, \epsilon_\Omega=4',...
				  'Local, \epsilon_\Omega=15','location','southeast');
