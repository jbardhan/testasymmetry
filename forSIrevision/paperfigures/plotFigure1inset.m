figure;
R1=plot(E_omega_vec,(Enl_sweepEps(1,:)-Eref_sweepEps(1,:))/joulesPerCalorie,'b','linewidth',1.8);
set(gca,'FontSize',16);
hold on; 
R2=plot(E_omega_vec,(El_sweepEps(1,:)-Eref_sweepEps(1,:))/joulesPerCalorie,'m-.','linewidth',1.8);
R3=plot(E_omega_vec,(Enl_sweepEps(2,:)-Eref_sweepEps(2,:))/joulesPerCalorie,'r--','linewidth',1.8);
R4=plot(E_omega_vec,(El_sweepEps(2,:)-Eref_sweepEps(2,:))/ ...
		  joulesPerCalorie,'k','linewidth',1.8);
XH = xlabel('Solute Dielectric Constant   \epsilon_\Omega');
YH = ylabel('Solvation Free Energy, \Delta G_{solv} (kcal/mol)');
LH = legend('Nonlocal, z=0','Local, z=0', 'Nonlocal, z=12','Local, z=12','location','SouthEast');

Enl_solv_center = Enl_sweepEps(1,:)-Eref_sweepEps(1,:);
Enl_solv_surface = Enl_sweepEps(2,:)-Eref_sweepEps(2,:);
Enl_burial = Enl_solv_center-Enl_solv_surface;
El_solv_center=El_sweepEps(1,:)-Eref_sweepEps(1,:);
El_solv_surface=El_sweepEps(2,:)-Eref_sweepEps(2,:);
El_burial=El_solv_center-El_solv_surface;
figure;
S1 =plot(E_omega_vec,Enl_burial/joulesPerCalorie,'linewidth',1.8);
hold on;
set(gca,'FontSize',16);
S2=plot(E_omega_vec,El_burial/joulesPerCalorie,'r-.','linewidth',1.8); 
XH = xlabel('Solute Dielectric Constant   \epsilon_\Omega');
YH = ylabel('Charge Burial Penalty \Delta \Delta G_{bury} (kcal/mol)');
LH = legend('Nonlocal Model', 'Local Model', 'location', ...
				'NorthEast');
