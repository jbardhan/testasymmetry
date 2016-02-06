% Finding Delta S 
load('Params_NLBC_295K');
load('Params_NLBC_300K');
load('Params_NLBC_305K');

% Temperatures at which the model was parameterized in K
T = [295 300 305];

% Ion index for plotting purposes (Na, K, Rb, Cs, Cl)
Ion = [1 2 3 4 5];

% Parameters at each temperature in T
ParamsAt295K = MakeParamsStruct(x1);
ParamsAt300K = MakeParamsStruct(x2);
ParamsAt305K = MakeParamsStruct(x3);

% Calculated energies using the temperature dependent Parameters 
% calculatedE295K = CalculateEnergiesFromBEM(ParamsAt295K);
% calculatedE300K = CalculateEnergiesFromBEM(ParamsAt300K);
calculatedE305K = CalculateEnergiesFromBEM(ParamsAt305K);
return
% Matrix of all the energies
Energies = [calculatedE295K calculatedE300K calculatedE305K];

% Column vectors of Gibbs free energy for each Ion calculated by
% CalculateEnergiesFromBEM
Na_eng = [Energies(1,1) Energies(1,2) Energies(1,3)];
K_eng = [Energies(2,1) Energies(2,2) Energies(2,3)];
Rb_eng = [Energies(3,1) Energies(3,2) Energies(3,3)];
Cs_eng = [Energies(4,1) Energies(4,2) Energies(4,3)];
Cl_eng = [Energies(5,1) Energies(5,2) Energies(5,3)];

% Slope and intercept values for Gibbs free energy vs. Temperature in K for
% each Ion
Na_values = polyfit(T,Na_eng,1);
K_values = polyfit(T,K_eng,1);
Rb_values = polyfit(T,Rb_eng,1);
Cs_values = polyfit(T,Cs_eng,1);
Cl_values = polyfit(T,Cl_eng,1);

% Data points corresponding to Gibs free energy at a certain temperature
% for each Ion.
F1 = polyval(Na_values,T);
F2 = polyval(K_values,T);
F3 = polyval(Rb_values,T);
F4 = polyval(Cs_values,T);
F5 = polyval(Cl_values,T);

%% Results

% Array of the slopes of the Gibbs free energy vs. Temperature plots for
% each Ion.
Slopes = [Na_values(1) K_values(1) Rb_values(1) Cs_values(1) Cl_values(1)];

% Array of the calculated Delta S values fom our model in J/mol/K
Delta_S_Calc = (-1*4184).*Slopes;

% Array of the intercepts of the Gibbs free energy vs. Temperature plots
% for each Ion
Intercepts = [Na_values(2) K_values(2) Rb_values(2) Cs_values(2) Cl_values(2)];

% Array of the calculated Delta H values from our model in KJ/mol
Delta_H_Calc = Intercepts./.239;

% Array of the reference Delta S values as found by experimentation and
% reported in Fawcett Chp 3. J/mol/K
Ref_Delta_S = [-133 -96 -87 -81 -53];

%% Plots

% Plot Gibbs free energy vs. Temperature for each Ion
plot(T,F1,T,F2,T,F3,T,F4,T,F5)
title('Gibbs Free Energy vs. Temperature (K) for Several Metal Cations and Halide Anions')
xlabel('Temperature (K)')
ylabel('Gibbs Free Energy (KJ/mol)')
legend('Na +1', 'K +1', 'Rb +1', 'Cs +1', 'Cl -1')

figure 

% Compare the calculated Delta S values to those reported by Fawcett Chp 3.
plot(Ion,Delta_S_Calc,'ro',Ion,Ref_Delta_S,'k*')
title('A Compasrison of Calculated \Delta S Values to Experimentally Found \Delta S Values')
xlabel('Ion Index (Na K Rb Cs Cl)')
ylabel('\Delta S (J/mol/K)')
legend('Calculated \Delta S', 'Reported \Delta S')
    