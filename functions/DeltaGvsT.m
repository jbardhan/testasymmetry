addpath('..')
addpath('../functions/')
addpath('../born/')
addpath('../../pointbem')
n = 10 ; % Number of picard iterations
q1 = 1; % Positive charge
q2 = -1; % Negative charge
k = 1; % kappa
epsIn = 1;
rscale = 0.92;
conv_factor = 332.112;
T = [0 5 10 15 20 25 30 35 40]+273.15; % Temperature in Kelvin
A = [1.061 1.068 1.075 1.083 1.09 1.098 1.106 1.114 1.123]; % Values for alpha at each temp
B = -1.*[118.542 120.061 121.606 123.157 124.734 126.394 128.008 129.657 131.343]; % Values for beta at each temp
G = -1.*[.857 .872 .888 .904 .92 .937 .953 .969 .986]; % Values for gamma at each temp
Ion = [1 2 3 4 5 6 7 8 9]; % Ion index for plotting purposes
%   [Li Na K Rb Cs F Cl Br I]
for i = 1:9
epsOut(i) = epsilon_t(T(i)); % Epsilon out as a function of temp
Params(i) = struct('alpha', A(i), 'beta', B(i), 'EfieldOffset', G(i)); % Structured array of alpha, beta, gamma as a function of temp
end

%% lithiumPlus i = 1
lithiumRminOver2 = 1.025; % All radii from Joung08 
for l = 1:9
    L_eng(l) = bornPicardNoStern(lithiumRminOver2*rscale,q1,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
L1 = polyfit(T,L_eng,1); % Slope and intercept for Li
F1 = polyval(L1,T); % Fitted line for Li
%% sodiumPlus i = 2
sodiumRminOver2 = 1.369; 
for l = 1:9
    NA_eng(l) = bornPicardNoStern(sodiumRminOver2*rscale,q1,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
NA1 = polyfit(T,NA_eng,1);
F2 = polyval(NA1,T);
%% potassiumPlus  i = 3
potassiumRminOver2 = 1.705; 
for l = 1:9
    K_eng(l) = bornPicardNoStern(potassiumRminOver2*rscale,q1,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
K1 = polyfit(T,K_eng,1);
F3 = polyval(K1,T);
%% rubidiumPlus  i = 4
rubidiumRminOver2 = 1.813;
for l = 1:9
    RB_eng(l) = bornPicardNoStern(rubidiumRminOver2*rscale,q1,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
RB1 = polyfit(T,RB_eng,1);
F4 = polyval(RB1,T);
%% cesiumPlus i = 5
cesiumRminOver2 = 1.976;
for l = 1:9
    CS_eng(l) = bornPicardNoStern(cesiumRminOver2*rscale,q1,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
CS1 = polyfit(T,CS_eng,1);
F5 = polyval(CS1,T);
%% fluorineMinus i = 6
fluorineRminOver2 = 2.303; 
for l = 1:9
    F_eng(l) = bornPicardNoStern(fluorineRminOver2*rscale,q2,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
FL1 = polyfit(T,F_eng,1);
F6 = polyval(FL1,T);
%% chlorideMinus  i = 7
chlorideRminOver2 = 2.513;
for l = 1:9
    CL_eng(l) = bornPicardNoStern(chlorideRminOver2*rscale,q2,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
CL1 = polyfit(T,CL_eng,1);
F7 = polyval(CL1,T);
%% bromineMinus i = 8
bromineRminOver2 = 2.608; 
for l = 1:9
    BR_eng(l) = bornPicardNoStern(bromineRminOver2*rscale,q2,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
BR1 = polyfit(T,BR_eng,1);
F8 = polyval(BR1,T);
%% iodineMinus i = 9
iodineRminOver2 = 2.860; 
for l = 1:9
    I_eng(l) = bornPicardNoStern(iodineRminOver2*rscale,q2,epsIn,epsOut(l),k,Params(l),conv_factor,n);
end
I1 = polyfit(T,I_eng,1);
F9 = polyval(I1,T);

%% Results
Slope = [L1(1) NA1(1) K1(1) RB1(1) CS1(1) FL1(1) CL1(1) BR1(1) I1(1)]; % Array of slopes
Our_Delta_S = (-1/.000239).*Slope; % Delta S in J/mol/K
Intercept = [L1(2) NA1(2) K1(2) RB1(2) CS1(2) FL1(2) CL1(2) BR1(2) I1(2)]; % Array of intercepts
Our_Delta_H = Intercept./.239; % Delta H in KJ/mol
Fawcett_Delta_S = [-199 -143 -101 -92 -78 -138 -89 -79 -66]; % Array of Delta S from Fawcett Chp. 3

%% Plots
plot(T,F1,T,F2,T,F3,T,F4,T,F5,T,F6,T,F7,T,F8,T,F9) % Plots Delta G vs. T for each ion
xlabel('Temp in K')
ylabel('Delta G Kcal/mol')
legend('Lithium +','Sodium +','Potassium +','Rubidium +','Cesium +','Fluorine -','Chloride -','Bromine -','Iodine -')
title('Fitted Delta G vs. Temperature')

figure
plot(Ion,Our_Delta_S,'o',Ion,Fawcett_Delta_S,'*') % A visual aid for comparing our Delta S values to those of Fawcett
axis([1 9 -210 -20])
xlabel('Ion Index')
ylabel('\Delta S in J/mol/K')
title('A Comparison of Fawcett \Delta S and Our \Delta S Values')
legend('Our \Delta S','Fawcett \Delta S','Location','Southeast')