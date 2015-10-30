% Calculate the free energies for each ion using bornPicardNoStern.m.  This
% script takes each ion and calculates the energy of solvation for a single
% atom of central charge q.  These Energies are then put into an array 'E'
% which is then fed into 'ObjectiveFunction.m' which takes the square of
% the difference between these calculated energies and the experimental MD
% data for each situation.
addpath('..')
addpath('../..')
addpath('../../Asymmetric/born')
addpath('../../pointbem')
loadconstants

rscale = 0.92;

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
Params = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

i = 1; 
k = 1;
n = 10;
q1 = 1;
q2 = -1;


sodiumRminOver2 = 1.41075; % new Roux toppar 1.36375;  % standard charmm
% sodiumPlus i = 1
% sodiumMinus i = 2

chlorideRminOver2 = 2.27;
% chloridePlus i = 3
% chlorideMinus  i = 4

potassiumRminOver2 =1.76375; % new Roux toppar
% potassiumPlus  i = 5
% potassiumMinus  i = 6

rubidiumRminOver2 = 1.90;
% rubidiumPlus  i = 7
% rubidiumMinus  i = 8

magnesiumRminOver2 = 1.185; % new Roux toppar
% magnesiumPlus i = 9
% magnesiumMinus  i = 10

cesiumRminOver2 = 2.1;
% cesiumPlus i = 11
% cesiumMinus i = 12

calciumRminOver2  = 1.367;
% calciumPlus  i = 13
% calciumMinus  i = 14

bariumRminOver2 = 1.89;
% bariumPlus i = 15
% bariumMinus i = 16

zincRminOver2 = 1.09;
% zincPlus  i = 17
% zincMinus  i = 18

cadmiumRminOver2 = 1.357;
% cadmiumPlus i = 19
% cadmiumMinus  i = 20

for i = 1:20

if i == 1 % Na plus
   
    E(i) = bornPicardNoStern(sodiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 2 % Na minus
  
    E(i) = bornPicardNoStern(sodiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 3 % Cl plus
  
    E(i) = bornPicardNoStern(chlorideRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 4 % Cl minus
  
    E(i) = bornPicardNoStern(chlorideRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 5 % K plus
   
    E(i) = bornPicardNoStern(potassiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 6 % K minus
  
    E(i) = bornPicardNoStern(potassiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 7 % Rb plus
   
    E(i) = bornPicardNoStern(rubidiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 8 % Rb minus
  
    E(i) = bornPicardNoStern(rubidiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 9 % Mg plus
   
    E(i) = bornPicardNoStern(magnesiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 10 % Mg minus
  
    E(i) = bornPicardNoStern(magnesiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 11 % Cs plus
   
    E(i) = bornPicardNoStern(cesiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 12 % Cs minus
  
    E(i) = bornPicardNoStern(cesiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 13 % Ca plus
   
    E(i) = bornPicardNoStern(calciumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 14 % Ca minus
  
    E(i) = bornPicardNoStern(calciumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 15 % Ba plus
   
    E(i) = bornPicardNoStern(bariumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 16 % Ba minus
  
    E(i) = bornPicardNoStern(bariumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 17 % Zn plus
   
    E(i) = bornPicardNoStern(zincRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 18 % Zn minus
  
    E(i) = bornPicardNoStern(zincRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
elseif i == 19 % Cd plus
   
    E(i) = bornPicardNoStern(cadmiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
elseif i == 20 % Cd minus
  
    E(i) = bornPicardNoStern(cadmiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
end
end


D = ObjectiveFunction(i,Params,E)
    
    
    
    
    