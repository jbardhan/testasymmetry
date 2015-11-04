function E = Level_2(Params)
%F           Returns the energy of an ion of charge q and radius R using
%            funtion bornPicardNoStern.
addpath('..')
addpath('../functions/')
addpath('../born/')
addpath('../../pointbem')


rscale = 0.92;

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
Params = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

i = 1; 
k = 1; % Kappa
n = 10; % numPicardIterations
q1 = 1; % plus charge
q2 = -1; % minus charge


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



% Na plus
   
    E(1) = bornPicardNoStern(sodiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Na minus
  
    E(2) = bornPicardNoStern(sodiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
% Cl plus
  
    E(3) = bornPicardNoStern(chlorideRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
    
% Cl minus
  
    E(4) = bornPicardNoStern(chlorideRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
   
% K plus
   
    E(5) = bornPicardNoStern(potassiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% K minus
  
    E(6) = bornPicardNoStern(potassiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
% Rb plus
   
    E(7) = bornPicardNoStern(rubidiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Rb minus
  
    E(8) = bornPicardNoStern(rubidiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
% Mg plus
   
    E(9) = bornPicardNoStern(magnesiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Mg minus
  
    E(10) = bornPicardNoStern(magnesiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
% Cs plus
   
    E(11) = bornPicardNoStern(cesiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Cs minus
  
    E(12) = bornPicardNoStern(cesiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
% Ca plus
   
    E(13) = bornPicardNoStern(calciumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Ca minus
  
    E(14) = bornPicardNoStern(calciumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
% Ba plus
   
    E(15) = bornPicardNoStern(bariumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Ba minus
  
    E(16) = bornPicardNoStern(bariumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
% Zn plus
   
    E(17) = bornPicardNoStern(zincRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Zn minus
  
    E(18) = bornPicardNoStern(zincRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
    
% Cd plus
   
    E(19) = bornPicardNoStern(cadmiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Cd minus
  
    E(20) = bornPicardNoStern(cadmiumRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
  
    
end
    
    
    