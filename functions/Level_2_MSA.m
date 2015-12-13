function E_MSA = Level_2_MSA(Params)
%F           Returns the energy of an ion of charge q and radius R using
%            funtion bornPicardNoStern.
addpath('..')
addpath('../functions/')
addpath('../born/')
addpath('../../pointbem')


rscale = 0.92;
t = 0.1:5:99;
epsIn  =  1;
epsOut = epsilon_t(30+273.15);
conv_factor = 332.112;


alpha = Params.alpha;
beta = Params.beta;
gamma = Params.EfieldOffset;
mu =  -alpha * tanh(-gamma);


i = 1; 
k = 1; % Kappa
n = 10; % numPicardIterations
q1 = 1; % plus charge
q2 = -1; % minus charge

lithiumRminOver2 = 1; %%%%%
% lithiumPlus i = 1

sodiumRminOver2 = 1.41075; % new Roux toppar 1.36375;  % standard charmm
% sodiumPlus i = 2

potassiumRminOver2 = 1.76375; % new Roux toppar
% potassiumPlus  i = 3

rubidiumRminOver2 = 1.90;
% rubidiumPlus  i = 4

cesiumRminOver2 = 2.1;
% cesiumPlus i = 5

fluorineRminOver2 = 1; %%%%%
% fluorineMinus i = 6

chlorideRminOver2 = 2.27;
% chlorideMinus  i = 7

bromineRminOver2 = 1; %%%%%
% bromineMinus i = 8

iodineRminOver2 = 1; %%%%%
% iodineMinus i = 9



% Li plus
    
    E_MSA(1) = bornPicardNoStern(lithiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);

% Na plus
   
    E_MSA(2) = bornPicardNoStern(sodiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% K plus
   
    E_MSA(3) = bornPicardNoStern(potassiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
  
% Rb plus
   
    E_MSA(4) = bornPicardNoStern(rubidiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
   
% Cs plus
   
    E_MSA(5) = bornPicardNoStern(cesiumRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
 
% F minus
    
    E_MSA(6) = bornPicardNoStern(fluorineRminOver2*rscale,q1,epsIn,epsOut,k,Params,conv_factor,n);
 
% Cl minus
  
    E_MSA(7) = bornPicardNoStern(chlorideRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);

% Br minus

    E_MSA(8) = bornPicardNoStern(bromineRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
 
% I minus

    E_MSA(9) = bornPicardNoStern(iodineRminOver2*rscale,q2,epsIn,epsOut,k,Params,conv_factor,n);
  
end
    
    
    