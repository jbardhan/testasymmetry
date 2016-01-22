function E_MSA = Level_2_MSA(Params,t)
%E_MSA       Returns the energy of an ion of charge q and radius R using
%            funtion bornPicardNoStern.
addpath('..')
addpath('../functions/')
addpath('../born/')
addpath('../../pointbem')


rscale = 0.92;
epsIn  =  1;
for j = 1:length(t)
     epsOut(j) = epsilon_t(t(j)+273.15);
end
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

lithiumRminOver2 = 1.025; %%%%%
% lithiumPlus i = 1

sodiumRminOver2 = 1.369; % new Roux toppar 1.36375;  % standard charmm
% sodiumPlus i = 2

potassiumRminOver2 = 1.705; % new Roux toppar
% potassiumPlus  i = 3

rubidiumRminOver2 = 1.813;
% rubidiumPlus  i = 4

cesiumRminOver2 = 1.976;
% cesiumPlus i = 5

fluorineRminOver2 = 2.303; %%%%%
% fluorineMinus i = 6

chlorideRminOver2 = 2.513;
% chlorideMinus  i = 7

bromineRminOver2 = 2.608; %%%%%
% bromineMinus i = 8

iodineRminOver2 = 2.860; %%%%%
% iodineMinus i = 9


for i = 1:length(t)
% Li plus
    
    E_MSA(i,1) = bornPicardNoStern(lithiumRminOver2*rscale,q1,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);

% Na plus
   
    E_MSA(i,2) = bornPicardNoStern(sodiumRminOver2*rscale,q1,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);
   
% K plus
   
    E_MSA(i,3) = bornPicardNoStern(potassiumRminOver2*rscale,q1,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);
  
% Rb plus
   
    E_MSA(i,4) = bornPicardNoStern(rubidiumRminOver2*rscale,q1,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);
   
% Cs plus
   
    E_MSA(i,5) = bornPicardNoStern(cesiumRminOver2*rscale,q1,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);
 
% F minus
    
    E_MSA(i,6) = bornPicardNoStern(fluorineRminOver2*rscale,q2,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);
 
% Cl minus
  
    E_MSA(i,7) = bornPicardNoStern(chlorideRminOver2*rscale,q2,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);

% Br minus

    E_MSA(i,8) = bornPicardNoStern(bromineRminOver2*rscale,q2,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);
 
% I minus

    E_MSA(i,9) = bornPicardNoStern(iodineRminOver2*rscale,q2,epsIn,...
                 epsOut(i),k,Params,conv_factor,n);
  
end
end    
    
    