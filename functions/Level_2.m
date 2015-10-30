% Calculate the free energies for each ion.
addpath('..')
addpath('../functions')
addpath('../Asymmetric')
addpath('../Asymmetric/born')
addpath('../pointbem')
loadconstants

rscale = 0.92;

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
Params = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);

i = 1; 
kappa_list = [0.000 0.001 0.02 0.125 0.200 0.500 1.0];

sodiumRminOver2 = 1.41075; % new Roux toppar 1.36375;  % standard charmm
sodiumPlus = -93.4 ;
sodiumMinus = -175.7;

chlorideRminOver2 = 2.27;
chloridePlus = -57.0;
chlorideMinus = -95.3 ;

potassiumRminOver2 =1.76375; % new Roux toppar
potassiumPlus  = -73.4;
potassiumMinus = -128.3895;

rubidiumRminOver2 = 1.90;
rubidiumPlus = -66.78 ;
rubidiumMinus = -114.1 ;

magnesiumRminOver2 = 1.185; % new Roux toppar
magnesiumPlus = -108.6 ;
magnesiumMinus = -218.5 ;

cesiumRminOver2 = 2.1;
cesiumPlus = -60.42 ;
cesiumMinus = -101.9 ; 

calciumRminOver2  = 1.367;
calciumPlus = -88.91 ;
calciumMinus = -163.4 ; 

bariumRminOver2 = 1.89;
bariumPlus = -67.03 ; 
bariumMinus = -115.1 ; 

zincRminOver2 = 1.09;
zincPlus = -99.05 ; 
zincMinus = -191.2 ; 

cadmiumRminOver2 = 1.357;
cadmiumPlus = -89.08; 
cadmiumMinus = -164.3 ; 

for i = 1:8

if i == 1 % NA plus
   
    E(i) = bornPicardNoStern(sodiumRminOver2*rscale,sodiumPlus,epsIn,epsOut,k,Params,conv_factor,k);
   
elseif i == 2 % NA minus
  
    E(i) = bornPicardNoStern(sodiumRminOver2*rscale,sodiumMinus,epsIn,epsOut,k,Params,conv_factor,k);
    
elseif i == 3 % CHLORIDE plus
  
    E(i) = bornPicardNoStern(chlorideRminOver2*rscale,chloridePlus,epsIn,epsOut,k,Params,conv_factor,k);
    
elseif i == 4 % CHLORIDE minus
  
    E(i) = bornPicardNoStern(chlorideRminOver2*rscale,chlorideMinus,epsIn,epsOut,k,Params,conv_factor,k);
   
elseif i == 5 % K plus
   
    E(i) = bornPicardNoStern(potassiumRminOver2*rscale,potassiumPlus,epsIn,epsOut,k,Params,conv_factor,k);
   
elseif i == 6 % K minus
  
    E(i) = bornPicardNoStern(potassiumRminOver2*rscale,potassiumMinus,epsIn,epsOut,k,Params,conv_factor,k);
    
elseif i == 7 % Rb plus
   
    E(i) = bornPicardNoStern(rubidiumRminOver2*rscale,rubidiumPlus,epsIn,epsOut,k,Params,conv_factor,k);
   
elseif i == 8 % Rb minus
  
    E(i) = bornPicardNoStern(rubidiumRminOver2*rscale,rubidiumMinus,epsIn,epsOut,k,Params,conv_factor,k);
    
end
end


D = ObjectiveFunction(i,Params,E)
    
    
    
    
    