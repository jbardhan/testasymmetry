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

for i = 1:4

if i == 1 % NA plus
   
    E(i) = bornPicardNoStern(sodiumRminOver2*rscale,sodiumPlus,epsIn,epsOut,k,Params,conv_factor,k);
   
elseif i == 2 % NA minus
  
    E(i) = bornPicardNoStern(sodiumRminOver2*rscale,sodiumMinus,epsIn,epsOut,k,Params,conv_factor,k);
    
elseif i == 3 % CHLORIDE plus
  
    E(i) = bornPicardNoStern(sodiumRminOver2*rscale,chloridePlus,epsIn,epsOut,k,Params,conv_factor,k);
    
elseif i == 4 % CHLORIDE minus
  
    E(i) = bornPicardNoStern(sodiumRminOver2*rscale,chlorideMinus,epsIn,epsOut,k,Params,conv_factor,k);
    
end
end


D = ObjectiveFunction(i,Params,E)
    
    
    
    
    