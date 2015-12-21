addpath('..')
addpath('../functions/')
addpath('../born/')
addpath('../../pointbem')
n = 10 ;
q1 = 1;
q2 = -1;
R_list = linspace(1,2.3,20);
k = 1;
epsIn = 1;
conv_factor = 332.112;

lithiumRminOver2 = 1.025; % All radii from Joung08 
% lithiumPlus i = 1

sodiumRminOver2 = 1.369; 
% sodiumPlus i = 2

potassiumRminOver2 = 1.705; 
% potassiumPlus  i = 3

rubidiumRminOver2 = 1.813;
% rubidiumPlus  i = 4

cesiumRminOver2 = 1.976;
% cesiumPlus i = 5

fluorineRminOver2 = 2.303; 
% fluorineMinus i = 6

chlorideRminOver2 = 2.513;
% chlorideMinus  i = 7

bromineRminOver2 = 2.608; 
% bromineMinus i = 8

iodineRminOver2 = 2.860; 
% iodineMinus i = 9

for i = 1:length(R_list)

for T = 0:5:40
    epsOut = epsilon_t(T+273.15);
end
for A = [1.061 1.068 1.075 1.083 1.09 1.098 1.106 1.114 1.123];
    B = -1.*[118.542 120.061 121.606 123.157 124.734 126.394 128.008 129.657 131.343];
    G = -1.*[.857 .872 .888 .904 .92 .937 .953 .969 .986];
Params = struct('alpha', A , 'beta', B,'EfieldOffset',G);
end
    
     
    Radii = R_list(i);
    
% E(i) = bornPicardNoStern(Radii,q1,epsIn,epsOut,k,Params,conv_factor,n);  
% F(i) = bornPicardNoStern(Radii,q2,epsIn,epsOut,k,Params,conv_factor,n)
end