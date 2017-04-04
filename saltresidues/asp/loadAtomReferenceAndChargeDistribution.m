chargeDistribution = zeros(length(pqrData.q),2*length(pqrData.q));
for i=1:length(pqrData.q)
  % charging i-th atom only
  chargeDistribution(i,i) = pqrData.q(i);

  % charging all BUT i-th atom
  chargeDistribution(:,length(pqrData.q)+i) = pqrData.q;
  chargeDistribution(i,length(pqrData.q)+i) = 0;
end

qBase =[   0.510000   
      -0.510000   
      -0.470000   
       0.310000   
      -0.110000   
       0.090000   
       0.090000   
       0.090000   
      -0.270000   
       0.090000   
       0.090000   
       0.090000   
       0.510000   
      -0.510000   
      -0.470000   
       0.310000   
       0.070000   
       0.090000   
      -0.210000   
       0.090000   
       0.090000   
       0.620000   
      -0.760000   
      -0.760000   ];
qzero = 0 * qBase;

E_uncorrected = [-6.59304	 
-25.9008	 
-17.4589	 
-2.60901	 
 -1.6926	 
0.238975	 
0.210328	 
0.235826	 
-6.85052    
   0.248414
   0.236573
   0.223705
   -7.71865
   -26.9773
   -18.7485
   -5.03835
   0.392946
   0.275239
   -3.99528
 0.238099 
 0.238929      
 -12.6442     
 -60.388      
-59.8957];      


phi_static = 10.7;
ksi =2.837;
conv_factor = 332.112;
box_length = 40.0;
correction = 0.5 * ksi * conv_factor / box_length;
E_charging_corrected = E_uncorrected - correction * (qBase.^2 - qzero.^2 );

chargeDistribution = chargeDistribution(:,1:length(pqrData.q));
referenceE = E_charging_corrected;

return;

%%%%% BELOW HERE, THE CALCULATIONS HAVE SYSTEMATIC ERROR THAT WE
%HAVE NOT FIXED YET.

E_allbutone_uncorrected = [
    9    -74.1935 ; % from backward.fepout
    8   -88.7565
    7   -88.4678
    6   -89.0825
    5   -75.9911
    4   -116.577
    3   -56.6756
    24  -23.4348
    23  -26.6767
    22  -185.155
    21   -90.3962 ; % from backward.fepout
    20  -90.7311
    2   -61.1742
    19  -61.7508
    18  -90.8791
    17  -88.4998
    16  -106.298
    15  -65.052
    14  -63.2827
    13  -131.192
    12  -87.3746
    11  -87.3007
    10  -87.3733
    1   -137.513
];

[junk,I]=sort(E_allbutone_uncorrected(:,1));
E_allbutone_uncorrected = E_allbutone_uncorrected(I,2);
deltaQSquared = qzero;
for i = 1:length(E_allbutone_uncorrected)
  deltaQSquared(i) = (sum(qBase)-qBase(i))^2 - sum(qzero)^2;
end

E_allbutone_corrected = E_allbutone_uncorrected - correction * deltaQSquared;

referenceE = [E_charging_corrected; E_allbutone_corrected];
