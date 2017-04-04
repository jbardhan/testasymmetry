chargeDistribution = zeros(length(pqrData.q),2*length(pqrData.q));
for i=1:length(pqrData.q)
  % charging i-th atom only
  chargeDistribution(i,i) = pqrData.q(i);

  % charging all BUT i-th atom
  chargeDistribution(:,length(pqrData.q)+i) = pqrData.q;
  chargeDistribution(i,length(pqrData.q)+i) = 0;
end

qBase = [   0.510000  
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
       -0.180000  
        0.090000  
        0.090000  
       -0.180000  
        0.090000  
        0.090000  
        0.200000  
        0.090000  
        0.090000  
       -0.700000  
        0.440000  
        0.640000  
       -0.800000  
        0.460000  
        0.460000  
       -0.800000  
        0.460000  
        0.460000  ];
qzero = 0 * qBase;

E_uncorrected=[ -6.84128
 -25.7549
 -18.7281
 -4.76642
 -1.65669
 0.228188
 0.230453
 0.246972
 -6.97647
 0.240295
 0.227364
 0.239682
 -7.13662
 -25.2817
 -18.0041
 -2.58958
 0.378089
 0.236512
 -3.18873
 0.228117
 0.26574
 -3.12232
 0.29188
 0.243501
 -0.110505
 0.230285
 0.24772
 -38.9521
 -13.3932
 -12.8606
 -57.2946
 -11.5657
 -15.2679
 -61.0791
 -15.319
 -15.3959];


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
    1  -64.2118
    10  -52.4662
    11  -54.5752
    12  -53.8913
    13  -66.4498
    14  -65.9664
    15  -75.7858
    16  -56.452
    17  -52.998
    18  -53.1851
    19  -64.2656
    2   -66.2528
    20  -50.3934
    21  -50.5863
    22  -66.1901
    23  -49.4473
    24  -50.4237
    25  -43.5777
    26  -49.0618
    27  -49.6192
    28  -129.118
    29  -33.7936
    3   -75.8733
    30  -29.8699
    31 -151.705
    32  -35.5586
    33  -34.234
    34  -158.213
    35  -33.1607
    36  -33.8967
    4   -49.9055
    5   -56.7101
    6   -54.9939
    7   -53.2833
    8   -54.6948
    9   -61.8066
    ];

[junk,I]=sort(E_allbutone_uncorrected(:,1));
E_allbutone_uncorrected = E_allbutone_uncorrected(I,2);
deltaQSquared = qzero;
for i = 1:length(E_allbutone_uncorrected)
  deltaQSquared(i) = (sum(qBase)-qBase(i))^2 - sum(qzero)^2;
end

E_allbutone_corrected = E_allbutone_uncorrected - correction * deltaQSquared;

referenceE = [E_charging_corrected; E_allbutone_corrected];
