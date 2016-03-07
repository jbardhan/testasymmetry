chargeDistribution = zeros(length(pqrData.q),2*length(pqrData.q));
for i=1:length(pqrData.q)
  % charging i-th atom only
  chargeDistribution(i,i) = pqrData.q(i);

  % charging all BUT i-th atom
  chargeDistribution(:,length(pqrData.q)+i) = pqrData.q;
  chargeDistribution(i,length(pqrData.q)+i) = 0;
end


qBase = [     0.5100
   -0.5100
   -0.4700
    0.3100
   -0.1100
    0.0900
    0.0900
    0.0900
   -0.2700
    0.0900
    0.0900
    0.0900
    0.5100
   -0.5100
   -0.4700
    0.3100
    0.0700
    0.0900
   -0.1800
    0.0900
    0.0900
         0
   -0.1150
    0.1150
   -0.1150
    0.1150
   -0.1150
    0.1150
   -0.1150
    0.1150
   -0.1150
    0.1150
  ];
qzero = 0 * qBase;

E_uncorrected=[  -6.67354
 -25.7141
 -18.2164
 -5.02364
 -1.67479
 0.268575
 0.231123
 0.240413
 -6.88561
 0.236076
 0.239537
 0.234325
 -7.03194
 -25.4544
 -18.0491
 -3.24099
 0.389057
 0.272644
 -3.13542
 0.238167
 0.25344
 -8.3766e-13
 -1.67777
 0.105415
 -1.74802
 0.0776271
 -1.71333
 0.0867925
 -1.71373
 0.101804
 -1.76353
 0.101094];

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
    1    -44.9508
    10   -19.7939
    11   -21.8541
    12   -21.0374
    13   -49.6887
    14   -12.4668
    15   -21.1646
    16   -26.524
    17   -19.7665
    18   -20.3333
    19   -19.1584
    2    -13.0063
    20   -19.0892
    21   -18.7206
    22   -18.5089
    23   -18.271
    24   -19.2292
    25   -18.183
    26   -19.9462
    27   -18.2636
    28   -20.1759
    29   -18.5493
    3    -25.5236
    30   -19.2954
    31   -17.9449
    32   -19.9246
    4    -22.8079
    5    -17.8739
    6    -20.5895
    7    -19.1916
    8    -20.5607
    9    -15.7878
    ];
[junk,I]=sort(E_allbutone_uncorrected(:,1));
E_allbutone_uncorrected = E_allbutone_uncorrected(I,2);
deltaQSquared = qzero;
for i = 1:length(E_allbutone_uncorrected)
  deltaQSquared(i) = (sum(qBase)-qBase(i))^2 - sum(qzero)^2;
end

E_allbutone_corrected = E_allbutone_uncorrected - correction * deltaQSquared;

referenceE = [E_charging_corrected; E_allbutone_corrected];
