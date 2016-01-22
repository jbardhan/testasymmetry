% this script is specialized for the testProblem.m script! do not change
% it without understanding what it's doing!

%  this is the list of atomic charges for protonated ASP, following
%  the protonated residue definition JD1 in the repo (somewhere).
q =[   0.510000   
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
       0.750000   
      -0.550000   
      -0.170000   ];

% it's helpful to have a vector of zero'd out charges.
qzero = 0 * q;

% these are the results directly from MD
E = [-6.59304
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
 -19.8237
 -31.8047
  -3.61322];

%staticpotential. should be deleted.
phi_static = 10.7;  

% this is the correction scale factor to remove the artifacts
% associated with FEP charging simulations when using PBC in NAMD.
% for more information see the NAMD FEP tutorial (Chipot)
ksi =2.837;
conv_factor = 332.112;
box_length = 40.0;
correction = 0.5 * ksi * conv_factor / box_length;

% each of the per-atom energies in E needs its own correction,
% which is the correction factor scaled by the difference in
% squared net charges at the start and the end.  THIS LINE WILL NOT
% WORK FOR "CHARGING ALL BUT ONE"
referenceE = E - correction * (q.^2 - qzero.^2 );

% for the ProblemSet structure, create a charge vector for each
% test.  Here, that means Natoms simulations.  Each one (the ith)
% has the q(i) charge and no others.
chargeDistribution = zeros(length(q));
for i=1:length(q)
  chargeDistribution(i,i) = q(i);
end
