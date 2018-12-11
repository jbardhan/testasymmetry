function [E,electrostatic,nonpolar] = runTestSA3(params, problem, chargeDistribution1, chargeDistribution2)
global UsefulConstants3 ProblemSet3

%alpha = params.alpha;
%beta  = params.beta;
%gamma = params.EfieldOffset;
%mu    =  -alpha * tanh(-gamma);

epsIn1  = UsefulConstants3.epsIn1;    % 1
epsIn2  = UsefulConstants3.epsIn2;    % 1
epsOut = UsefulConstants3.epsOut;   % 80
kappa  = UsefulConstants3.kappa;    % 0.0 for non-ionic
conv_factor = UsefulConstants3.conv_factor; % 332.112
staticpotential = UsefulConstants3.staticpotential;  % 10.7 kcal/mol/e,
                                                    % according to our
                                                    % earlier work

%% make sure all of our BEM stuff is loaded
index = initializeProblem3(problem);

%% assign current test's charge distribution to the PQR
pqrData1 = problem.pqrData1;
pqrData2 = problem.pqrData2;
pqrData1.q = chargeDistribution1;  
pqrData2.q = chargeDistribution2;  
				 
[phiReac1, phiReac2, phiBndy1, dphiDnBndy1,phiBndy2, dphiDnBndy2] = ...
    solvePanelConsistentSternAsym3(ProblemSet3(index).srf1SternData.dielBndy(1),ProblemSet3(index).srf2SternData.dielBndy(1), ...
				  ProblemSet3(index).srf3SternData.sternBndy(1), ...
				  pqrData1, pqrData2, ProblemSet3(index).bemYoonStern, ...
				  epsIn1, epsIn2, epsOut, kappa, ...
				  conv_factor, params, ...
				  ProblemSet3(index).asymBemPcm1,ProblemSet3(index).asymBemPcm2);

% phiReac holds the vector of reaction potentials at the charge
% locations (in kcal/mol/e).  This is the response due to
% electrostatic polarization
dG_asym = 0.5 * (pqrData1.q'* phiReac1 + pqrData2.q'* phiReac2);

% the additional term is due to the work done against the static
% potential field that arises due to water structure around even
% uncharged solutes.  we model it as a constant field so the extra
% free energy = potential * totalCharge
electrostatic = dG_asym + params.phiStatic*sum(sum(pqrData1.q) + sum(pqrData2.q));

% now account for the nonpolar solvation term
nonpolar = params.surfAreaConstant * problem.surfArea + params.NPoffset;
E = electrostatic + nonpolar;
