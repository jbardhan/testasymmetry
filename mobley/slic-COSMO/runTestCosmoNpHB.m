function [E,nonpolar,hb,disp,...
          disp_sl_sl,disp_sv_sl,disp_sv_sv,cav,comb] = runTestCosmoNpHB(params, problem, chargeDistribution)
global UsefulConstants ProblemSetNp

%alpha = params.alpha;
%beta  = params.beta;
%gamma = params.EfieldOffset;
%mu    =  -alpha * tanh(-gamma);

%all_atom_types = UsefulConstants.all_atom_types;
epsIn  = UsefulConstants.epsIn;    % 1
epsOut = UsefulConstants.epsOut;   % 80
kappa  = UsefulConstants.kappa;    % 0.0 for non-ionic
conv_factor = UsefulConstants.conv_factor; % 332.112
staticpotential = UsefulConstants.staticpotential;  % 10.7 kcal/mol/e,
                                                    % according to our
                                                    % earlier work
kB = 0.001987;
%% make sure all of our BEM stuff is loaded
index = initializeProblemNp(problem);

% phiReac holds the vector of reaction potentials at the charge
% locations (in kcal/mol/e).  This is the response due to
% electrostatic polarization 

% the additional term is due to the work done against the static
% potential field that arises due to water structure around even
% uncharged solutes.  we model it as a constant field so the extra
% free energy = potential * totalCharge

% now account for the nonpolar solvation and hbonding terms;
% dG_np = dG_disp + dG_comb
% defining variables for Cosmo np and hb calc
solute_data = problem.soluteAtomAreas;
solute_atom_types = problem.soluteAtomTypes;
solute_vdw_v = problem.solute_VdWV;
solute_vdw_a = problem.solute_VdWA;
solute_hbond_data = problem.soluteHbondData;
solvent_vdw_v = problem.solvent_VdWV;
solvent_vdw_a = problem.solvent_VdWA;
solvent_data = problem.solventAtomAreas;
solvent_atom_types = problem.solventAtomTypes;
solvent_hbond_data = problem.solventHbondData;
newhb = problem.newHB;
temp = problem.temperature;
disp_coeffs = params.dispCoeffs;
z_comb = params.zComb;
q_s = params.q_s;
hbond_coeffs = params.hbondCoeffs;
cavity_coeff = params.cavity_coeff;
atom_vols = problem.atom_vols;
area_solute = solute_data(3);
area_solvent = solvent_data(3);
vol_solute = solute_data(4);
vol_solvent = solvent_data(4);
lambda=2;

% Cavity
cav = kB*temp*CalcCavityE(area_solute, vol_solute, atom_vols, cavity_coeff, temp);

% Dispersion
disp_sv_sv = (vol_solute/vol_solvent)*...
             CalcDispersionE(solvent_data, solvent_atom_types, ...
                             solvent_data, solvent_atom_types, ...
                             disp_coeffs, q_s, lambda, temp);

disp_sv_sl = CalcDispersionE(solute_data, solute_atom_types, ...
                             solvent_data, solvent_atom_types, ...
                             disp_coeffs, q_s, lambda, temp);
                
disp_sl_sl = (CalcDispersionE(solute_data, solute_atom_types, ...
                               solute_data, solute_atom_types, ...
                               disp_coeffs, q_s, lambda, temp));
                           
disp = 2*disp_sv_sl - disp_sv_sv + disp_sl_sl;

%Combinatorial
comb = -kB*temp*CalcCombinatorialE(solute_vdw_a, solute_vdw_v, ...
                    solvent_vdw_a, solvent_vdw_v, z_comb);

%% Hydrogen bonding
hbond_coeffs = params.hbondCoeffs;
hb = CalcHBondE(solute_hbond_data, solvent_hbond_data, hbond_coeffs, temp);

nonpolar = comb + disp + cav;

E = nonpolar + hb;
