function [E,electrostatic,nonpolar,hb,disp,...
          disp_sl_sl,disp_sv_sl,disp_sv_sv,comb] = runTestCosmo(params, problem, chargeDistribution)
global UsefulConstants ProblemSet

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
index = initializeProblem(problem);
%% assign current test's charge distribution to the PQR
pqrData = problem.pqrData;
pqrData.q = chargeDistribution;  
				 
[phiReac, phiBndy, dphiDnBndy] = ...
    solvePanelConsistentSternAsym(ProblemSet(index).srfSternData.dielBndy(1), ...
				  ProblemSet(index).srfSternData.sternBndy(1), ...
				  pqrData, ProblemSet(index).bemYoonStern, ...
				  epsIn, epsOut, kappa, ...
				  conv_factor, params, ...
				  ProblemSet(index).asymBemPcm);

% phiReac holds the vector of reaction potentials at the charge
% locations (in kcal/mol/e).  This is the response due to
% electrostatic polarization 
dG_asym = 0.5 * pqrData.q' * phiReac;

% the additional term is due to the work done against the static
% potential field that arises due to water structure around even
% uncharged solutes.  we model it as a constant field so the extra
% free energy = potential * totalCharge
electrostatic = dG_asym + params.phiStatic*sum(pqrData.q);

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
%all_atom_types = problem.allAtomTypess;
temp = problem.temperature;
disp_coeffs = params.dispCoeffs;
hbond_coeffs = params.hbondCoeffs;
z_comb = params.zComb;
q_s = params.q_s;
area_solute = solute_data(3);
area_solvent = solvent_data(3);
vol_solute = solute_data(4);
vol_solvent = solvent_data(4);

% hbond term
if newhb==1
    hb = CalcHbondE(solute_hbond_data,solvent_hbond_data,hbond_coeffs)-...
    area_solute/area_solvent*CalcHbondE(solvent_hbond_data,solvent_hbond_data,hbond_coeffs);
else
    hb = CalcHbondE(solute_hbond_data,solvent_hbond_data,hbond_coeffs);
end

% dispersion term

disp_sv_sv = (vol_solute/vol_solvent)*CalcDispE(solvent_data,solvent_atom_types,...
              solvent_data,solvent_atom_types,disp_coeffs,q_s);

disp_sv_sl = CalcDispE(solute_data,solute_atom_types,...
                    solvent_data,solvent_atom_types,...
                    disp_coeffs,q_s);
                
disp_sl_sl = -(CalcDispE(solute_data,solute_atom_types,...
                    solute_data,solute_atom_types,...
                    disp_coeffs,q_s));
disp_sl = disp_sl_sl;
disp_sv = (disp_sv_sv - 2*disp_sv_sl);              
disp = disp_sl - disp_sv;

% combinatorial term
comb = -kB*temp*CalcCombE(solute_vdw_a,solute_vdw_v,...
                    solvent_vdw_a,solvent_vdw_v,z_comb);
nonpolar = comb + disp;

E = electrostatic + nonpolar + hb;
