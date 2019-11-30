close all
clear all
clc

% Path information

Home = getenv('HOME');
repo_path=sprintf('%s/repos',Home);
dropbox_path=sprintf('%s/Dropbox',Home);
addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));                    
addpath(sprintf('%s/testasymmetry/mobley/reference-data',repo_path));                    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%                                Set these values before running the code                      %%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ionflag=0; % ionflag=0 : ions data are not included in the testset 
           % ionflag=1 : ions data are included in the testset

temp_min= 4.85; % lower bound of the temperature interval
temp_max=44.85; % upper bound in the temperature interval
tempdiv=5; % number of temperatures that we calculate solvation free energies 
temp_C=linspace(temp_min,temp_max,tempdiv); % temperatures in Celcius    
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Note: we calculate everything at 298 K which is equal to 24.85C, We
%%%% will calculate \Delta G of ions which is available at 25C at 24.85 C
%%%% to use them in our objective function. But in the postprocessing
%%%% script (whaterthermo.m) we will compare 
                   
% a bunch of useful variables and constants. also defining the global
% variable "ProblemSet" which we'll use to hold the BEM systems.

for j=1:tempdiv
    clear global
    loadConstants
    convertKJtoKcal = 1/joulesPerCalorie;
    global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
    saveMemory = 0;
    writeLogfile = 0;
    logfileName = 'junklogfile';
    epsIn  =  1; % dielectric constant of the solutes
    KelvinOffset = 273.15;
    conv_factor = 332.112;
    staticpotential = 0.0; % this only affects charged molecules;
    kappa = 0.0;  % should be zero, meaning non-ionic solutions!
    epsOut = (-1.410e-6)*temp_C(j)^3+(9.398e-4)*temp_C(j)^2-0.40008*temp_C(j)+87.740;  
    % temperature dependence of the dielectric constant of Water T in C
    % from Mollerup15 (Modeling the permittivity of electrolyte solutions)
    
    UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa',kappa,...
                             'conv_factor',conv_factor,'staticpotential',staticpotential);
    Data = readtable('reference-data/thermo.csv');
    % Experimental data for 502 neutral molecules and 9 monovalent ions
    % from multiple references - See 'reference-data/thermo_expt.xls' for more details
    % Var1 = solutes | 2 = dG_expt | 3 = Surface_Area | 4 = Volume 
    %    5 = dG_np   | 6 = dG_es   | 7 = dG_Mobley

    Data_analogs = readtable('reference-data/analogs.csv');
    % Experimental data for 12 Amino Acid side-chain analog neutral molecules 
    % See 'reference-data/thermo_expt.xls' for more details
    % Var1 = analogs_list | 2 = dG | 3 = dH | 4 = -TdS | 5 = Cp 
    
    Data_aa = readtable('reference-data/aa.csv');
    % Experimental data for 12 Amino Acid side-chains from
    % Hess06 - Hydration Thermodynamic Properties of Amino Acid Analogues
    % Var1 = aa_list | 2 = dG | 3 = dH | 4 = -TdS | 5 = Cp 
    
    Data_full = readtable('reference-data/thermo_full.csv');
    % Full thermodynamic experimental data (dG,dS,dH and Cp) for 111 neutral compounds
    % from multiple references - See 'reference-data/thermo_expt.xls' for more details
    % Var1 = molecules_list | 2 = dG | 3 = dH | 4 = -TdS | 5 = Cp 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ionflag==1

    % Training-set including amino acid side-chain analogs and monovalent ions
    % See 'reference-data/thermo_expt.xls' for more details
    
        testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane',...
                    '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol',...
                    'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane',...
                    'Li','Na','K','Rb','Cs','Cl','Br','I'}; 

    % We can use: testset = [Data_analogs.Var1;Data.Var1(503:511)];
    
    elseif ionflag==0
        testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane',...
                    '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol',...
                    'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane'};

    % We could use: testset = Data_analogs.Var1;            
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    all_solutes = Data.Var1;
    all_dG = Data.Var2;
    all_surfAreas = Data.Var3;
    [m, index] = ismember(testset,all_solutes);
    surfArea_list = all_surfAreas(index);
    
    t_ref_aca=24.85; %reference tempereture for amino acid analogues
    t_ref_ion=25;  %reference tempereture for ions

    dG_ref_298 = Data_analogs.Var2; % kcal/mol
    dH_ref_298 = Data_analogs.Var3; % kcal/mol
    dS_ref_298 =-Data_analogs.Var4/(t_ref_aca + KelvinOffset); % dS = -(-TdS)/T ; kcal/mol/K
    Cp_ref_298 = Data_analogs.Var5/1000;  % kcal/mol/K

    % or dG_ref_298 = [1.99,-9.71,-1.24,2.07,2.32,-1.5,-0.89,-5.1,-5,-5.88,-6.13,1.96]'; 
    %    dH_ref_298=[-2.53,-16.27,-5.22,-5.59,-5.18,-7.97,-8.08,-10.19,-11.97,-14.46,-14.13,-4.87]';
    %    dS_ref_298=(dH_ref_298-dG_ref_298(1:12))/298; 
    %    Cp_ref_298=1e-3*[52.58,26,69.7,92.02,88.58,39.8,85.8,27.34,47.56,99.1,62.02,78.87]';
    % See reference-data/thermo_expt.xls for more details about the experimental data
    
    dG_aca = dG_ref_298 - dS_ref_298*(temp_C(j)-t_ref_aca) + ...
             Cp_ref_298 * ((temp_C(j) - t_ref_aca) - ...
             (temp_C(j) + KelvinOffset) * log(((temp_C(j) + KelvinOffset))/((t_ref_aca + KelvinOffset))));
    % Using the temperature dependence of the standard-state free energy of solvation
    % See Chamberlin08 - Extension of a Temperature-Dependent Aqueous Solvation Model ...
        
    aca_num=length(dG_ref_298);
    ion_num=0;
    
    if ionflag==1
      
        dG_ref_ion_298_15=[-529;-424;-352;-329;-306;-304;-278;-243]./joulesPerCalorie; 
        dS_ref_ion_298_15=[-0.164;-0.133;-0.096;-0.087;-0.081;-0.053;-0.037;-0.014]./joulesPerCalorie;  
        % Data from Fawcett Book Chapter 3 (Note: Data in Fawcett are at 25C which is 298.15K)

        Cp_ref_ion_298_15=1e-3*[-9;-28;-58;-80;-94;-56;-60;-50]./joulesPerCalorie;
        % Data from Abraham86 - The Thermodynamics of Solvation of Ions

        dG_ion = dG_ref_ion_298_15 - dS_ref_ion_298_15*(temp_C(j)-t_ref_ion) ...
               + Cp_ref_ion_298_15 * ((temp_C(j)-t_ref_ion) ...
               - (temp_C(j)+KelvinOffset)*log(((temp_C(j)+KelvinOffset))/((t_ref_ion+KelvinOffset))));
        % Using the temperature dependence of the standard-state free energy of solvation
        % See Chamberlin08 - Extension of a Temperature-Dependent Aqueous Solvation Model ...
        
        dG=[dG_aca;dG_ion];
        dS=[dS_ref_298;dS_ref_ion_298_15];
        Cp=[Cp_ref_298;Cp_ref_ion_298_15];
        
        ion_num=length(dG_ref_ion_298_15);
        
    elseif ionflag==0
        
        dG=dG_aca;
        dS=dS_ref_298;
        Cp=Cp_ref_298;
        
    end
    
    curdir=pwd;
    for i=1:length(testset)
      dir=sprintf('%s/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',dropbox_path,testset{i});
      chdir(dir);
      pqrData = loadPqr('test.pqr');
      pqrAll{i} = pqrData;
      srfFile{i} = sprintf('%s/test_2.srf',dir);
      chargeDist{i} = pqrData.q;%chargeDistribution;
      foo = strcmp(all_solutes,testset{i});
      index = find(foo);
      if length(index) ~= 1
        fprintf('error finding refdata!\n');
        keyboard
      end
      referenceData{i} = dG(i);
      surfArea{i} = surfArea_list(i);
      chdir(curdir);
      addProblemSA(testset{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
    end


    % The following script is specialized to this example.  We'll
    % handle generating others.  Not complicated, but it's not self-explanatory.


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % alpha beta gamma mu phi_stat np_a np_b
    
    
    if ionflag==0
        x0 = [0.476 -121.47 -1.392 -0.59 0.001 0.005 0.079]; 
        % A good initial guess for water from previous calculations
        lb = [0  -200 -100 -20  -0.1  -0.1  -4];
        ub = [+2 +200 +100 +20  +0.1  +0.1  +4];
    
    elseif ionflag==1
        x0 = [1.23 -15.344 0.634 0.498 2.82 0.003 1.644]; 
        % A good initial guess for water from previous calculations
        lb = [0  -200 -100 -20  -20  -0.1  -4];
        ub = [+2 +200 +100 +20  +20  +0.1  +4];
    end

    options = optimoptions('lsqnonlin','MaxIter',8);
    options = optimoptions(options,'Display', 'iter');

    y = @(x)ObjectiveFromBEMSA(x);
    [x,resnorm,residual,exitflag,output] = lsqnonlin(y,x0,lb,ub,options);
    [err,calc,ref,es,np]=ObjectiveFromBEMSA(x);
    [err0,calc0,ref0,es0,np0]=ObjectiveFromBEMSA(x0);

    xvec(j,:)=x;refvec(j,:)=ref;calcvec(j,:)=calc;
    esvec(j,:)=es;npvec(j,:)=np;x0vec(j,:)=x0;
    calc0vec(j,:)=calc0;es0vec(j,:)=es0;np0vec(j,:)=np0;
    tempvec(j,:)=temp_C(j);
    rsme(j,:) = rms(calc-ref);
end
if ionflag==0
    save('OptWater_thermo_wo_ion','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','testset','dS','Cp','tempvec','ionflag','aca_num','ion_num','t_ref_aca','t_ref_ion','rsme');
elseif ionflag==1
    save('OptWater_thermo','xvec','refvec','calcvec','esvec','npvec','x0vec','calc0vec','es0vec','np0vec','testset','dS','Cp','tempvec','ionflag','aca_num','ion_num','t_ref_aca','t_ref_ion','rmse');
end
