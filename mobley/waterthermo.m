clear all
%close all

Home = getenv('HOME');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%                                Set these values before running the code                      %%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

repo_path=sprintf('%s/Research',Home);

ionflag=0;          % if ionflag=0, ions data are not included in the testset 
                    % if ionflag=1, ions data are included in the testset
paramboundflag=1;   % if paramboundflag=0 there is no bound for parameters in the optimization process 
                    % if paramboundflag=1 there is abound
                    
temp_min=5;     % lower bound of the temperature interval 
temp_max=45;    % upper bound in the temperature interval
tempdiv=5;      % number of divisions in the temperature interval                     
                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));


TEMP=linspace(5,45,tempdiv);        % create the temperature vector                   

if ionflag==1
    testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane','Li','Na','K','Rb','Cs','F','Cl','Br','I'};
elseif ionflag==0
    testset  = {'methane', 'ethanamide', 'methanethiol', 'n_butane', '2_methylpropane', 'methyl_ethyl_sulfide', 'toluene', 'methanol', 'ethanol', '3_methyl_1h_indole', 'p_cresol', 'propane'};
end


dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures

loadConstants
convertKJtoKcal = 1/joulesPerCalorie;

if ionflag==1
    if paramboundflag==0
        ParamWatInfo = load('OptWater_w_ion_wo_bound');
    elseif paramboundflag==1
        ParamWatInfo = load('OptWater_w_ion');
    end
elseif ionflag==0
    if paramboundflag==0
        ParamWatInfo = load('OptWater_wo_ion_wo_bound');
    elseif paramboundflag==1
        ParamWatInfo = load('OptWater_wo_ion');
    end
end
    
        
x = ParamWatInfo.xvec;
calcvec = ParamWatInfo.calcvec;

for i=1:length(testset)
    f=fit(TEMP',calcvec(:,i),'poly1'); % the linear function that fits to the calculated dG  at differet temperatures
    dGfunc(i).name=testset(i); 
    dGfunc(i).func=f;
end

for i=1:length(testset)
    dgvec(i)=dGfunc(i).func(24.85);  % the evaluation of the dG function at 298K = 24.85C in kcal/mol/K
    p=[dGfunc(i).func.p1,dGfunc(i).func.p2];
%     p=[thermofunc(i).func.p1,thermofunc(i).func.p2,thermofunc(i).func.p3,thermofunc(i).func.p4,thermofunc(i).func.p5];
    pder=polyder(p);  % derivative of the linear function dG
    dsvec(i)=-polyval(pder,24.85)*1000;    % Evaluationg entropy, dS at 298K = 24.85C in cal/mol/K
%     figure(i);
%     plot(dGfunc(i).func,TEMP',calcvec(:,i),'o')   
end


if ionflag==1
    dG_list_ref_at_298=[8.1,-40.5,-5.2,9,9.5,-6.2,-3.2,-21.2,-20.4,-24.6,-25.6,8.3]'./joulesPerCalorie; %Hess in kcal/mol
    H_list_ref_at_298=[-8.3,-67.0,-23.9,-17.1,-17.1,-34.6,-25.3,-43.0,-45,-58.8,-57.4,-13.7]'./joulesPerCalorie;  %Hess in kcal/mol
    dS_list_ref_at_298=(H_list_ref_at_298-dG_list_ref_at_298(1:12))/298*1000;  %  in cal/mol/K

    dG_list_ref_ion_at_298=[-529;-424;-352;-329;-306;-429;-304;-278;-243]./joulesPerCalorie;         %Fawcett (Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference
    dS_list_ref_ion_at_298=[-0.164;-0.133;-0.096;-0.087;-0.081;-0.115;-0.053;-0.037;-0.014]./joulesPerCalorie*1000;   %Fawcett(Data in Fawcett are at 25C which is 298.15K. I ignored that 0.15K difference

    dG_list_ref_at_298=[dG_list_ref_at_298;dG_list_ref_ion_at_298]';
    dS_list_ref_at_298=[dS_list_ref_at_298;dS_list_ref_ion_at_298]';


elseif ionflag==0
    dG_list_ref_at_298=[8.1,-40.5,-5.2,9,9.5,-6.2,-3.2,-21.2,-20.4,-24.6,-25.6,8.3]./joulesPerCalorie;%Hess in kcal/mol
    H_list_ref_at_298=[-8.3,-67.0,-23.9,-17.1,-17.1,-34.6,-25.3,-43.0,-45,-58.8,-57.4,-13.7]./joulesPerCalorie;%Hess  in kcal/mol
    dS_list_ref_at_298=(H_list_ref_at_298-dG_list_ref_at_298)/298*1000;  % in cal/mol/K

end

dG_list_ref_at_298
dgvec

dS_list_ref_at_298
dsvec

dgerr=abs(dG_list_ref_at_298-dgvec)
dserr=abs(dS_list_ref_at_298-dsvec)



figure()
scatter(dG_list_ref_at_298,dS_list_ref_at_298,100*ones(1,length(testset)),'x','linewidth',2);
hold on
scatter(dgvec,dsvec,100*ones(1,length(testset)),'o','linewidth',2);
hold on
%line([-100 100],[0 0])
ax=gca;
ax.LineWidth=1;
ax.Box='on';
set(gca,'fontsize',30);
xlabel('\DeltaG (kcal/mol)');
ylabel('\DeltaS (cal/mol K)');
legend('Experiment', 'Calculated');
if ionflag==1
    if paramboundflag==1
        title('\DeltaS VS. \DeltaG with ions and with bounded parameters in optimization')
    elseif paramboundflag==0
        title('\DeltaS VS. \DeltaG with ions but without bounded parameters  in optimization')
    end
elseif ionflag==0
    if paramboundflag==1
        title('\DeltaS VS. \DeltaG without ions but with bounded parameters  in optimization')
    elseif paramboundflag==0
        title('\DeltaS VS. \DeltaG without ions and without bounded parameters in optimization')
    end
end

figure()
scatter(1:length(testset),dgerr,100*ones(1,length(testset)),'x','linewidth',2);
ax=gca;
ax.LineWidth=1;
ax.Box='on';
set(gca,'fontsize',30);
xlabel('solute');
ylabel('Error in \DeltaG (kcal/mol)');
if ionflag==1
    if paramboundflag==1
        title('Absolute error in \DeltaG  for each solute with ions and with bounded parameters in optimization')
    elseif paramboundflag==0
        title('Absolute error in \DeltaG  for each solute with ions and with out bounded parameters in optimization')
    end
elseif ionflag==0
    if paramboundflag==1
        title('Absolute error in \DeltaG  for each solute with out ions and with bounded parameters in optimization')
    elseif paramboundflag==0
        title('Absolute error in \DeltaG  for each solute with out ions and with out bounded parameters in optimization')
    end
end

figure()
scatter(1:length(testset),dserr,100*ones(1,length(testset)),'s','linewidth',2);
ax=gca;
ax.LineWidth=1;
ax.Box='on';
set(gca,'fontsize',30);
xlabel('solute');
ylabel('Error in \DeltaS (cal/mol K)');
if ionflag==1
    if paramboundflag==1
        title('Absolute error in \DeltaS for each solute with ions and with bounded parameters in optimization')
    elseif paramboundflag==0
        title('Absolute error in \DeltaS for each solute with ions and with out bounded parameters in optimization')
    end
elseif ionflag==0
    if paramboundflag==1
        title('Absolute error in \DeltaS for each solute with out ions and with bounded parameters in optimization')
    elseif paramboundflag==0
        title('Absolute error in \DeltaS for each solute with out ions and with out bounded parameters in optimization')
    end
end
        
    
    
    