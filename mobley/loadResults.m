close all; clear all
addpath('export_fig/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Set your toggles and define the solute and solvent lists %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploton = 1;

%solvents = {'water','octanol','dichloroethane', 'propanone', 'dimethylsulfoxide',...
%             'propanol', 'dimethylformamide', 'ethanol', 'methanol', 'acetonitrile'}; 
solvents = {'water','octanol','dichloroethane'};
%,'octanol','dichloroethane'}; 
%solvents_1 = {'octanol'};
       
%common_solutes = {'ethanol','butanone','n_octane','pyrene','cyclohexane',};
        
testset  = {'butanone','n_octane','ethanol','benzene','cyclohexane','pyrene','n_heptane'};
%testset_ions = {'Li','Na','K','Cl','Br','I'};
%test_ion = {'Rb','Cs'};


%Solute families:
esters = {'11_diacetoxyethane','12_diacetoxyethane','phenyl_formate','trimethoxymethylbenzene',...
    'diethyl_malonate','diethyl_succinate','ethyl_acetate','ethyl_benzoate','ethyl_formate',...
    'isoamyl_acetate','isobutyl_formate','isopropyl_formate','methyl_acetate','methyl_benzoate',...
    'methyl_butanoate','methyl_formate','methyl_hexanoate','methyl_p_methoxybenzoate',...
    'methyl_pentanoate','methyl_propanoate','n_butyl_acetate','n_propyl_acetate',...
    'n_propyl_formate','triacetyl_glycerol'};

alcohols = {'111_trifluoropropan_2_ol','2_methylpropan_2_ol','butan_1_ol','ethanol','heptan_1_ol',...
    'hexan_1_ol','methanol','octan_1_ol','pentan_1_ol','propan_1_ol'};

amines = {'dimethylamine','ethylamine','hydrazine','methylamine','n_butylamine','n_propylamine',...
    '2_naphthylamine','1_naphthylamine','aniline','N_methylaniline','NN_dimethylaniline',...
    'p_toluidine','4_nitroaniline','piperidine','piperazine'};

phenol_derivatives ={'4_methoxyacetophenone','heptan_2_one','33_dimethylbutan_2_one',...
    'acetophenone','butanone','hexan_2_one','octan_2_one','pentan_2_one','pentan_3_one'};

acids = {'3_methylbutanoic_acid','acetic_acid','butanoic_acid','hexanoic_acid','pentanoic_acid','propanoic_acid'};

ketones = {'4_methoxyacetophenone','heptan_2_one','33_dimethylbutan_2_one','acetophenone',...
    'butanone','hexan_2_one','octan_2_one','pentan_2_one','pentan_3_one'};

haloalkanes = {'1_bromo_2_chloroethane','111_trifluoro_222_trimethoxyethane','12_dibromoethane',...
    'dibromomethane','pentachloroethane','tetrachloromethane','tribromomethane'};

nitroalkanes = {'nitromethane','nitroethane','2_nitropropane','2_nitrotoluene',...
    '1_nitrobutane','1_nitropentane','1_nitropropane'};

halobenzenes = {'1234_tetrachlorobenzene','1235_tetrachlorobenzene','1245_tetrachlorobenzene','12_dichlorobenzene','m_bis_trifluoromethyl__benzene','p_dibromobenzene'};
% This routine will plot all of the calulated vs. experimental results if 
% plot is set to 1.  Futhermore, it creates structures that contain all of
% the results for calculated solvation free energies as a function of each 
% solvent and again as a function of solute
for i = 1:length(solvents)
           
    fid = fopen(['mnsol/',lower(solvents{i}),'.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};
    [~,m] = ismember(testset,mol_list);
    [~,est] = ismember(esters,mol_list);
    est(est==0) = [];
    [~,alc] = ismember(alcohols,mol_list);
    alc(alc==0) = [];
    [~,amn] = ismember(amines,mol_list);
    amn(amn==0) = [];
    [~,phn] = ismember(phenol_derivatives,mol_list);
    phn(phn==0) = [];
    [~,acd] = ismember(acids,mol_list);
    acd(acd==0) = [];
    [~,ktn] = ismember(ketones,mol_list);
    ktn(ktn==0) = [];
    [~,nlk] = ismember(nitroalkanes,mol_list);
    nlk(nlk==0) = [];
    [~,hlk] = ismember(haloalkanes,mol_list);
    hlk(hlk==0) = [];
    [~,hbz] = ismember(halobenzenes,mol_list);
    hbz(hbz==0) = [];
    
    %[~,n] = ismember(testset_ions,mol_list);
    
    allCases = 1:length(mol_list);
    %if ~ismember(solvents{i},solvents_1)
    %    [~,p] = ismember(test_ion,mol_list);
    %    allCases(p) = 0;
    %end 
    
    allCases(m) = 0;
    allCases(est) = 0;
    allCases(alc) = 0;
    allCases(acd) = 0;
    allCases(ktn) = 0;
    allCases(nlk) = 0;
    allCases(hlk) = 0;
    allCases(phn) = 0;
    allCases(amn) = 0;
    %allCases(n) = 0;
    
    temp = load(['Run',solvents{i}]);
%     solute_errors(i,:) = temp.errfinal(n);
%    sorted_errors = sort(abs(temp.errfinal));
    results(i) = struct('Solvent', solvents{i},'Num_Solutes',round(length(mol_list),3),'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE),...
                        'Mean_Abs_error',mean(abs(temp.errfinal)));
    if ploton
        figure
        plot(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)),'bo','markers',2)
        %plot(temp.refE,temp.calcE,'bo','markers',12)
        set(gca,'FontSize',15)
        hold on
        plot(temp.refE(m),temp.calcE(m),'rs','markers',12)
        hold on
        plot(temp.refE(est),temp.calcE(est),'r+','markers',12)
        hold on
        plot(temp.refE(alc),temp.calcE(alc),'go','markers',12)
        hold on
        plot(temp.refE(amn),temp.calcE(amn),'b*','markers',12)
        hold on
        plot(temp.refE(acd),temp.calcE(acd),'cx','markers',12)
        hold on
        plot(temp.refE(ktn),temp.calcE(ktn),'ms','markers',12)
        hold on
        plot(temp.refE(nlk),temp.calcE(nlk),'yd','markers',12)
        hold on
        plot(temp.refE(hlk),temp.calcE(hlk),'kv','markers',12)
        hold on
        plot(temp.refE(phn),temp.calcE(phn),'gd','markers',12)
        hold on
        plot(temp.refE(hbz),temp.calcE(hbz),'k*','markers',12)
        
        minax = round(min(min(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)))));
        maxax = round(max(max(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)))));
        axis([minax-2 maxax+2 minax-2 maxax+2]);
        foo = refline(1,0);
        set(foo,'Linewidth',2,'color','k');
        xlabel(['\Delta','G_{expt}^{solv, ',solvents{i},'}',' (kcal.mol^{-1})'])
        ylabel(['\Delta','G_{calc}^{solv, ',solvents{i},'}',' (kcal.mol^{-1})'])
        legend('Predictions ','Training Set','Esters','Alcohols','Amines','Acids','Ketones','Nitroalkanes','Haloalkanes','Phenol derivatives','Halobenzenes','Location','southeast')
      filename = sprintf('Output/Figures/DeltaG-%s.PDF',solvents{i});
      export_fig(filename,'-painters','-transparent');

    end
    
end
clear all;
clear all
addpath('export_fig/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Set your toggles and define the solute and solvent lists %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ploton = 1;

