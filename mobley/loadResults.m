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
solvents = {'water','octanol','hexadecane','diethylether','hexane','carbontet','benzene','toluene','dichloroethane','methanol','butanol','methanol','butanol','xylene','cyclohexane','dibutylether','chloroform'};
%,'octanol','dichloroethane'}; 
%solvents_1 = {'octanol'};
       
%common_solutes = {'ethanol','butanone','n_octane','pyrene','cyclohexane',};
        
testset  = {'butanone','n_octane','ethanol','toluene','nitromethane','14_dioxane','phenol'};
flag_testset = 0;
%testset_ions = {'Li','Na','K','Cl','Br','I'};
%test_ion = {'Rb','Cs'};


%Solute families:

alkanes = {'1_methylcyclohexene','2_methylbutane','2_methylhexane','2_methylpentane','2_methylpropane','2_methylpropene',...
           '22_dimethylbutane','22_dimethylpentane','22_dimethylpropane','224_trimethylpentane','225_trimethylhexane','23_dimethylbutane',...
            '23_dimethylpentane','234_trimethylpentane','24_dimethylpentane','3_methylheptane','3_methylhexane','3_methylpentane',...
            '33_dimethylpentane','cis_12_dimethylcyclohexane','cyclohexane','cyclopentane','cyclopropane','ethane','methane','methylcyclohexane',...
            'methylcyclopentane','n_butane','n_decane','n_heptane','n_hexane','n_nonane',...
            'n_octane','n_pentane','n_pentylcyclopentane','n_propylcyclopentane','propane','trans_14_dimethylcyclohexane'};
flag_lka =1;

%%%%%%%%%%%%%%%%%%%

alkenes = {'2_methyl_but_2_ene','2_methylbut_2_ene','2_methylbuta_13_diene','2_methylpent_1_ene',...
           '23_dimethylbuta_13_diene','3_methyl_but_1_ene','3_methylbut_1_ene','but_1_ene','buta_13_diene','cyclohepta_135_triene','cyclohexene',...
           'cyclopentene','E_hept_2_ene','ethene','hept_1_ene','hex_1_ene','hexa_15_diene',...
           'non_1_ene','oct_1_ene','pent_1_ene','penta_14_diene','propene','styrene','Z_pent_2_ene'};
flag_lke = 1;

%%%%%%%%%%%%%%%%%%%

aromatic_hydrocarbons ={'1_ethylnaphthalene','1_methylnaphthalene','123_trimethylbenzene','124_trimethylbenzene','13_dimethylnaphthalene',...
                        '135_trimethylbenzene','14_dimethylnaphthalene','2_ethyltoluene','23_dimethylnaphthalene','26_dimethylnaphthalene','4_ethyltoluene',...
                        '4_isopropyltoluene','acenaphthene','alpha_methylstyrene','anthracene','benzene','biphenyl','ethylbenzene','fluorene',...
                        'indane','isobutylbenzene','isopropylbenzene','m_xylene','n_butylbenzene','n_hexylbenzene','n_pentylbenzene','n_propylbenzene',...
                        'naphthalene','o_xylene','p_xylene','phenanthrene','pyrene','sec_butylbenzene','tert_butylbenzene','toluene'};
flag_ahc = 1;

%%%%%%%%%%%%%%%%%%%

esters = {'11_diacetoxyethane','12_diacetoxyethane','diethyl_malonate','diethyl_succinate','ethyl_acetate','ethyl_benzoate',...
          'ethyl_butanoate','ethyl_formate','ethyl_hexanoate','ethyl_pentanoate','ethyl_phenyl_ether','ethyl_propanoate','isoamyl_acetate',...
          'isoamyl_formate','isobutyl_acetate','isobutyl_formate','isobutyl_isobutanoate','isopropyl_acetate','isopropyl_formate','methyl_acetate',...
          'methyl_benzoate','methyl_butanoate','methyl_cyanoacetate','methyl_cyclohexanecarboxylate','methyl_cyclopropanecarboxylate','methyl_formate','methyl_hexanoate',...
          'methyl_methanesulfonate','methyl_octanoate','methyl_pentanoate','methyl_p_methoxybenzoate','methyl_propanoate','methyl_trimethylacetate','n_butyl_acetate','n_hexyl_acetate',...
          'n_pentyl_acetate','n_pentyl_propanoate','n_propyl_acetate','n_propyl_butyrate','n_propyl_formate','n_propyl_propanoate','phenyl_formate'};
flag_est = 0;

%%%%%%%%%%%%%%%%%%%

alcohols = {'2_methylbutan_1_ol','2_methylbutan_2_ol','2_methylpentan_2_ol','2_methylpentan_3_ol','2_methylpropan_1_ol',...
            '2_methylpropan_2_ol','3_methylbutan_1_ol','4_methylpentan_2_ol','butan_1_ol','butan_2_ol','cycloheptanol',...
            'cyclohexanol','cyclopentanol','decan_1_ol','ethanol','heptan_1_ol','hexan_1_ol',...
            'hexan_3_ol	methanol','nonan_1_ol','octan_1_ol','pentan_1_ol','pentan_2_ol',...
            'pentan_3_ol','prop_2_en_1_ol','propan_1_ol','propan_2_ol','3_phenylpropanol','2_phenylethanol','benzyl_alcohol'};
flag_alc = 0;

%%%%%%%%%%%%%%%%%%%

aldehydes = {'acetaldehyde','butyraldehyde','E_but_2_enal','E_hex_2_enal','E_oct_2_enal','formaldehyde','heptanal',...
             'hexanal','isobutyraldehyde','nonanal','octanal','pentanal','propionaldehyde','4_methylbenzaldehyde','benzaldehyde'};
flag_ald = 0;

%%%%%%%%%%%%%%%%%%%

ethers = {'11_diethoxyethane','111_trimethoxyethane','12_diethoxyethane','12_dimethoxyethane','14_dioxane','2_methoxy_111_trimethoxyethane',...
          'di_n_butyl_ether','di_n_propyl_ether','diethyl_ether','diisopropyl_ether','dimethoxymethane','dimethyl_ether',...
          'methyl_ethyl_ether','methyl_isopropyl_ether','methyl_propyl_ether','methyl_t_butyl_ether','methyl_tert_butyl_ether','tetrahydrofuran','tetrahydropyran','trimethoxy_methane'};
flag_eth = 0;

%%%%%%%%%%%%%%%%%%%

amides = {'benzamide','ethanamide','n_butylacetamide','N_methylacetamide','NN_dimethyl_p_methylbenzamide','NN_dimethylbenzamide','NN_dimethylformamide'};
flag_amd = 0;

%%%%%%%%%%%%%%%%%%%

amines_1 = {'1_naphthylamine','2_naphthylamine','aniline','cyclohexylamine','ethylamine','hydrazine',...
            'methylamine','n_butylamine','n_heptylamine','n_hexylamine','n_octylamine','n_pentylamine','n_propylamine',...
            'o_toluidine','p_toluidine'};
flag_am1 = 0;

%%%%%%%%%%%%%%%%%%%

amines_2 = {'1_methyl_pyrrole','14_dimethyl_piperazine','4_methyl_1h_imidazole','azetidine','di_n_butylamine',...
            'di_n_propylamine','diethylamine','diisopropylamine','dimethylamine','imidazole','morpholine','N_methylaniline',...
            'N_methylpiperazine','piperazine','piperidine','pyrrole','pyrrolidine'};
flag_am2 = 0;

%%%%%%%%%%%%%%%%%%%

amines_3 = {'1_methyl_imidazole','2_ethylpyrazine','2_ethylpyridine','2_isobutylpyrazine','2_methylpyrazine','2_methylpyridine','23_dimethylpyridine',...
            '24_dimethylpyridine','25_dimethylpyridine','26_dimethylpyridine','3_cyanopyridine','3_ethylpyridine','3_methylpyridine',...
            '34_dimethylpyridine','35_dimethylpyridine','4_cyanopyridine','4_ethylpyridine','4_methylpyridine',...
            'N_methylmorpholine','N_methylpiperidine','NN_dimethylaniline','pyridin','quinoline','triethylamine','trimethylamine','N_acetylpyrrolidine'};
flag_am3 = 0;

%%%%%%%%%%%%%%%%%%%

ether_amines = {'2_methoxyaniline','2_methoxyethanamine','3_methoxyaniline','333_trimethoxypropionitrile','4_methoxyaniline'};
flag_eam = 0;

%%%%%%%%%%%%%%%%%%%

ketone_amines = {'3_acetylpyridine','3_formylpyridine','4_acetylpyridine','4_formylpyridine'};
flag_kam = 0;

%%%%%%%%%%%%%%%%%%%

halo_amines = {'2_chloroaniline','2_chloropyridine','3_chloroaniline','3_chloropyridine','4_chloroaniline','N_methyl_N__222_trifluoroethyl__aniline'};
flag_ham = 0;

%%%%%%%%%%%%%%%%%%%

nitro_compounds = {'1_nitrobutane','1_nitropentane','1_nitropropane','2_nitropropane','nitroethane','nitromethane','2_nitroaniline',...
                '3_nitroaniline','4_nitroaniline','2_nitrotoluene','3_nitrotoluene','3_nitrophenol','4_nitrophenol','2_nitrophenol','NN_dimethyl_p_nitrobenzamide'};
flag_nit = 1;

%%%%%%%%%%%%%%%%%%%

nitriles = {'benzonitrile','acetonitrile','butanenitrile','pentanenitrile','propanenitrile'};
flag_ntr = 0;

%%%%%%%%%%%%%%%%%%%

phenols ={'1_naphthol','2_chlorophenol','2_fluorophenol','2_iodophenol','2_naphthol','23_dimethylphenol',...
          '24_dimethylphenol','25_dimethylphenol','26_dimethylphenol','3_ethylphenol','34_dimethylphenol',...
          '35_dimethylphenol','4_ethylphenol','4_n_propylphenol','4_tert_butylphenol','m_cresol','o_cresol',...
          '2_methoxyphenol','3_methoxyphenol','p_cresol','phenol'};
flag_phn = 0;

%%%%%%%%%%%%%%%%%%%

acids = {'3_methylbutanoic_acid','acetic_acid','butanoic_acid','hexanoic_acid','pentanoic_acid','propanoic_acid'};
flag_acd = 1;

%%%%%%%%%%%%%%%%%%%

ketones = {'24_dimethylpentan_3_one','3_methylbutan_2_one','33_dimethylbutan_2_one','4_methylpentan_2_one','butanone','cyclohexanone',...
           'cyclopentanone','decan_2_one','heptan_2_one','heptan_4_one','hexan_2_one','methyl_cyclohexyl_ketone','methyl_cyclopropyl_ketone',...
           'nonan_2_one','nonan_5_one','octan_2_one','pentan_2_one','pentan_3_one','propanone','undecan_2_one'};
flag_ktn = 0;

%%%%%%%%%%%%%%%%%%% 

haloalkanes = {'1_bromo_2_chloroethane','1_bromo_2_methylpropane','1_bromobutane','1_bromoheptane','1_bromohexane','1_bromooctane',...
               '1_bromopentane','1_bromopropane','1_chloro_222_trifluoroethane','1_chlorobutane','1_chloroheptane','1_chlorohexane',...
               '1_chloropentane','1_chloropropane','1_iodobutane','1_iodoheptane','1_iodohexane','1_iodopentane','1_iodopropane','11_dichloroethane',...
               '11_difluoroethane','111_trichloroethane','1112_tetrachloroethane','112_trichloro_122_trifluoroethane','112_trichloroethane',...
               '1122_tetrachloroethane','12_dibromoethane','12_dichloroethane','12_dichloropropane','13_dichloropropane','14_dichlorobutane',...
               '2_bromo_2_methylpropane','2_bromopropane','2_chloro_111_trimethoxyethane','2_chloro_2_methylpropane','2_chlorobutane',...
               '2_chloropropane','2_iodopropane','bromoethane','bromoethane','bromomethane','bromotrifluoromethane','chlorodifluoromethane',...
               'chloroethane','chlorofluoromethane','chloromethane','dibromomethane','dichloromethane','diiodomethane','fluoromethane',...
               'halothane','iodoethane','iodoethane','iodomethane','pentachloroethane','teflurane','tetrachloromethane','tetrafluoromethane',...
               'tribromomethane','trichloroethene','trichloromethane'};
flag_hlk = 0;

%%%%%%%%%%%%%%%%%%%

haloalkenes = {'11_dichloroethene','3_chloroprop_1_ene','chloroethylene','E_12_dichloroethene','hexafluoropropene','tetrachloroethene','Z_12_dichloroethene'};
flag_hle = 0;

%%%%%%%%%%%%%%%%%%%

halobenzenes = {'12_dichlorobenzene','123_trichlorobenzene','1234_tetrachlorobenzene','1235_tetrachlorobenzene','124_trichlorobenzene',...
                '1245_tetrachlorobenzene','13_dichlorobenzene','135_trichlorobenzene','14_dichlorobenzene','2_chlorotoluene','4_bromotoluene',...
                'benzotrifluoride','benzyl_chloride','bromobenzene','chlorobenzene','fluorobenzene','iodobenzene','p_dibromobenzene','m_bis_trifluoromethyl__benzene',};
flag_hbz = 0;

%%%%%%%%%%%%%%%%%%%

halo_phenols = {'3_chlorophenol','4_bromophenol','4_chloro_3_methylphenol','4_chlorophenol','4_fluorophenol'};
flag_hlf = 0;

%%%%%%%%%%%%%%%%%%%

sulfides = {'di_isopropyl_sulfide','di_n_propyl_sulfide','diethyl_disulfide','diethyl_sulfide','dimethyl_disulfide',...
            'dimethyl_sulfide','hydrogen_sulfide','methyl_ethyl_sulfide','phenyl_methyl_sulfide'};
flag_slf = 0;

%%%%%%%%%%%%%%%%%%%

ether_glycols = {'12_ethanediol','2_ethoxyethanol','2_methoxyethanol','2_propoxyethanol','2_butoxyethanol'};
flag_etg = 0;

%%%%%%%%%%%%%%%%%%%

aromatic_ethers = {'anisole','diethoxymethoxybenzene','trimethoxymethylbenzene','phenyl_trifluoroethyl_ether'};
flag_are = 0;

%%%%%%%%%%%%%%%%%%%

thiols = {'thiophenol','ethanethiol','methanethiol','n_butanethiol','n_propanethiol'};
flag_thi = 0;

%%%%%%%%%%%%%%%%%%%

test = {'N_methyl_N__222_trifluoroethyl__aniline','phenyl_trifluoroethyl_ether','m_bis_trifluoromethyl__benzene'};
flag_tst = 1;

%%%%%%%%%%%%%%%%%%%
% This routine will plot all of the calulated vs. experimental results if 
% plot is set to 1.  Futhermore, it creates structures that contain all of
% the results for calculated solvation free energies as a function of each 
% solvent and again as a function of solute
for i = 1:length(solvents)
           
    fid = fopen(['mnsol/',lower(solvents{i}),'.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};
    allCases = 1:length(mol_list);
    
    [~,m] = ismember(testset,mol_list);
    if flag_testset    
        allCases(m) = 0;
    end
    
    [~,est] = ismember(esters,mol_list);
    if flag_est
        est(est==0) = [];
        allCases(est) = 0;
    end
    
    [~,eth] = ismember(ethers,mol_list);
    if flag_eth
        eth(eth==0) = [];
        allCases(eth) = 0;
    end
    
    [~,alc] = ismember(alcohols,mol_list);
    if flag_alc
        alc(alc==0) = [];
        allCases(alc) = 0;
    end
    
    [~,amd] = ismember(amides,mol_list);
    if flag_amd
        amd(amd==0) = [];
        allCases(amd) = 0;
    end
    
    [~,am1] = ismember(amines_1,mol_list);
    if flag_am1
        am1(am1==0) = [];
        allCases(am1) = 0;
    end
    
    [~,am2] = ismember(amines_2,mol_list);
    if flag_am2
        am2(am2==0) = [];
        allCases(am2) = 0;
    end
    
    [~,am3] = ismember(amines_3,mol_list);
    if flag_am3
        am3(am3==0) = [];
        allCases(am3) = 0;
    end
    
    [~,kam] = ismember(ketone_amines,mol_list);
    if flag_kam
        kam(kam==0) = [];
        allCases(kam) = 0;
    end
    
    [~,ham] = ismember(halo_amines,mol_list);
    if flag_ham
        ham(ham==0) = [];
        allCases(ham) = 0;
    end
    
    [~,nit] = ismember(nitro_compounds,mol_list);
    if flag_nit
        nit(nit==0) = [];
        allCases(nit) = 0;
    end
    
    [~,ntr] = ismember(nitriles,mol_list);
    if flag_ntr
        ntr(ntr==0) = [];
        allCases(ntr) = 0;
    end
    
    [~,phn] = ismember(phenols,mol_list);
    if flag_phn
        phn(phn==0) = [];
        allCases(phn) = 0;
    end
    
    [~,acd] = ismember(acids,mol_list);
    if flag_acd
        acd(acd==0) = [];
        allCases(acd) = 0;
    end
    
    [~,ktn] = ismember(ketones,mol_list);
    if flag_ktn
        ktn(ktn==0) = [];
        allCases(ktn) = 0;
    end
    
    if flag_hlk
    [~,hlk] = ismember(haloalkanes,mol_list);
        hlk(hlk==0) = [];
        allCases(hlk) = 0;
    end
    
    [~,hle] = ismember(haloalkenes,mol_list);
    if flag_hle 
        hle(hle==0) = [];
        allCases(hle) = 0;
    end
    
    [~,hlf] = ismember(halo_phenols,mol_list);
    if flag_hlf 
        hlf(hlf==0) = [];
        allCases(hlf) = 0;
    end
    
    [~,hbz] = ismember(halobenzenes,mol_list);
    if flag_hbz   
        hbz(hbz==0) = [];
        allCases(hbz) = 0;
    end
    
    [~,ahc] = ismember(aromatic_hydrocarbons,mol_list);
    if flag_ahc
        ahc(ahc==0) = [];
        allCases(ahc) = 0;
    end
    
    [~,lka] = ismember(alkanes,mol_list);
    if flag_lka
        lka(lka==0) = [];
        allCases(lka) = 0;
    end
    
    [~,lke] = ismember(alkenes,mol_list);
    if flag_lke
        lke(lke==0) = [];
        allCases(lke) = 0;
    end
    
    [~,slf] = ismember(sulfides,mol_list);
    if flag_slf
        slf(slf==0) = [];
        allCases(slf) = 0;
    end
    
    [~,ald] = ismember(aldehydes,mol_list);
    if flag_ald
        ald(ald==0) = [];
        allCases(ald) = 0;
    end
    
    [~,eam] = ismember(ether_amines,mol_list);
    if flag_eam
        eam(eam==0) = [];
        allCases(eam) = 0;
    end
    
    [~,etg] = ismember(ether_glycols,mol_list);
    if flag_etg
        etg(etg==0) = [];
        allCases(etg) = 0;
    end
    
    [~,are] = ismember(aromatic_ethers,mol_list);
    if flag_are
        are(are==0) = [];
        allCases(are) = 0;
    end
    
    [~,thi] = ismember(thiols,mol_list);
    if flag_thi
        thi(thi==0) = [];
        allCases(thi) = 0;
    end
    
    [~,tst] = ismember(test,mol_list);
    if flag_tst
        tst(tst==0) = [];
        allCases(tst) = 0;
    end
    
    temp = load(['Run',solvents{i}]);
%     solute_errors(i,:) = temp.errfinal(n);
%    sorted_errors = sort(abs(temp.errfinal));
    results(i) = struct('Solvent', solvents{i},'Num_Solutes',round(length(mol_list),3),'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE),...
                        'Mean_Abs_error',mean(abs(temp.errfinal)));
    if ploton
        figure
        flag_pred = 0;
        if flag_pred
            hold on
            plot(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)),'bo','markers',2)
        end
        set(gca,'FontSize',15)
        
        if flag_testset
            hold on
            plot(temp.refE(m),temp.calcE(m),'gs','markers',12)
        end
        
        if flag_est
            hold on
            plot(temp.refE(est),temp.calcE(est),'r+','markers',12)
        end
        
        if flag_eth
            hold on
            plot(temp.refE(eth),temp.calcE(eth),'rd','markers',12)
        end
        
        if flag_alc
            hold on
            plot(temp.refE(alc),temp.calcE(alc),'ro','markers',12)
        end
        
        if flag_amd
            hold on
            plot(temp.refE(amd),temp.calcE(amd),'b+','markers',12)
        end
        
        if flag_am1
            hold on
            plot(temp.refE(am1),temp.calcE(am1),'b*','markers',12)
        end
        
        if flag_am2
            hold on
            plot(temp.refE(am2),temp.calcE(am2),'bd','markers',12)
        end
        
        if flag_am3
            hold on
            plot(temp.refE(am3),temp.calcE(am3),'bx','markers',12)
        end
        
        if flag_ham
            hold on
            plot(temp.refE(ham),temp.calcE(ham),'b<','markers',12)
        end
        
        if flag_eam
            hold on
            plot(temp.refE(eam),temp.calcE(eam),'b>','markers',12)
        end
        
        if flag_kam
            hold on
            plot(temp.refE(kam),temp.calcE(kam),'b^','markers',12)
        end
        
        if flag_nit
            hold on
            plot(temp.refE(nit),temp.calcE(nit),'bo','markers',12)
        end
        
        if flag_ntr
            hold on
            plot(temp.refE(ntr),temp.calcE(ntr),'bd','markers',12)
        end
        
        if flag_acd    
            hold on
            plot(temp.refE(acd),temp.calcE(acd),'rx','markers',12)
        end
        
        if flag_ktn
            hold on
            plot(temp.refE(ktn),temp.calcE(ktn),'rs','markers',12)
        end
        
        if flag_hlk
            hold on
            plot(temp.refE(hlk),temp.calcE(hlk),'m*','markers',12)
        end
        
        if flag_hle
            hold on
            plot(temp.refE(hle),temp.calcE(hle),'mx','markers',12)
        end
        
        if flag_hlf
            hold on
            plot(temp.refE(hlf),temp.calcE(hlf),'ks','markers',12)
        end
        
        if flag_phn
            hold on
            plot(temp.refE(phn),temp.calcE(phn),'kd','markers',12)
        end
        
        if flag_hbz
            hold on
            plot(temp.refE(hbz),temp.calcE(hbz),'m+','markers',12)
        end
        
        if flag_ahc
            hold on
            plot(temp.refE(ahc),temp.calcE(ahc),'k*','markers',12)
        end
        
        if flag_ald
            hold on
            plot(temp.refE(ald),temp.calcE(ald),'r<','markers',12)
        end
        
        if flag_lka
            hold on
            plot(temp.refE(lka),temp.calcE(lka),'ko','markers',12)
        end
        
        if flag_lke
            hold on
            plot(temp.refE(lke),temp.calcE(lke),'kx','markers',12)
        end
        
        if flag_slf
            hold on
            plot(temp.refE(slf),temp.calcE(slf),'co','markers',12)
        end
        
        if flag_etg
            hold on
            plot(temp.refE(etg),temp.calcE(etg),'r>','markers',12)
        end
        
        if flag_are
            hold on
            plot(temp.refE(are),temp.calcE(are),'r<','markers',12)
        end
        
        if flag_thi
            hold on
            plot(temp.refE(thi),temp.calcE(thi),'k<','markers',12)
        end
        
        if flag_tst
            hold on
            plot(temp.refE(tst),temp.calcE(tst),'g*','markers',12)
        end
        
        A = [flag_pred flag_testset flag_est flag_eth flag_alc flag_amd flag_am1 flag_am2 flag_am3 flag_ham flag_eam flag_kam flag_nit flag_ntr flag_acd flag_ktn flag_hlk flag_hle flag_hlf flag_phn flag_hbz flag_ahc flag_ald flag_lka flag_lke flag_slf flag_etg flag_are flag_thi flag_tst];
        Leg = {'Predictions','Training Set','Esters','Ethers','Alcohols','Amides','Primary amines','Secondary amines','Teritiary amines','Halo amines','Ether + amines','Ketone + amines','Nitro compounds','Nitriles','Acids','Ketones','Haloalkanes','Haloalkenes','Halo phenols','Phenol derivatives','Halobenzenes','Aromatic hydrocarbons','Aldehydes','Alkanes','Alkenes','Sulfides','Ether glycols','Aromatic ethers','Thiols','test'};
        Leg(A==0)=[];
        minax = round(min(min(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)))));
        maxax = round(max(max(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)))));
        axis([minax-2 maxax+2 minax-2 maxax+2]);
        foo = refline(1,0);
        set(foo,'Linewidth',2,'color','k');
        xlabel(['\Delta','G_{expt}^{solv, ',solvents{i},'}',' (kcal.mol^{-1})'])
        ylabel(['\Delta','G_{calc}^{solv, ',solvents{i},'}',' (kcal.mol^{-1})'])
        legend(Leg,'Location','southeast')
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

