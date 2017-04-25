common_solutes = {'ethanol','butanone','toluene','n_octane','nitromethane',...
            '14_dioxane','phenol','acetic_acid','methanol','propanoic_acid',...
            'propan_1_ol','ethyl_acetate','aniline','ethylamine','n_butyl_acetate',...
            'methyl_acetate','hexan_1_ol','pentan_1_ol','n_propyl_acetate','butan_1_ol',...
            'heptan_1_ol','methyl_pentanoate','p_cresol','methyl_propanoate',...
            'o_cresol','propanone','pyridine'};

num_solutes = 12;
test_set = randsample(common_solutes,num_solutes);
save('testSet','test_set');