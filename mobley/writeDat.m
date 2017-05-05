function writeDat(fileID,run_water_training)
    
    testset=run_water_training.testset;
    index_training=run_water_training.index;  % index of 298K =24.85C in the temp vector
    dg_ref_training=run_water_training.refE(index_training,:);   % expaerimental Delta_G of the training set in kcal/mol
    ds_ref_training=run_water_training.refS;  % expaerimentalDelta S of the training set in cal/mol/K
    cp_ref_training=run_water_training.refCP;  % expaerimentalDelta S of the training set in cal/mol/K
    dg_calc_training=run_water_training.calcE(index_training,:); % calculated Delta_G of the training set in kcal/mol
    ds_calc_training=run_water_training.calcS; % calculated S of the training set in cal/mol/K
    cp_calc_training=run_water_training.calcCP; % calculated S of the training set in cal/mol/K
    
    fid = fopen(fileID,'w');
    fprintf(fid,'\\documentclass{amsart} \n');
    fprintf(fid,'\\usepackage{float} \n');
    fprintf(fid,'\\usepackage{caption} \n');
    fprintf(fid,'\\usepackage{longtable} \n');
    fprintf(fid,'\\usepackage{multirow} \n');
    fprintf(fid,'\\usepackage{array} \n');
    fprintf(fid,'\\newcolumntype{C}{>{\\centering\\arraybackslash}p{5em}}  \n');
    fprintf(fid,'\\addtolength{\\oddsidemargin}{-1.6in} \n');
    fprintf(fid,'\\addtolength{\\evensidemargin}{-1.6in} \n');
    fprintf(fid,'\\addtolength{\\topmargin}{-1in} \n');


    fprintf(fid,'\\begin{document} \n');
    fprintf(fid,'\\captionof{table}{} \n');
    
    fprintf(fid,'\\begin{tabular}{l|C|C|C|C|C|C} \n');
    fprintf(fid,'\\centering \n');
    
    fprintf(fid,'\\textbf{Solute Name}  &   \\multicolumn{2}{c|}{\\textbf{$\\Delta G^{water}$}($\\frac{kcal}{mol}$)}  &   \\multicolumn{2}{c|}{\\textbf{$\\Delta S^{water}$}($\\frac{cal}{mol.K}$)} &   \\multicolumn{2}{c}{\\textbf{$Cp^{water}$}($\\frac{cal}{mol.K}$)} \\\\ \n');
    fprintf(fid,'& SLIC & Expt. & SLIC & Expt. & SLIC & Expt. \\\\ \n');
    fprintf(fid,'\\hline \n');
    
    for i=1:length(run_water_training.testset)
        fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\ \n',replace(string(testset(i)),'_','\_'),dg_calc_training(i),dg_ref_training(i),ds_calc_training(i),ds_ref_training(i),cp_calc_training(i),cp_ref_training(i));
    end
    
    fprintf(fid,'\\end{tabular} \n');
    fprintf(fid,'\\end{document} \n');
    fclose(fid);


end