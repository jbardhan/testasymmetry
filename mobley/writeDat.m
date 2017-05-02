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
    fprintf(fid,'\\addtolength{\\oddsidemargin}{-1.6in} \n');
    fprintf(fid,'\\addtolength{\\evensidemargin}{-1.6in} \n');
    fprintf(fid,'\\addtolength{\\topmargin}{-1in} \n');
    fprintf(fid,'\\begin{document} \n');
    
    
    fprintf(fid,'\\begin{center}\n');
    fprintf(fid,'\\captionof{table}{} \n');   
    fprintf(fid,'\\begin{tabular}{l|c|c|c} \n');
    fprintf(fid,'\\textbf{Solute Name}&\\textbf{$\\Delta G$(kcal/mol)}&\\textbf{$\\Delta S$(cal/mol$^\\circ$K)}&\\textbf{$Cp$(cal/mol$^\\circ$K)} \\\\ \n');
    fprintf(fid,'\\hline \\\\ \n');
    for i=1:length(run_water_training.testset)
        fprintf(fid,'%s & %.2f(\\textit{%.2f}) & %.2f(\\textit{%.2f}) & %.2f(\\textit{%.2f}) \\\\ \n',replace(string(testset(i)),'_','\_'),dg_calc_training(i),dg_ref_training(i),ds_calc_training(i),ds_ref_training(i),cp_calc_training(i),cp_ref_training(i));
    end
    
    
    fprintf(fid,'\\end{tabular} \n');
    fprintf(fid,'\\end{center} \n');
    fprintf(fid,'\\end{document} \n');
    fclose(fid);

%     fprintf(fid,'\\hline \n');
%     for i = 1:length(dataStructure)
%         if strcmp(first,'Solvent')
%             fprintf(fid,'\\textbf{%s}&%3.0f&%f&%f&%f&%f&%f&%f \\\\ \\hline \n',replace(string(dataStructure(i).Solvent),'_','\_'),dataStructure(i).Num_Solutes,round(dataStructure(i).RMS,2),round(dataStructure(i).Mean_Abs_error,2),round(dataStructure(i).RMS_Training,2),round(dataStructure(i).Mean_Abs_error_Training,2),round(dataStructure(i).RMS_Cons,2),round(dataStructure(i).Mean_Abs_error_Cons,2)); 
%         else
%             fprintf(fid,'\\textbf{%s}&%f&%f \\\\ \\hline \n',replace(string(dataStructure(i).Solute_Name),'_','\_'),round(dataStructure(i).RMS,2),round(dataStructure(i).Mean_Abs_error,2));
%         end
%     end
%     fprintf(fid,'\\end{tabular} \n');
%     fprintf(fid,'\\end{center} \n');
%     fprintf(fid,'\\end{document} \n');
%     fclose(fid);
% elseif numel(size(dataStructure))==2
%     fid = fopen(fileID,'w');
%     fprintf(fid,'\\documentclass{amsart} \n');
%     fprintf(fid,'\\usepackage{float} \n');
%     fprintf(fid,'\\usepackage{caption} \n');
%     fprintf(fid,'\\usepackage{longtable} \n');
%     fprintf(fid,'\\addtolength{\\oddsidemargin}{-1.6in} \n');
%     fprintf(fid,'\\addtolength{\\evensidemargin}{-1.6in} \n');
%     fprintf(fid,'\\addtolength{\\topmargin}{-1in} \n');
%     fprintf(fid,'\\begin{document} \n');
%     fprintf(fid,'\\begin{center}\n');
%     fprintf(fid,'\\captionof{table}{} \n');
%     fprintf(fid,'\\begin{tabular}{l|l|l|l|l|l|l|l|l|l} \n');
%     fprintf(fid,strcat('&',join(solvents,'&'),' \\\\ \\hline \n'));
%     
%     for i = 1:length(solvents)
%         fprintf(fid,strcat(['\\textbf{',solvents{i},'}'],'&',join(string(round(dataStructure(i,:),2)),' &'),'\\\\ \\hline \n'));
%     end
%     
%     fprintf(fid,'\\end{tabular} \n');
%     fprintf(fid,'\\end{center} \n');
%     fprintf(fid,'\\end{document} \n');
%     fclose(fid);
% elseif numel(size(dataStructure))==3
%     dataSize = size(dataStructure);
%     numDocs = dataSize(3);
%     for j = 1:numDocs
%         temp1 = strsplit(fileID,'.');
%         temp2 = strcat(temp1(1),string(j));
%         fname = strcat(temp2,'.',temp1(2));
%         fid = fopen(fname,'w');
%         fprintf(fid,'\\documentclass{amsart} \n');
%         fprintf(fid,'\\usepackage{float} \n');
%         fprintf(fid,'\\usepackage{caption} \n');
%         fprintf(fid,'\\usepackage{longtable} \n');
%         fprintf(fid,'\\addtolength{\\oddsidemargin}{-1.6in} \n');
%         fprintf(fid,'\\addtolength{\\evensidemargin}{-1.6in} \n');
%         fprintf(fid,'\\addtolength{\\topmargin}{-1in} \n');
%         fprintf(fid,'\\begin{document} \n');
%         fprintf(fid,'\\begin{center}\n');
%         fprintf(fid,'\\captionof{table}{} \n');
%         fprintf(fid,'\\begin{tabular}{l|l|l|l|l|l|l|l|l|l} \n');
%         fprintf(fid,strcat('&',join(solvents,'&'),' \\\\ \\hline \n'));
% 
%         for i = 1:length(solvents)
%             fprintf(fid,strcat(['\\textbf{',solvents{i},'}'],'&',join(string(round(dataStructure(i,:,j),2)),' &'),'\\\\ \\hline \n'));
%         end
% 
%         fprintf(fid,'\\end{tabular} \n');
%         fprintf(fid,'\\end{center} \n');
%         fprintf(fid,'\\end{document} \n');
%         fclose(fid);
%     end
end