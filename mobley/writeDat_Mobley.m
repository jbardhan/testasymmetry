function writeDat_Mobley(fileID,data)
    
    mol_list=data.mol_list;
    index=data.index;
    dg_ref=data.refE(index,:);
    dg_calc_slic=data.calcE(index,:);
    ds_calc_slic=data.dsvec;
    cp_calc_slic=data.cpvec;

    
    dg_calc_mobley=data.calc_mobley;
    dg_rms_slic_mobley=data.dg_rms_298;
    dg_rms_mobley=data.dg_rms_298_MD;
    
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
    fprintf(fid,'\\begin{longtable}{l|c|c|c} \n');
    fprintf(fid,'\\textbf{Solute Name}&\\textbf{$\\Delta G$(kcal/mol)}&\\textbf{$\\Delta S$(cal/mol$^\\circ$K)}&\\textbf{$Cp$(cal/mol$^\\circ$K)} \\\\ \n');
    fprintf(fid,'\\hline \\\\ \n');
    for i=1:length(mol_list)
        fprintf(fid,'%s & %.2f(\\textit{%.2f}) & %.2f & %.2f \\\\ \n',replace(string(mol_list(i)),'_','\_'),dg_calc_slic(i),dg_ref(i),ds_calc_slic(i),cp_calc_slic(i));
    end
    
    
    fprintf(fid,'\\end{longtable} \n');
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