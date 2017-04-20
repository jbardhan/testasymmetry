function writeDat(fileID,dataStructure)
fnames = fieldnames(dataStructure);
first = sprintf(fnames{1});

if strcmp(first,'Solute_Name')
    solute = 1;
    solvent = 0;
else
    solvent = 1;
    solute = 0;
end

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
if solute
    fprintf(fid,'\\begin{tabular}{|l|l|l|} \n');
elseif solvent
    fprintf(fid,'\\begin{tabular}{|l|l|l|l|l|l|l|l|} \n');
end
fprintf(fid,'\\hline  \n');
if solute
    fprintf(fid,'\\textbf{Solute Name}&\\textbf{RMS Error}&\\textbf{Mean Absolute Error} \\\\ \\hline\n');
elseif solvent
    fprintf(fid,'\\textbf{Solvent Name}&\\textbf{Number of Solutes}&\\textbf{RMS Error}&\\textbf{Mean Absolute Error}&\\textbf{RMS$_{T}$}&\\textbf{Absolute$_{T}$}&\\textbf{RMS$_{C}$}&\\textbf{Absolute$_{C}$} \\\\ \\hline\n'); 
end
fprintf(fid,'\\hline \n');
for i = 1:length(dataStructure)
    if strcmp(first,'Solvent')
        fprintf(fid,'\\textbf{%s}&%3.0f&%f&%f&%f&%f&%f&%f \\\\ \\hline \n',replace(string(dataStructure(i).Solvent),'_','\_'),dataStructure(i).Num_Solutes,dataStructure(i).RMS,dataStructure(i).Mean_Abs_error,dataStructure(i).RMS_Training,dataStructure(i).Mean_Abs_error_Training,dataStructure(i).RMS_Cons,dataStructure(i).Mean_Abs_error_Cons); 
    else
        fprintf(fid,'\\textbf{%s}&%f&%f \\\\ \\hline \n',replace(string(dataStructure(i).Solute_Name),'_','\_'),dataStructure(i).RMS,dataStructure(i).Mean_Abs_error);
    end
end
fprintf(fid,'\\end{tabular} \n');
fprintf(fid,'\\end{center} \n');
fprintf(fid,'\\end{document} \n');
fclose(fid);