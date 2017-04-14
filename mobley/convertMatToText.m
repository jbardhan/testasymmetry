solvents = {'Water', 'Octanol', 'Hexadecane', 'Chloroform', 'Cyclohexane',...
            'Carbontet', 'Hexane', 'Toluene', 'Xylene'}; 
        
for i = 1:length(solvents)
   data = load(['Run',solvents{i}]);
   fid = fopen([solvents{i},'-Data.txt'],'w');
   fprintf(fid,'calcE, errfinal, es, np, refE\n');
   fprintf(fid,'%f, %f, %f, %f, %f \n',[calcE,errfinal,es,np,refE].');
   fclose(fid);
end