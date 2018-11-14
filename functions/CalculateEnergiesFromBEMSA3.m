function [calculatedE, referenceE,electrostatic,nonpolar ] = CalculateEnergiesFromBEMSA3(Params)
global UsefulConstants3 ProblemSet3 saveMemory writeLogfile logfileName

% define empty vectors so we can extend them easily
calculatedE  = [];
referenceE = [];
electrostatic  =[];
nonpolar = [];
% how many problem geometries do we have ?  note that this can be
% DIFFERENT FROM the number of overall "BEM calculations" we do.

numProblems = length(ProblemSet3);

for i=1:numProblems
    i
  curProblem = ProblemSet3(i);
  [newCalculatedE, newReferenceE,newElectrostatic,newNonpolar]  = calculateProblemSA3(curProblem, Params);

if saveMemory
  ProblemSet3(i).bemPcm1 = [];
  ProblemSet3(i).bemPcm2 = [];
  ProblemSet3(i).bemYoonStern = [];
  ProblemSet3(i).asymBemPcm1 = [];
  ProblemSet3(i).asymBemPcm2 = [];
  ProblemSet3(i).uninitialized = 1;
end

if writeLogfile
  fid = fopen(logfileName,'a');
  fprintf(fid,'%s,%s,%f,%f,%f,%f\n',curProblem.name1,curProblem.name2,newReferenceE(1), ...
	  newCalculatedE(1),newElectrostatic(1),newNonpolar(1));
  fclose(fid);
end

  % append calculated and reference results to total list  
  calculatedE = [calculatedE; newCalculatedE];
  referenceE = [referenceE; newReferenceE];
  electrostatic= [electrostatic; newElectrostatic];
  nonpolar = [nonpolar; newNonpolar];
end

% there's no need for a trailing "end" at the end of a Matlab
% function, when it's in a .m file by itself!