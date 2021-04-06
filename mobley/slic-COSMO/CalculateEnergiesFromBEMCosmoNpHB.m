function [calculatedE, referenceE,nonpolar,...
          hb,disp,disp_sl_sl,disp_sv_sl,disp_sv_sv,cav,comb] = CalculateEnergiesFromBEMCosmoNpHB(Params)
global UsefulConstants ProblemSetNp saveMemory writeLogfile logfileName

% define empty vectors so we can extend them easily
calculatedE  = [];
referenceE = [];
nonpolar = [];
hb = [];
disp = [];
disp_sl_sl= [];
disp_sv_sl= [];
disp_sv_sv= [];
cav = [];
comb = [];

% how many problem geometries do we have ?  note that this can be
% DIFFERENT FROM the number of overall "BEM calculations" we do.

numProblems = length(ProblemSetNp);

for i=1:numProblems
    i
  curProblem = ProblemSetNp(i);
  [newCalculatedE, newReferenceE,newNonpolar,newHB,newDisp,...
   newDisp_sl_sl,newDisp_sv_sl,newDisp_sv_sv,newCav,newComb] = ...
   calculateProblemCosmoNpHB(curProblem, Params);

if saveMemory
  ProblemSetNp(i).uninitialized = 1;
end

if writeLogfile
  fid = fopen(logfileName,'a');
  fprintf(fid,'%s,%f,%f,%f\n',curProblem.name,newReferenceE(1), ...
	  newCalculatedE(1),newNonpolar(1));
  fclose(fid);
end

  % append calculated and reference results to total list  
  calculatedE = [calculatedE; newCalculatedE];
  referenceE = [referenceE; newReferenceE];
  nonpolar = [nonpolar; newNonpolar];
  hb = [hb; newHB];
  disp = [disp; newDisp];
  disp_sl_sl = [disp_sl_sl; newDisp_sl_sl];
  disp_sv_sl = [disp_sv_sl; newDisp_sv_sl];
  disp_sv_sv = [disp_sv_sv; newDisp_sv_sv];
  cav = [cav; newCav];
  comb = [comb; newComb];
  
end

% there's no need for a trailing "end" at the end of a Matlab
% function, when it's in a .m file by itself!