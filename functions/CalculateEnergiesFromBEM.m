function [calculatedE, referenceE ] = CalculateEnergiesFromBEM(Params)
global UsefulConstants ProblemSet 

% define empty vectors so we can extend them easily
calculatedE  = [];
referenceE = [];

% how many problem geometries do we have ?  note that this can be
% DIFFERENT FROM the number of overall "BEM calculations" we do.

numProblems = length(ProblemSet);

for i=1:numProblems
  curProblem = ProblemSet(i);
  [newCalculatedE, newReferenceE]  = calculateProblem(curProblem, Params);

%  ProblemSet(i).bemPcm = [];
%  ProblemSet(i).bemYoonStern = [];
%  ProblemSet(i).asymBemPcm = [];
%  ProblemSet(i).uninitialized = 1;
  
  % append calculated and reference results to total list  
  calculatedE = [calculatedE; newCalculatedE];
  referenceE = [referenceE; newReferenceE];
end

% there's no need for a trailing "end" at the end of a Matlab
% function, when it's in a .m file by itself!