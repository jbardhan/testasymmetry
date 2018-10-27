function problemIndex =  initializeProblem2(problem)
global UsefulConstants2 ProblemSet2

% this use of problemIndex is SUBTLE! please ask if the motivation/use
% is not clear, after you think about it for a while.
problemIndex = -1;
for i=1:length(ProblemSet2)
  if ProblemSet2(i).index == problem.index
    problemIndex = i;
  end
end
if problemIndex==-1
  fprintf('Error in initializeProblem: Could not identify problem in ProblemSet2!\n');
  keyboard
end

% Key detail: notice that we use ProblemSet(problemIndex) throughout,
% INSTEAD OF "problem" itself.

if ProblemSet2(problemIndex).uninitialized

  % useful constants
  epsIn1  = UsefulConstants2.epsIn1;    % e.g., 1
  epsIn2  = UsefulConstants2.epsIn2;    % e.g., 1
  epsOut = UsefulConstants2.epsOut;   % e.g., 80
  kappa  = UsefulConstants2.kappa;    % e.g., 0.0 for non-ionic
  
  % loadSternSrfIntoPanels goes through our custom, BEM for
  % molecules .srf file.  It's worth reading to understand about
  % the file format, and the complications of dealing with meshes.
  ProblemSet2(problemIndex).srf1SternData = ...
      loadSternSrfIntoPanels(ProblemSet2(problemIndex).srfFile1);
  
  ProblemSet2(problemIndex).srf2SternData = ...
      loadSternSrfIntoPanels(ProblemSet2(problemIndex).srfFile2);

  % the operators from makePanelBemEcfQualMatrices are used for calculating the
  % electric field just inside the interface
  ProblemSet2(problemIndex).bemPcm1 = ...
      makePanelBemEcfQualMatrices(ProblemSet2(problemIndex).srf1SternData.dielBndy(1), ...
				  ProblemSet2(problemIndex).pqrData1, epsIn1, epsOut);

  ProblemSet2(problemIndex).bemPcm2 = ...
      makePanelBemEcfQualMatrices(ProblemSet2(problemIndex).srf2SternData.dielBndy(1), ...
				  ProblemSet2(problemIndex).pqrData2, epsIn2, epsOut);

  % here is the main integral equation operator, preconditioner, etc.
  ProblemSet2(problemIndex).bemYoonStern = makePanelBemSternMatrices2(ProblemSet2(problemIndex).srf1SternData, ProblemSet2(problemIndex).srf2SternData, ...
						  ProblemSet2(problemIndex).pqrData1,ProblemSet2(problemIndex).pqrData2, epsIn1, epsIn2, ...
						  epsOut, kappa);

  % here are the asymmetry-related perturbations to the problem
  % that are independent of the electric field.  have a look.
  ProblemSet2(problemIndex).asymBemPcm1 = makePanelAsymEcfCollocMatrices(ProblemSet2(problemIndex).srf1SternData.dielBndy(1), ...
						  ProblemSet2(problemIndex).bemPcm1, ...
						  ProblemSet2(problemIndex).pqrData1);
                      
  ProblemSet2(problemIndex).asymBemPcm2 = makePanelAsymEcfCollocMatrices(ProblemSet2(problemIndex).srf2SternData.dielBndy(1), ...
						  ProblemSet2(problemIndex).bemPcm2, ...
						  ProblemSet2(problemIndex).pqrData2);

  % unset uninitialized field, so it's understood that the next
  % "test" calculation won't repeat these expensive (and memory
  % hungry!!) steps
  ProblemSet2(problemIndex).uninitialized = 0;
end

