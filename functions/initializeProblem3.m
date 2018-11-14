function problemIndex =  initializeProblem3(problem)
global UsefulConstants3 ProblemSet3

% this use of problemIndex is SUBTLE! please ask if the motivation/use
% is not clear, after you think about it for a while.
problemIndex = -1;
for i=1:length(ProblemSet3)
  if ProblemSet3(i).index == problem.index
    problemIndex = i;
  end
end
if problemIndex==-1
  fprintf('Error in initializeProblem: Could not identify problem in ProblemSet3!\n');
  keyboard
end

% Key detail: notice that we use ProblemSet(problemIndex) throughout,
% INSTEAD OF "problem" itself.

if ProblemSet3(problemIndex).uninitialized

  % useful constants
  epsIn1  = UsefulConstants3.epsIn1;    % e.g., 1
  epsIn2  = UsefulConstants3.epsIn2;    % e.g., 1
  epsOut = UsefulConstants3.epsOut;   % e.g., 80
  kappa  = UsefulConstants3.kappa;    % e.g., 0.0 for non-ionic
  
  % loadSternSrfIntoPanels goes through our custom, BEM for
  % molecules .srf file.  It's worth reading to understand about
  % the file format, and the complications of dealing with meshes.
  ProblemSet3(problemIndex).srf1SternData = ...
      loadSternSrfIntoPanels(ProblemSet3(problemIndex).srfFile1);
  
  ProblemSet3(problemIndex).srf2SternData = ...
      loadSternSrfIntoPanels(ProblemSet3(problemIndex).srfFile2);

  ProblemSet3(problemIndex).srf3SternData = ...
      loadSternSrfIntoPanels(ProblemSet3(problemIndex).srfFile3);

  % the operators from makePanelBemEcfQualMatrices are used for calculating the
  % electric field just inside the interface
  ProblemSet3(problemIndex).bemPcm1 = ...
      makePanelBemEcfQualMatrices(ProblemSet3(problemIndex).srf1SternData.dielBndy(1), ...
				  ProblemSet3(problemIndex).pqrData1, epsIn1, epsOut);

  ProblemSet3(problemIndex).bemPcm2 = ...
      makePanelBemEcfQualMatrices(ProblemSet3(problemIndex).srf2SternData.dielBndy(1), ...
				  ProblemSet3(problemIndex).pqrData2, epsIn2, epsOut);

  % here is the main integral equation operator, preconditioner, etc.
  ProblemSet3(problemIndex).bemYoonStern = makePanelBemSternMatrices3(ProblemSet3(problemIndex).srf1SternData, ...
                          ProblemSet3(problemIndex).srf2SternData, ProblemSet3(problemIndex).srf3SternData, ...
						  ProblemSet3(problemIndex).pqrData1,ProblemSet3(problemIndex).pqrData2, epsIn1, epsIn2, ...
						  epsOut, kappa);

  % here are the asymmetry-related perturbations to the problem
  % that are independent of the electric field.  have a look.
  ProblemSet3(problemIndex).asymBemPcm1 = makePanelAsymEcfCollocMatrices(ProblemSet3(problemIndex).srf1SternData.dielBndy(1), ...
						  ProblemSet3(problemIndex).bemPcm1, ...
						  ProblemSet3(problemIndex).pqrData1);
                      
  ProblemSet3(problemIndex).asymBemPcm2 = makePanelAsymEcfCollocMatrices(ProblemSet3(problemIndex).srf2SternData.dielBndy(1), ...
						  ProblemSet3(problemIndex).bemPcm2, ...
						  ProblemSet3(problemIndex).pqrData2);

  % unset uninitialized field, so it's understood that the next
  % "test" calculation won't repeat these expensive (and memory
  % hungry!!) steps
  ProblemSet3(problemIndex).uninitialized = 0;
end

