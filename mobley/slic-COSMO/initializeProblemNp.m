function problemIndex =  initializeProblemNp(problem)
global UsefulConstants ProblemSetNp

% this use of problemIndex is SUBTLE! please ask if the motivation/use
% is not clear, after you think about it for a while.
problemIndex = -1;
for i=1:length(ProblemSetNp)
  if ProblemSetNp(i).index == problem.index
    problemIndex = i;
  end
end
if problemIndex==-1
  fprintf('Error in initializeProblemNp: Could not identify problem in ProblemSetNp!\n');
  keyboard
end

% Key detail: notice that we use ProblemSetNp(problemIndex) throughout,
% INSTEAD OF "problem" itself.

if ProblemSetNp(problemIndex).uninitializedNp

  % useful constants
  epsIn  = UsefulConstants.epsIn;    % e.g., 1
  epsOut = UsefulConstants.epsOut;   % e.g., 80
  kappa  = UsefulConstants.kappa;    % e.g., 0.0 for non-ionic
  
% unset uninitialized field, so it's understood that the next
  % "test" calculation won't repeat these expensive (and memory
  % hungry!!) steps
  ProblemSetNp(problemIndex).uninitialized = 0;
end

