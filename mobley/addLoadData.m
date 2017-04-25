function addLoadData(solvent, loadFile)
global LoadData

newIndex = length(LoadData)+1;
newproblem = struct('Solvent',solvent,'loadFile',loadFile);

if newIndex == 1
  LoadData = newproblem;
else
  LoadData(newIndex) = newproblem;
end
