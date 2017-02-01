function [indicesIn, indicesOut] = crosscheck(list1,list2);

counter = 1;
indices = zeros(length(list1),1);
for i=1:length(list2)
  foo = strcmp(list1,list2{i});
  index = find(foo);
  if length(index) == 0
%    fprintf('No match for solute %s in solvent 2\n',set1_list{i});
  elseif length(index) == 1
    indices(index) = 1;
  end
end

indicesIn = find(indices==1);
indicesOut=find(indices==0);
