function [matched_list,table] = calcTransferEnergies(file1,file2)

[set1_list,set1_ref,set1_calc,set1_es,set1_np]=...
    textread(file1,'%s %f %f %f %f','delimiter',',');
[set2_list,set2_ref,set2_calc,set2_es,set2_np]=...
    textread(file2,'%s %f %f %f %f','delimiter',',');

counter = 1;
ref1=[];
ref2=[];
for i=1:length(set1_list)
  foo = strcmp(set2_list,set1_list{i});
  index = find(foo);
  if length(index) == 0
%    fprintf('No match for solute %s in solvent 2\n',set1_list{i});
  elseif length(index) == 1
    matched_list{counter}=set1_list{i};

    ref1(counter) = set1_ref(i);
    calc1(counter) = set1_calc(i);
    es1(counter) = set1_es(i);
    np1(counter) = set1_np(i);
    
    ref2(counter) = set2_ref(index);
    calc2(counter) = set2_calc(index);
    es2(counter) = set2_es(index);
    np2(counter) = set2_np(index);

    counter = counter+1;

  else
    fprintf('Multiple matches for solute %s in solvent 2?',set1_list{i});
    keyboard
  end
end

table = struct('ref1',ref1,'calc1',calc1,...
	       'es1',es1,'np1',np1,...
	       'ref2',ref2,'calc2',calc2,...
	       'es2',es2,'np2',np2);
