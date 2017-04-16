solvents = {'Carbontet','Chloroform','Cyclohexane','Hexadecane', ...
	    'Hexane','Octanol','Toluene','Water','Xylene'};

for i=1:length(solvents)
  filename = sprintf('Run%s',solvents{i});
  load(filename);
  allResults = [refE; calcE];
  figure;
  plot([refE],[calcE],'bx','linewidth',2);
  hold on;
  set(gca,'FontSize',16);
  plot([min(allResults) max(allResults)],[min(allResults) max(allResults)],'k','linewidth',2);
  hold on
  plot([min(allResults) max(allResults)],1+[min(allResults) max(allResults)],'r--','linewidth',2);
  plot([min(allResults) max(allResults)],-1+[min(allResults) max(allResults)],'r--','linewidth',2);
  title(solvents{i});
end
