function crg = readCrg(crgFile)
% this file expects that the CRG and PDB have atoms listed in
% exactly the same order!
fid = fopen(crgFile,'r');
line = fgetl(fid);

n = 1;
while 1
  line = fgetl(fid);
  if ~ischar(line), break, end
  atomname(n) = cellstr(strtrim(line(1:6)));
  resname(n)  = cellstr(strtrim(line(7:9)));
  resnum(n)   = cellstr(strtrim(line(10:13)));
  segid(n)    = cellstr(strtrim(line(14)));
  charge(n)   = str2num(strtrim(line(15:22)));
  n = n+1;
end
crg = charge';
return;
crg = struct('atomname',atomname,...
	     'resname',resname,...
	     'resnum', resnum,...
	     'segid', segid,...
	     'charge',charge);
fclose(fid);