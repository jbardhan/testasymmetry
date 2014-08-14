scaleFactor = 1.02;

fid = fopen('radii_roux.siz');
fgetl(fid);
count = 1;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  [x1, tline] = strtok(tline); 
  if length(x1) == 0, break, end
  [x2, tline] = strtok(tline);
  [x3, tline] = strtok(tline);
  atomname{count} = x1;
  if length(x3) == 0
    resname{count} = '';
    charge{count} = str2double(x2);
  else
    resname{count} = x2;
    charge{count} = str2double(x3);
  end

  count = count + 1;
end
fclose(fid);

for count = 1:length(atomname)
  newcharge{count} = scaleFactor * charge{count};
end


fid = fopen('roux_ses.siz','w');
fprintf(fid,'atom__res_radius_\n');
for count=1:length(atomname)
  spacestr = repmat(' ',1,5-length(atomname{count}));
  atm = sprintf('%s%s',spacestr,atomname{count});
  spacestr = repmat(' ',1,5-length(resname{count}));
  res = sprintf(' %s%s',resname{count},spacestr);
  fprintf(fid,'%s%s%5.4f\n',atm,res,newcharge{count});
end
fclose(fid);