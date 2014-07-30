function [crd, crg, atype, resname,segid,resid,atomnamelist] = readpdb(pdbfile, Nmax)

% get coordinates and atomic form factor from pdb file
%
% crd = atomic coordinates (N*3)
% atype = atom type (N*1)

crd = zeros(Nmax,3);  atype = char(zeros(Nmax,1));
crg = zeros(Nmax,1);
fid = fopen(pdbfile, 'r');

n = 0;
while 1
  line = fgetl(fid);
  if ~ischar(line), break, end
  if length(line) > 6 ...
     && ( strcmp(line(1:4),'ATOM') || strcmp(line(1:6),'HETATM') )
    n = n + 1;
    atomname = strtrim(line(13:16));
    resname(n) = cellstr(strtrim(line(18:21)));
	 segid(n)   = cellstr(strtrim(line(22:23)));
	 resid(n)   = cellstr(strtrim(line(24:28)));
	 atomnamelist(n) = cellstr(atomname);
    if strcmp(atomname, 'FE')
      atype(n) = 'f';
    elseif '0' <= atomname(1) && atomname(1) <= '9'
      atype(n) = atomname(2);
    else
      atype(n) = atomname(1);
    end
%    if strcmp(resname(1:3),'TIP') || strcmp(resname(1:3),'HOH')
%      atype(n) = lower(atype(n));
%    end
    x = str2num(line(31:38));
    y = str2num(line(39:46));
    z = str2num(line(47:54));
    crg(n)   = str2num(line(61:66));
    crd(n,:) = [x y z];
  end
end

crd = crd(1:n,:);  
crg = crg(1:n);
atype = atype(1:n,:);
fclose(fid);
