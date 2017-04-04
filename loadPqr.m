function pqr = loadPqr(filename)
fh = fopen(filename,'r');
C = textscan(fh, '%s%d%s%s%d%f%f%f%f%f');
fclose(fh);
xyz = [C{6} C{7} C{8}];
q = C{9};
r = C{10};
resnum = C{5};
resid = C{4};
atomid = C{3};
atomnum = C{2};
pqr = struct('xyz', xyz, 'q', q, 'r', r, 'resnum',resnum,'resid',0, ...
	     'atomid', 0, 'atomnum',atomnum);
pqr.resid = resid;
pqr.atomid = atomid;
