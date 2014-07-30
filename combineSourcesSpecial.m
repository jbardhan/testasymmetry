function [pqrNew, RNew] = combineSourcesSpecial(pqr1, R1, pqr2, R2, ...
						ind)
pqrNew = struct('xyz', [pqr1.xyz; pqr2.xyz(ind,:)],...
		'q', [pqr1.q; 1],...
		'R', [pqr1.R; pqr2.R(ind)]);
RNew = [R1; R2(ind)];