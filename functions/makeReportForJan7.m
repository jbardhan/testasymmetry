%Optimize_Params_at_temp

delS = [-199 -143 -101 -92 -78 -138 -89 -79 -66].*(0.239/10^3)'; % Delta S at 25 C from Fawcett Ch. 3
delG = [-483 -403 -333 -316 -288 -396 -310 -291 -264].*0.239'; % Delta G at 25 C from Fawcett Ch. 3


for l = 1:length(t)
del_t(l) = (t(l)-25);

end



for i = 1:length(t)
MD(i,:) = [delG(1)-del_t(i)*delS(1)    % Li+
       delG(2)-del_t(i)*delS(2)    % Na+
       delG(3)-del_t(i)*delS(3)    % K+
       delG(4)-del_t(i)*delS(4)    % Rb+
       delG(5)-del_t(i)*delS(5)    % Cs+
       delG(6)-del_t(i)*delS(6)    % F-
       delG(7)-del_t(i)*delS(7)    % Cl-
       delG(8)-del_t(i)*delS(8)    % Br-
       delG(9)-del_t(i)*delS(9)    % I-
       ]'; % The constant is Delta G at 25 C from Fawcett Ch. 3
Params = struct('alpha',x(i,1), 'beta', x(i,2),'EfieldOffset',x(i, ...
						  3));
F(i,:) = Level_2_MSA(Params,t(i)); % Calls the function Level_2_MSA.m which calculates the energy using bornPicardNoStern.
end 

for i=1:length(delG)
  junk = polyfit(t', F(:,i),1);
  delS_NLBC(i,1) = junk(1);
  junk = polyfit(t', MD(:,i),1);
  delS_MD(i,1) = junk(1);
end

