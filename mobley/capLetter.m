function str = capLetter(strToCap)

for i = 1:length(strToCap)
   if i == 1
      str(i) = upper(strToCap(i));
   else
      str(i) = strToCap(i);
   end
end