function [RMSdGTrans, MeanAbsdGTrans, dGTransCalc,...
    dGTransRef] = calcRMSTrans(fromSolvent, toSolvent)

fidFrom = fopen(strcat('mnsol/',lower(fromSolvent),'.csv'),'r'); 
DataFrom = textscan(fidFrom,'%s %f %f','delimiter',',');
fclose(fidFrom);
mol_listFrom = DataFrom{1};

fidTo = fopen(strcat('mnsol/',lower(toSolvent),'.csv'),'r'); 
DataTo = textscan(fidTo,'%s %f %f','delimiter',',');
fclose(fidTo);
mol_listTo = DataTo{1};


solvents = {capLetter(fromSolvent), capLetter(toSolvent)};
solvents = strcat('Run',solvents,'.mat');

filesInDir = dir(pwd);

there1 = 0;
there2 = 0;
for i = 1:length(filesInDir)
    if filesInDir(i).isdir == 0
       if strcmp(filesInDir(i).name,solvents{1})
           fromSolventDat = load(filesInDir(i).name);
           there1 = 1;
       end
       
       if strcmp(filesInDir(i).name,solvents{2})
           toSolventDat = load(filesInDir(i).name);
           there2 = 1;
       end
       
       if (there1 && there2)
           break
       end
       
    end
end 



fromdGCalc = fromSolventDat.calcE;
fromdGRef = fromSolventDat.refE;
todGCalc = toSolventDat.calcE;
todGRef = toSolventDat.refE;

fromLen = length(fromdGCalc);
toLen = length(todGCalc);

if fromLen > toLen
   [m1,n1] = ismember(mol_listFrom,mol_listTo);
   [m2,n2] = ismember(mol_listTo,mol_listFrom);
   if length(n2(n2~=0))~=length(n1(n1~=0))
       fprintf('Error in loading data! \n');
       keyboard
   end
   dGTransCalc = fromdGCalc((n2(n2~=0)),:) - todGCalc((n1(n1~=0)),:);
   dGTransRef = fromdGRef((n2(n2~=0)),:) - todGRef((n1(n1~=0)),:);
elseif toLen > fromLen
   [m1,n1] = ismember(mol_listTo,mol_listFrom);
   [m2,n2] = ismember(mol_listFrom,mol_listTo);
   if length(n2(n2~=0))~=length(n1(n1~=0))
       fprintf('Error in loading data! \n');
       keyboard
   end
   dGTransCalc = fromdGCalc((n1(n1~=0)),:) - todGCalc((n2(n2~=0)),:);
   dGTransRef = fromdGRef((n1(n1~=0)),:) - todGRef((n2(n2~=0)),:);
elseif toLen == fromLen
   dGTransCalc = fromdGCalc - todGCalc;
   dGTransRef = fromdGRef - todGRef;
end

RMSdGTrans = rms(dGTransCalc-dGTransRef);
MeanAbsdGTrans = mean(abs(dGTransCalc-dGTransRef));
foo = 9;


