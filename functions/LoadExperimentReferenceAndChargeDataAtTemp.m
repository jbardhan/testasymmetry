delS = [-133; -96; -87; -81; -53].*(0.239/10^3);
delG = [-424; -352; -329; -306; -304].*0.239;
referencetemp = 25+KelvinOffset; % Fawcett data is at 25 C
delt = mytemp - referencetemp; % Kelvin


delG_new = delG-delt*delS;

NaReference = delG_new(1);
KReference = delG_new(2);
RbReference = delG_new(3);
CsReference = delG_new(4);
ClReference = delG_new(5);

CationChargePlusOne = 1;
AnionChargeMinusOne = -1;