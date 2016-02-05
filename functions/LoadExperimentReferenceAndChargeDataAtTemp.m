delS = [-133 -96 -87 -81 -53].*(0.239/10^3)';
delG = [-424 -352 -329 -306 -304].*0.239';
t = mytemp-273.15; 
delt = t-25; % Celcius


delG_new = [delG(1)-delt*delS(1)
            delG(2)-delt*delS(2)
            delG(3)-delt*delS(3)
            delG(4)-delt*delS(4)
            delG(5)-delt*delS(5)];

NaReference = delG_new(1);
KReference = delG_new(2);
RbReference = delG_new(3);
CsReference = delG_new(4);
ClReference = delG_new(5);

CationChargePlusOne = 1;
AnionChargeMinusOne = -1;