% Finding Delta S 
load('Params_NLBC_295K');
load('Params_NLBC_300K');
load('Params_NLBC_305K');

T = [295 300 305];

ParamsAt295K = MakeParamsStruct(x1);
ParamsAt300K = MakeParamsStruct(x2);
ParamsAt305K = MakeParamsStruct(x3);

calculatedE295K = CalculateEnergiesFromBEM(ParamsAt295K);
calculatedE300K = CalculateEnergiesFromBEM(ParamsAt300K);
calculatedE305K = CalculateEnergiesFromBEM(ParamsAt305K);

