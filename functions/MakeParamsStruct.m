function Params = MakeParamsStruct(x)
Params = struct('alpha',x(1),'beta',x(2),'EfieldOffset',x(3),'mu',x(4),'phiStatic',x(5));
fprintf('Params: alpha = %f\tbeta = %f\tgamma = %f\n',x(1),x(2),x(3));