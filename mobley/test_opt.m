 clear all
 clc
 dG_ref=8.1; %Hess in kcal/mol
 H_ref=-8.3;  %Hess . in kcal/mol
 CP_ref=1e-3*142;
 dS_ref=(H_ref-dG_ref)/298;  % in kcal/mol/K
 
 TEMP_K=linspace(14.85,34.85,3)+273.15;
 
 dG=dG_ref-dS_ref*(TEMP_K-298)+CP_ref*((TEMP_K-298)-TEMP_K.*log(TEMP_K./298));
 
 r = (0.1).*rand(3,1) + 0.95;
 
 dG=r'.*dG;
 
 f = @(R) (R(1)-R(2)*(TEMP_K-298)+R(3)*((TEMP_K-298)-TEMP_K.*log(TEMP_K./298)))-dG;
 R0=[0,0,0];
 options=optimoptions('lsqnonlin','StepTolerance',1e-12);
 options=optimoptions(options,'OptimalityTolerance',1e-12);
 options=optimoptions(options,'FunctionTolerance',1e-12);
 [R,resnorm,residual,exitflag,output]=lsqnonlin(f,R0,[],[],options);