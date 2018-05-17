clear all
close all
clc

Home = getenv('HOME');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%                                Set these values before running the code                      %%%%%%%%%
%%%%%%                                                                                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

repo_path=sprintf('%s/repos',Home);
dropbox_path=sprintf('%s/Dropbox',Home);
                   
dataset='mobley';   % options are mobley or mnsol , in the mnsol case we use mobley syrface areas               
                    
calcflag=1;     % if calcflag=1 the code actually calculate the /delta G 's using BEM if it is zero, it means 
                % delta G's has been calculated before and all we need is
                % to load the data. 
                    
temp_min=4.85;     % lower bound of the temperature interval 
temp_max=44.85;    % upper bound in the temperature interval
tempdiv=5;      % number of divisions in the temperature interval                                                   11      
                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));


if calcflag==1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load data from parameterization


    %%% define the new tempereture vector 
    new_temp=linspace(temp_min,temp_max,tempdiv);
    
    ionflag=0;      % if ionflag=0, ions data are not included in the testset 
                    % if ionflag=1, ions data are included in the testset


    %%% Finding the values of parameters from the quadratic fit at the new
    %%% temperatures
    
    if ionflag==1
        ParamWatInfo=load('OptWater_thermo');
    elseif ionflag==0
        ParamWatInfo=load('OptWater_thermo_wo_ion');
    end
        
    x = ParamWatInfo.xvec;
    %x=x(2:4,:);
    
    for j=1:length(new_temp)
        clear global
        loadConstants

        convertKJtoKcal = 1/joulesPerCalorie;
        global UsefulConstants ProblemSet saveMemory writeLogfile logfileName
        saveMemory = 1;
        writeLogfile = 1;
        logfileName = 'junklogfile';
        
        epsIn  =  1;
        Tbase = 300; 
        %epsOut = 78.36; % from MNSol
        
        mytemp=Tbase;
        KelvinOffset = 273.15;
        conv_factor = 332.112;
        staticpotential = 0.0; % this only affects charged molecules;
        kappa = 0.0;  % should be zero, meaning non-ionic solutions!
        
        
        temp=new_temp(j);
        epsOut = (-1.410e-6)*temp^3+(9.398e-4)*temp^2-0.40008*temp+87.740;


        % the staticpotential below should not be used any more, please check
        UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
                     kappa,'conv_factor',conv_factor,...
                     'staticpotential',staticpotential);

        
        if strcmp(dataset,'mnsol')
                 
                 
            fid = fopen('mnsol/water.csv','r'); 
            Data = textscan(fid,'%s %f %f','delimiter',',');
            fclose(fid);
            mol_list = Data{1};
            dG_list = Data{2};

            fid = fopen('mnsol/mobley_dG_AND_sa_and_vol_fixed.csv','r');
            Data = textscan(fid,'%s %f  %f  %f  %f  %f  %f  %f','delimiter',',');
	    fprintf('Dont use this option for now!\n');
            fclose(fid);
            all_solutes = Data{1};
            all_surfAreas = Data{3};
            [m, index] = ismember(mol_list,all_solutes);
            surfArea_list = all_surfAreas(index);
            
        elseif strcmp(dataset,'mobley')
            
            fid = fopen('mnsol/mobley_dG_AND_sa_and_vol_fixed.csv','r');
            Data = textscan(fid,'%s %f  %f  %f  %f  %f  %f  %f','delimiter',',');
            fclose(fid);
            
            if ionflag==1
                
                mol_list = Data{1};
                dG_list = Data{2};
                surfArea_list = Data{3};
                calc_mobley= Data{8};
            
            elseif ionflag==0
                
                mol_list = Data{1}(1:502);
                dG_list = Data{2}(1:502);
                surfArea_list = Data{3}(1:502);
                calc_mobley= Data{8}(1:502);
                
            end
        end
        

                 
        curdir = pwd;
        for i=1:length(mol_list)
          dir=sprintf('%s/lab/projects/slic-jctc-mnsol/nlbc-mobley/nlbc_test/%s',dropbox_path,mol_list{i});
          chdir(dir);
          pqrData = loadPqr('test.pqr');
          pqrAll{i} = pqrData;
          srfFile{i} = sprintf('%s/test_2.srf',dir);
          chargeDist{i} = pqrData.q;
          referenceData{i} = dG_list(i);
          surfArea{i} = surfArea_list(i);
          chdir(curdir);
          addProblemSA(mol_list{i},pqrAll{i},srfFile{i},chargeDist{i},referenceData{i},surfArea{i});
        end         

        [errfinal(j,:),calcE(j,:),refE(j,:),es(j,:),np(j, :)]=ObjectiveFromBEMSA(x(j,:));
    end
      
    if strcmp(dataset,'mnsol')
        save('RunWater_MNnsol_mobleysurf_184_ini_run','errfinal','calcE','refE','es','np','new_temp','x','mol_list');
    elseif strcmp(dataset,'mobley')
        save('RunWater_mobley_513_ini_run','errfinal','calcE','refE','es','np','new_temp','x','mol_list','calc_mobley');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(dataset,'mnsol')
    RunWater=load('RunWater_MNnsol_mobleysurf_184_ini_run.mat');
elseif strcmp(dataset,'mobley')
    RunWater=load('RunWater_mobley_513_ini_run');
end

dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures


loadConstants
convertKJtoKcal = 1/joulesPerCalorie;
calc_mobley=0;

errfinal=RunWater.errfinal;
calcE = RunWater.calcE; %calculated values for \Delta G at different temperatures
refE = RunWater.refE; %calculated values for \Delta G at different temperatures
%refS = RunWater_mobleyset_184.dS_list; %calculated values for \Delta G at
%different temperatures
if strcmp(dataset,'mobley')
%     calcE = RunWater.calcE(:,1:502); %calculated values for \Delta G at different temperatures
%     refE = RunWater.refE(:,1:502); %calculated values for \Delta G at different temperatures
%     calcE_ion = RunWater.calcE(:,503:end); %calculated values for \Delta G at different temperatures
%     refE_ion = RunWater.refE(:,503:end); %calculated values for \Delta G at different temperatures
     calcE = RunWater.calcE; %calculated values for \Delta G at different temperatures
     refE = RunWater.refE; %calculated values for \Delta G at different temperatures
     calc_mobley=RunWater.calc_mobley(1:502);
end
    
es=RunWater.es;
np=RunWater.np;
x = RunWater.x;   % parameters at different temperetures
TEMP=RunWater.new_temp';        % create the temperature vector '
TEMP_K=TEMP+273.15;
[m,index]=ismember(24.85,TEMP);
mol_list=RunWater.mol_list;
if strcmp(dataset,'mobley')
%     mol_list=RunWater.mol_list(1:502);
    mol_list=RunWater.mol_list;
end

for i=1:length(mol_list)
    f = @(R) (R(1)-R(2)*(TEMP_K-298)+R(3)*((TEMP_K-298)-TEMP_K.*log(TEMP_K./298)))-calcE(:,i);
    R0=[refE(index,i),1,1];
    options=optimoptions('lsqnonlin','StepTolerance',1e-6);
    options=optimoptions(options,'OptimalityTolerance',1e-6);
    options=optimoptions(options,'FunctionTolerance',1e-6);
    [R,resnorm,residual,exitflag,output]=lsqnonlin(f,R0,[],[],options);
%      figure()
%      plot(TEMP_K,calcE(:,i),'ko',TEMP_K,R(1)-R(2)*(TEMP_K-298)+R(3)*((TEMP_K-298)-TEMP_K.*log(TEMP_K./298)),'b-')
    dGfunc(i).name=mol_list(i); 
    dGfunc(i).dg=R(1);
    dGfunc(i).ds=R(2);
    dGfunc(i).cp=R(3); 
    dsvec(i)=dGfunc(i).ds*1000;
    cpvec(i)=dGfunc(i).cp*1000;
    resnorm(i)=resnorm;
    exitflag(i)=exitflag;
    output(i)=output;
end


 dg_rms_298_MD=0;
 dg_rms_298=rms(calcE-refE);
 
 if strcmp(dataset,'mobley')
    dg_rms_298_mol=rms(calcE(index,1:502)-refE(index,1:502));
    dg_rms_298_MD_mol=rms(calc_mobley(1:502)'-refE(index,1:502));
    if ionflag == 1
        dg_rms_298_ion=rms(calcE(index,503:end)-refE(index,503:end));
 
    end
 end


if strcmp(dataset,'mnsol')
    output_name='RunWater_mnsol_mobleySurf_thermo';
elseif strcmp(dataset,'mobley')
    if ionflag == 1
        output_name='RunWater_mobley_thermo';
    elseif ionflag == 0
        output_name='RunWater_mobley_thermo_wo_ion';
    end
end
if ionflag == 1
    save(output_name,'errfinal','calcE','refE','es','np','TEMP','x','mol_list','dGfunc','dsvec','cpvec','dg_rms_298_ion','dg_rms_298_mol','dg_rms_298_MD_mol','index','calc_mobley','resnorm','residual','output','exitflag');
elseif ionflag == 0
    save(output_name,'errfinal','calcE','refE','es','np','TEMP','x','mol_list','dGfunc','dsvec','cpvec','dg_rms_298_mol','dg_rms_298_MD_mol','index','calc_mobley','resnorm','residual','output','exitflag');
end




    