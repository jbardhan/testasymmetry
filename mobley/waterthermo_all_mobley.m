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

repo_path=sprintf('%s/Research',Home);
dropbox_path=sprintf('%s/Dropbox-NEU/Dropbox',Home);
                   
dataset='mobley';   % options are mobley or mnsol , in the mnsol case we use mobley syrface areas               
                    
calcflag=0;     % if calcflag=1 the code actually calculate the /delta G 's using BEM if it is zero, it means 
                % delta G's has been calculated before and all we need is
                % to load the data. 
                    
temp_min=14.85;     % lower bound of the temperature interval 
temp_max=34.85;    % upper bound in the temperature interval
tempdiv=3;      % number of divisions in the temperature interval                     
                    
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

addpath(sprintf('%s/pointbem',repo_path));
addpath(sprintf('%s/panelbem',repo_path));
addpath(sprintf('%s/testasymmetry',repo_path));
addpath(sprintf('%s/testasymmetry/functions',repo_path));
addpath(sprintf('%s/testasymmetry/mobley',repo_path));
addpath(sprintf('%s/testasymmetry/born',repo_path));


for mm=1:5
    
    if calcflag==1
    
   


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Load data from parameterization


        %%% define the new tempereture vector 
        new_temp=linspace(temp_min,temp_max,tempdiv);


        %%% Finding the values of parameters from the quadratic fit at the new
        %%% temperatures

        mat_name=sprintf('OptWater_thermo_rand_%d',mm);

        ParamWatInfo=load(mat_name);

        x = ParamWatInfo.xvec;
        x=x(2:4,:);

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

            mol_list=ParamWatInfo.testset_total;
            dG_list =ParamWatInfo.dG_list_total;
            fid = fopen('mnsol/mobley_dG_AND_sa_and_vol.csv','r');
            Data = textscan(fid,'%s %f  %f  %f  %f  %f  %f  %f','delimiter',',');
            fclose(fid);
            all_solutes = Data{1};
            all_surfAreas = Data{3};
            [m, index] = ismember(mol_list,all_solutes);
            surfArea_list = all_surfAreas(index);

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

        run_out__ini_name=sprintf('RunWater_total_ini_run_%d',mm);

        save(run_out__ini_name,'errfinal','calcE','refE','es','np','new_temp','x','mol_list');



    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    run_in_name=sprintf('RunWater_total_ini_run_%d',mm);
    
    RunWater=load(run_in_name);


    dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures


    loadConstants
    convertKJtoKcal = 1/joulesPerCalorie;
    calc_mobley=0;

    errfinal=RunWater.errfinal;
    calcE = RunWater.calcE; %calculated values for \Delta G at different temperatures
    refE = RunWater.refE; %calculated values for \Delta G at different temperatures
    %refS = RunWater_mobleyset_184.dS_list; %calculated values for \Delta G at
    %different temperatures

    es=RunWater.es;
    np=RunWater.np;
    x = RunWater.x;   % parameters at different temperetures
    TEMP=RunWater.new_temp';        % create the temperature vector '
    TEMP_K=TEMP+273.15;
    [m,index]=ismember(24.85,TEMP);
    mol_list=RunWater.mol_list;


    for i=1:length(mol_list)
        f = @(R) (R(1)-R(2)*(TEMP_K-298)+R(3)*((TEMP_K-298)-TEMP_K.*log(TEMP_K./298)))-calcE(:,i);
        R0=[0,0,0];
        options=optimoptions('lsqnonlin','StepTolerance',1e-12);
        options=optimoptions(options,'OptimalityTolerance',1e-12);
        options=optimoptions(options,'FunctionTolerance',1e-12);
        [R,resnorm,residual,exitflag,output]=lsqnonlin(f,R0,[],[],options);
%           figure()
%           plot(TEMP_K,calcE(:,i),'ko',TEMP_K,R(1)-R(2)*(TEMP_K-298)+R(3)*((TEMP_K-298)-TEMP_K.*log(TEMP_K./298)),'b-')
        dGfunc(i).name=mol_list(i); 
        dGfunc(i).dg=R(1);
        dGfunc(i).ds=R(2);
        dGfunc(i).cp=R(3); 
        dsvec(i)=dGfunc(i).ds*1000;
        cpvec(i)=dGfunc(i).cp*1000;
    end


   dg_rms_298=rms(calcE-refE);

    output_name=sprintf('RunWater_total_thermo_%d',mm);


    save(output_name,'errfinal','calcE','refE','es','np','TEMP','x','mol_list','dGfunc','dsvec','cpvec','dg_rms_298','index');
    
end

    