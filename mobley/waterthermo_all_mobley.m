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
dropbox_path=sprintf('%s/Dropbox',Home);

ionflag=1;          % if ionflag=0, ions data are not included in the testset 
                    % if ionflag=1, ions data are included in the testset
paramboundflag=1;   % if paramboundflag=0 there is no bound for parameters in the optimization process 
                    % if paramboundflag=1 there is abound
                    
calcflag=0;     % if calcflag=1 the code actually calculate the /delta G 's using BEM if it is zero, it means 
                % delta G's has been calculated before and all we need is
                % to load the data. 
                    
temp_min=288;     % lower bound of the temperature interval 
temp_max=308;    % upper bound in the temperature interval
tempdiv=3;      % number of divisions in the temperature interval                     
                    
                    
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
    

    %%% Finding the values of parameters from the quadratic fit at the new
    %%% temperatures
    ParamWatInfo=load('OptWater_w_ion_wo_florine');
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
        
        new_temp=new_temp-KelvinOffset;
        temp=new_temp(j);
        epsOut = (-1.410e-6)*temp^3+(9.398e-4)*temp^2-0.40008*temp+87.740;


        % the staticpotential below should not be used any more, please check
        UsefulConstants = struct('epsIn',epsIn,'epsOut',epsOut,'kappa', ...
                     kappa,'conv_factor',conv_factor,...
                     'staticpotential',staticpotential);

        
                 
        fid = fopen('mnsol/water.csv','r'); 
        Data = textscan(fid,'%s %f %f','delimiter',',');
        fclose(fid);
        mol_list = Data{1};
        dG_list = Data{2};

        fid = fopen('mnsol/mobley_sa.csv','r');
        Data = textscan(fid,'%s %f','delimiter',',');
        fclose(fid);
        all_solutes = Data{1};
        all_surfAreas = Data{2};
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

    save('RunWater_mobleyset_184','errfinal','calcE','refE','es','np','new_temp','x','mol_list');
    
end


RunWater_mobleyset_184=load('RunWater_mobleyset_184.mat');

dGfunc=struct();  % structure that has the information of the solutes and the linear function that fits the calculated values of dG at different temperatures


loadConstants
convertKJtoKcal = 1/joulesPerCalorie;

x = RunWater_mobleyset_184.x;   % parameters at different temperetures
calcE = RunWater_mobleyset_184.calcE; %calculated values for \Delta G at different temperatures
refE = RunWater_mobleyset_184.refE; %calculated values for \Delta G at different temperatures
%refS = RunWater_mobleyset_184.dS_list; %calculated values for \Delta G at different temperatures
TEMP=RunWater_mobleyset_184.new_temp;        % create the temperature vector 
mol_list=RunWater_mobleyset_184.mol_list;

for i=1:length(mol_list)
    dGfunc(i).name=mol_list(i);
    dGfunc(i).func=fit(TEMP',calcE(:,i),'poly1'); % the linear function that fits to the calculated dG  at differet temperatures
end

for i=1:length(mol_list)
    p=[dGfunc(i).func.p1,dGfunc(i).func.p2];
    pder=polyder(p);  % derivative of the linear function dG
    dsvec(i)=-polyval(pder,24.85)*1000;    % Evaluationg entropy, dS at 298K = 24.85C in cal/mol/K
%     figure(i);
%     plot(dGfunc(i).func,TEMP',calcvec(:,i),'o')   
end


% Reference_dG =refE
% Calculated_dG=calcE
% 
% dG_err=abs(Reference_dG-Calculated_dG)
% 
% Reference_dS =refS'*1000
% Calculated_ds=dsvec
% 
% dS_err=abs(Reference_dS-Calculated_ds)

figure()
plot(refE(2,:),calcE(2,:),'o')
hold on
line([min(refE(2,:)), max(refE(2,:))],[min(refE(2,:)) ,max(refE(2,:))])
rms(calcE(2,:)-refE(2,:))

mol_list_num=[1:length(mol_list)];
scatter(mol_list_num,dsvec,'rx')

    