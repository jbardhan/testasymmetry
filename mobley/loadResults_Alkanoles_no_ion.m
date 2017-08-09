close all; clear all
addpath('export_fig/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Set your toggles and define the solute and solvent lists %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploton = 1;

solvents = {'water','methanol','ethanol','propanol','butanol','octanol'}; 
solvents_1 = {'butanol'};
solvents_2 = {'octanol'};       
%common_solutes = {'ethanol','butanone','n_octane','pyrene','cyclohexane',};
        
testset  = {'4_nitroaniline','anthracene','benzene','butan_1_ol','butanone','cyclohexane','ethanol','methanol','nitromethane','phenanthrene','propan_1_ol','toluene'};
%testset_ions = {'Li','Na','K','Cl','Br','I'};
%but_ions= {'Rb','Na','K','Cl','Br','I'};
%test_ion = {'Rb','Cs'};


% This routine will plot all of the calulated vs. experimental results if 
% plot is set to 1.  Futhermore, it creates structures that contain all of
% the results for calculated solvation free energies as a function of each 
% solvent and again as a function of solute
for i = 1:length(solvents)
          
        fid = fopen(['mnsol/',lower(solvents{i}),'.csv'],'r'); 
        Data = textscan(fid,'%s %f %f','delimiter',',');
        fclose(fid);
        mol_list = Data{1};
        allCases = 1:length(mol_list);
        [~,m] = ismember(testset,mol_list);
        allCases(m) = 0;
        %if ~ismember(solvents{i},solvents_1)
        %    [~,n] = ismember(testset_ions,mol_list);
        %    allCases(n) = 0;
        %else
        %    [~,n] = ismember(but_ions,mol_list);
        %   allCases(n) = 0;
        %end
        %if ~ismember(solvents{i},solvents_2)
        %    [~,p] = ismember(test_ion,mol_list);
        %    allCases(p) = 0;
        %end
        temp = load(['Run',solvents{i}]);
%     solute_errors(i,:) = temp.errfinal(n);
%    sorted_errors = sort(abs(temp.errfinal));
        results(i) = struct('Solvent', solvents{i},'Num_Solutes',round(length(mol_list),3),'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE),...
                        'Mean_Abs_error',mean(abs(temp.errfinal)));
        if ploton
            figure
            plot(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)),'bo','markers',12)
           %plot(temp.refE,temp.calcE,'bo','markers',12)
            set(gca,'FontSize',15)
            hold on
            plot(temp.refE(m),temp.calcE(m),'rs','markers',12)
            minax = round(min(min(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)))));
            maxax = round(max(max(temp.refE(allCases(allCases~=0)),temp.calcE(allCases(allCases~=0)))));
            axis([minax-2 maxax+2 minax-2 maxax+2]);
            foo = refline(1,0);
            set(foo,'Linewidth',2,'color','k');
            xlabel(['\Delta','G_{expt}^{solv, ',solvents{i},'}',' (kcal.mol^{-1})'])
            ylabel(['\Delta','G_{calc}^{solv, ',solvents{i},'}',' (kcal.mol^{-1})'])
            legend('Predictions','Training Set','Location','southeast')
            filename = sprintf('Output/Figures/DeltaG-%s.PDF',solvents{i});
            export_fig(filename,'-painters','-transparent');

        end
end
%{
clear 
addpath('export_fig/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Set your toggles and define the solute and solvent lists %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ploton = 1;
solvents = {'water','methanol','ethanol','propanol','butanol','octanol'}; 
solvents_1 = {'butanol'};
solvents_2 = {'octanol'};       
%common_solutes = {'ethanol','butanone','n_octane','pyrene','cyclohexane',};
        
testset_ions = {'Li','Na','K','Cl','Br','I'};
but_ions= {'Rb','Na','K','Cl','Br','I'};
test_ion = {'Rb','Cs'};
test_Cs={'Cs'};

for i = 1:length(solvents)
           
    fid = fopen(['mnsol/',lower(solvents{i}),'_ions.csv'],'r'); 
    Data = textscan(fid,'%s %f %f','delimiter',',');
    fclose(fid);
    mol_list = Data{1};
    
    
    allCases = 1:length(mol_list);
    if ~ismember(solvents{i},solvents_1)
        [~,m] = ismember(testset_ions,mol_list);
        if ~ismember(solvents{i},solvents_2)
            [~,n] = ismember(test_ion,mol_list);
        end
    else
        [~,m] = ismember(but_ions,mol_list);
        [~,n] = ismember(test_Cs,mol_list);
    end 
    
    allCases(~m) = 0;
    allCases(~n) = 0;
    
    temp = load(['Run',solvents{i}]);
%     solute_errors(i,:) = temp.errfinal(n);
%    sorted_errors = sort(abs(temp.errfinal));
    results(i) = struct('Solvent', solvents{i},'Num_Solutes',round(length(mol_list),3),'CalcE',temp.calcE,...
                        'errfinal',temp.errfinal,'es',temp.es,'np',...
                        temp.np,'refE',temp.refE,'RMS',rms(temp.calcE-temp.refE),...
                        'Mean_Abs_error',mean(abs(temp.errfinal)));
    if ploton
        figure
        plot(temp.refE(allCases(n)),temp.calcE(allCases(n)),'bo','markers',20,'Linewidth',3)
        hold on
        plot(temp.refE(allCases(m)),temp.calcE(allCases(m)),'rs','markers',20,'Linewidth',3)
        set(gca,'FontSize',15)
        hold on
        %plot(temp.refE,temp.calcE,'bo','markers',12)
        
        minax = round(min(min(temp.refE(allCases(m)),temp.calcE(allCases(m)))));
        maxax = round(max(max(temp.refE(allCases(m)),temp.calcE(allCases(m)))));
        axis([minax-2 maxax+2 minax-2 maxax+2]);
        foo = refline(1,0);
        set(foo,'Linewidth',2,'color','k');
        xlabel(['\Delta','G_{expt}^{solv, ',solvents{i},'}',' (kcal.mol^{-1})'])
        ylabel(['\Delta','G_{calc}^{solv, ',solvents{i},'}',' (kcal.mol^{-1})'])
        legend('Predictions','Training Set','Location','southeast')
        filename = sprintf('Output/Figures/DeltaG-%s.PDF',solvents{i});
        export_fig(filename,'-painters','-transparent');

    end

end


%Outliers = readErr(max_err);
%save('Solvent Outliers','Outliers');
%}