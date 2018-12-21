
curdir = pwd;
i=20;
fileName = sprintf('%s/convergence_neg/%3.1f/%s%3.1f%s',curdir,i,'RunMembrane_SLIC_',i,'.mat');
S = load(fileName);
fileNamePB = sprintf('%s/convergence_neg/%3.1f/%s%3.1f%s',curdir,i,'RunMembrane_PB_',i,'.mat');
SPB = load(fileNamePB);
B = {S.es.',SPB.es.'};   
x=[25.0,30.0,35.0,40.0,42.5,45.0,46.0,47.0,47.5,48.0,48.5];
colRed = colormap(hot(8));
colGreen = colormap(winter(8));

hL =zeros(2,1);
face = {colGreen(6,:),colRed(4,:)};
style={'-o','-s'};
for i=1:2
hL(i) = plot(50-x,B{i},style{i},...
        'Color',face{i},...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',face{i},...
        'MarkerSize',10);
hold on
end

for i=1:2
hL(i) = plot(x,B{i},style{i},...
        'Color',face{i},...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',face{i},...
        'MarkerSize',10);
hold on
end

xlabel('Membrane depth (A)')
ylabel('Electrostatic energy (kcal/mol)') 
hLeg = legend(hL,'SLIC','PB');

%{
Radii = [5.0,7.5,10.0,15.0,20.0,30.0,40.0];
x = x(6:11);
style={'-o','-s','-d','-x','-*','+','^'};
face = {colGreen(6,:),colRed(4,:),'b','c','m','y','k'};
i = 1;
for R=Radii
    SLfileName = sprintf('%s/convergence_neg/%3.1f/%s%3.1f%s',curdir,R,'RunMembrane_SLIC_',R,'.mat');
    SLICdata = load(SLfileName);
    hL(i) = plot(x,SLICdata.es(6:11),style{i},...
        'Color',face{i},...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',face{i},...
        'MarkerSize',10);
    i =+ 1;
hold on
end
xlabel('Membrane depth (A)')
ylabel('Electrostatic energy (kcal/mol)') 
%hLeg = legend(hL,'SLIC','PB');
%}

