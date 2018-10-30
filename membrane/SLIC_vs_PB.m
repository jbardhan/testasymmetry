
curdir = pwd;
i=10;
fileName = sprintf('%s/%d/%s%d%s',curdir,i,'RunMembrane_',i,'.mat');
S = load(fileName);
fileNamePB = sprintf('%s/%d/%s%d%s',curdir,i,'RunMembranePB_',i,'.mat');
SPB = load(fileNamePB);
B = {S.es.',SPB.es.'};   
x=(25:-2.5:0);
x(11) = 0.1;
colRed = colormap(hot(8));
colGreen = colormap(summer);

hL =zeros(2,1);
face = {colGreen(6,:),colRed(4,:)};
style={'-o','-s'};
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
