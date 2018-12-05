curdir = pwd;
Radii = [5.0,7.5,10.0,15.0,20.0,30.0,40.0];
leg = {'5.0','7.5','10.0','15.0','20.0','30.0','40.0'};
x = [35.0,40.0,42.5,45.0,46.0,47.0,47.5,48.0,48.5];
style={'-o','-s','-d','-x','-*','-<','-^'};
colRed = colormap(lines(8));
colGreen = colormap(winter(8));
%face = {colGreen(6,:),colRed(4,:),'b','c','m','y','k'};
i = 1;
for R=Radii
SLfileName = sprintf('%s/convergence/%3.1f/%s%3.1f%s',curdir,R,'RunMembrane_',R,'.mat');
SLICdata = load(SLfileName);
hL(i) = plot(x,SLICdata.es(3:11),style{i},...
        'Color',colRed(i,:),...
        'LineWidth',2,...
        'MarkerEdgeColor',colRed(i,:),...
        'MarkerFaceColor',colRed(i,:),...
        'MarkerSize',10);
i = i + 1;   
hold on
end 
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',16)

xlabel('Membrane depth (A)', 'FontSize', 14)
ylabel('Electrostatic energy (kcal/mol)', 'FontSize', 14)
legend(leg);



