curdir = pwd;
Radii = [5.0,7.5,10.0,15.0,20.0,30.0,40.0];
leg = {'R = 5.0 A','R = 7.5 A','R = 10.0 A','R = 15.0 A',...
       'R = 20.0 A','R = 30.0 A','R = 40.0 A'};
x = [25.0,35.0,40.0,42.5,45.0,46.0,47.0,47.5,48.0,48.5];
style={'-o','-s','-d','-x','-*','-<','-^'};
colRed = colormap(lines(8));
colGreen = colormap(winter(8));
%face = {colGreen(6,:),colRed(4,:),'b','c','m','y','k'};
i = 1;
for R=Radii
SLfileName = sprintf('%s/convergence_neg/%3.1f/%s%3.1f%s',curdir,R,'RunMembrane_SLIC_',R,'.mat');
SLICdata = load(SLfileName);
hL(i) = plot(x,SLICdata.es(2:11),style{i},...
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
title_text = sprintf('Electrostatic energy at different charge locations (-1 charge) as the  \nradius of the membrane varies from %d to %d Angstroms',5,40);
title(title_text)
legend(leg);


