close all
for i=[1:10,11,13,14,16:20,22,24,25,26]
    fig = figure('Renderer', 'painters', 'Position', [8 6 800 600]);
    mat_file = sprintf('OptCosmoBondiiNp_%d.mat',i);
    load(mat_file);
    plot(cav_mob,cav+comb,'r*','markers',10)
    hold on
    plot(disp_mob,disp,'b*','markers',10)
    hold on
    plot(np_mob,np,'g*','markers',10)
    titl = sprintf("i = %d - q = %f - z = %f - cav_{factor} = %f - rmse = %f",i,x(18),x(19),x(20),rmse);
    title(titl)
    %hold off
end
