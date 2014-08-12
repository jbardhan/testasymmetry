leftstart = 0.06;
topstart = 0.54;
rightstart = 0.55;
botstart = 0.09;
imwidth = 0.30;
imheight = 0.45;

ionSurface_nl = reshape(phinl_part0(2:end,2), size(Ygrid,1), ...
								size(Ygrid,2));
figure;
subplot(2,2,1);
H_ionsurface_nl = contourf(Ygrid, Zgrid, ionSurface_nl,'linestyle','none');
axis equal;
colorbar;
set(gca,'fontsize',16,'position',[leftstart topstart imwidth imheight]);
title('(a)');
%xlabel('Y');
%ylabel('Z');

ionSurface_l = reshape(phil_part0(2:end,2), size(Ygrid,1), ...
								size(Ygrid,2));
subplot(2,2,2);
H_ionsurface_l = contourf(Ygrid, Zgrid, ionSurface_l,'linestyle','none');
axis equal;
colorbar;
set(gca,'fontsize',16,'position',[rightstart topstart imwidth imheight]);
title('(b)');
%xlabel('Y');
%ylabel('Z');

ionBuried_nl = reshape(phinl_part0(2:end,1), size(Ygrid,1), ...
								size(Ygrid,2));
subplot(2,2,3);
H_ionburied_nl = contourf(Ygrid, Zgrid, ionBuried_nl,'linestyle','none');
axis equal;
colorbar;
set(gca,'fontsize',16,'position',[leftstart botstart imwidth imheight]);
title('(c)');
%set(gca,'fontsize',16);
%xlabel('Y');
%ylabel('Z');

ionBuried_l = reshape(phil_part0(2:end,1), size(Ygrid,1), ...
								size(Ygrid,2));
subplot(2,2,4);
H_ionburied_l = contourf(Ygrid, Zgrid, ionBuried_l,'linestyle','none');
axis equal;
colorbar;
set(gca,'fontsize',16,'position',[rightstart botstart imwidth imheight]);
title('(d)');

%xlabel('Y');
%ylabel('Z');