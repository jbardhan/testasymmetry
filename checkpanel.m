jr1 =[522    -116.242751   -123.650472
804    -112.742302   -116.381854
1218   -111.318026   -113.428231
1644   -110.696043   -112.138193
2020   -110.521395   -111.701703
2418   -110.019520   -110.945229];


jr2 = [522    -45.064267     -49.630273
804    -40.661151     -42.919079
1218   -39.995102     -41.372785
1644   -39.639035     -40.505187
2020   -39.508701     -40.243670
2418   -39.266584     -39.853197];

np = jr1(:,1);

panelqual(1) =RichardsonExtrapolation(length(np)-1,length(np), np, jr1(:,2), -1);
panelqual(2) =RichardsonExtrapolation(length(np)-1,length(np), np, jr2(:,2), -1);

panelyl(1) =RichardsonExtrapolation(length(np)-1,length(np), np, jr1(:,3), -1);
panelyl(2) =RichardsonExtrapolation(length(np)-1,length(np), np, jr2(:,3), -1);

