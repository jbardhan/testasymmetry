//diel
R_diel = 4;
st_thickness = 2;
R = R_diel;
lc = 2;
lc_fine = 0.05;
hf = 50;
h0 = 0;
Point(1) = {R, 0, hf, lc};
Point(2) = {0, R, hf, lc};
Point(3) = {-R, 0, hf, lc};
Point(4) = {0, -R, hf, lc};
Point(5) = {0, 0, hf, lc_fine};

Circle(1) = {4, 5, 1};
Circle(2) = {1, 5, 2};
Circle(3) = {2, 5, 3};
Circle(4) = {3, 5, 4};
Extrude {0, 0, -hf + h0} {
  Line{4, 3, 2, 1};
}

Line Loop(23) = {3, 4, 1, 2};
Plane Surface(24) = {23};
Point{5} In Surface{24};
Line Loop(25) = {9, 5, 17, 13};
Plane Surface(26) = {25};
