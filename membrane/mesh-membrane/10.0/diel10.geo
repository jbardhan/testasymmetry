r = 10;
h = 50;
lc =2.214;
lc_fine =0.010;
st_thickness = 0;

h = h + st_thickness;
h0 = 0 - st_thickness;
Point(1) = {r, 0, h, lc};
Point(2) = {0, r, h, lc};
Point(3) = {-r, 0, h, lc};
Point(4) = {0, -r, h, lc};
Point(5) = {0, 0, h, lc_fine};

Circle(1) = {4, 5, 1};
Circle(2) = {1, 5, 2};
Circle(3) = {2, 5, 3};
Circle(4) = {3, 5, 4};
Extrude {0, 0, -h + h0} {
  Line{4, 3, 2, 1};
}

Line Loop(23) = {3, 4, 1, 2};
Plane Surface(24) = {23};
Point{5} In Surface{24};


Line Loop(25) = {5, 17, 13, 9};
Plane Surface(26) = {25};
