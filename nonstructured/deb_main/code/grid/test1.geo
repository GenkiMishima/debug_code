lc = 0.5;
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
lc2 = 0.1;
Point(5) = {0.3, 0.3, 0, lc2};
Point{5} In Surface{6};
Physical Surface(7) = {6};
Physical Line(8 ) = {1};
Physical Line(9 ) = {2};
Physical Line(10) = {3};
Physical Line(11) = {4};