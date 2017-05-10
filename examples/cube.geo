//
// radius of a critical sphere
rc = 7.428998;
// // cube half-length to maintain the same volume as the critical sphere
a = (4*Pi/(8*3))^(1/3) * rc;

// a = 50;
lc = a/10;        // element characteristic length

Point(1) = {-a, -a, -a, lc};
Point(2) = {+a, -a, -a, lc};
Point(3) = {+a, +a, -a, lc};
Point(4) = {-a, +a, -a, lc};

Point(11) = {-a, -a, +a, lc};
Point(12) = {+a, -a, +a, lc};
Point(13) = {+a, +a, +a, lc};
Point(14) = {-a, +a, +a, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 11};

Line(15) = {1, 11};
Line(16) = {2, 12};
Line(17) = {3, 13};
Line(18) = {4, 14};

Line Loop(19) = {-15, -11, +16, +1};
Plane Surface(20) = {19};
Line Loop(21) = {-16, -12, +17, +2};
Plane Surface(22) = {21};
Line Loop(23) = {-17, -13, +18, +3};
Plane Surface(24) = {23};
Line Loop(25) = {4, 15, -14, -18};
Plane Surface(26) = {25};
Line Loop(27) = {11, 12, 13, 14};
Plane Surface(28) = {27};
Line Loop(29) = {-1, -2, -3, -4};
Plane Surface(30) = {29};

Surface Loop(31) = {20, 26, 30, 24, 22, 28};
Volume(32) = {31};

Mesh.Algorithm = 5;
Mesh.RecombineAll = 0;

Physical Surface("external") = {20, 22, 24, 26, 30, 28};
Physical Volume("fuel") = {32};
