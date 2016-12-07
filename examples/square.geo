a = 100;
lc = a/20;

Point(1) = {0, 0, 0, lc};
Point(2) = {a, 0, 0, lc};
Point(3) = {a, a, 0, lc};
Point(4) = {0, a, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Physical Line("left") = {4};
Physical Line("right") = {2};
Physical Line("front") = {1};
Physical Line("back") = {3};

Physical Surface("fuel") = {6};

Mesh.Algorithm = 5;
//Mesh.RecombineAll = 1;
