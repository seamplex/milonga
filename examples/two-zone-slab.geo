a = 100;
b = 50;
lc = 2;

Point(1) = {0, 0, 0, lc};
Point(2) = {b, 0, 0, lc};
Point(3) = {a, 0, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Physical Line("A") = {1};
Physical Line("B") = {2};

Physical Point("left") = {1};
Physical Point("right") = {3};
