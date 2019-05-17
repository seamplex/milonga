lc = 1;
f = 0.1;

Point(1) = {0, 0, 0, f*lc};
Point(2) = {2, 0, 0, f*lc};
Point(3) = {3, 0, 0, f*lc};
Point(4) = {5, 0, 0, f*lc};
Point(5) = {6, 0, 0, f*lc};
Point(6) = {8, 0, 0, f*lc};

Point(11) = {1,   0, 0, lc};
Point(12) = {2.5, 0, 0, lc};
Point(13) = {4,   0, 0, lc};
Point(14) = {5.5, 0, 0, lc};
Point(15) = {7,   0, 0, lc};

Line(1) = {1, 11};
Line(11) = {11, 2};
Line(2) = {2, 12};
Line(12) = {12, 3};
Line(3) = {3, 13};
Line(13) = {13, 4};
Line(4) = {4, 14};
Line(14) = {14, 5};
Line(5) = {5, 15};
Line(15) = {15, 6};


Physical Line("source_abs")  = {1,11};
Physical Line("absorber")    = {2,12};
Physical Line("void")        = {3,13};
Physical Line("source_scat") = {4,14};
Physical Line("reflector")   = {5,15};

Physical Point("mirror") = {1};
Physical Point("vacuum") = {6};

