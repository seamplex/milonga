//
SetFactory("OpenCASCADE");


d = 1;
Box(1) = { -45-d, -35-d, -1.6-d, 2*(45+d), 2*(35+d), 1.6+87+2*d};
a() = ShapeFromFile("bunny.brep");

BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; };
Coherence;

Physical Surface("external") = {709, 710, 712, 707, 708, 711};
Physical Volume("fuel") = {2};
Physical Volume("vacuum") = {3};

lc = 20;
Mesh.CharacteristicLengthMin = 0.5 * lc;
Mesh.CharacteristicLengthMax = 2.0 * lc;

Mesh.Algorithm = 6;
