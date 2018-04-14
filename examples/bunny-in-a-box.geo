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
p2() = PointsOf{ Volume{2}; };
p3() = PointsOf{ Volume{3}; };
Characteristic Length { p2() } = lc;
Characteristic Length { p3() } = lc;
Mesh.CharacteristicLengthMin = 0.5 * lc;
Mesh.CharacteristicLengthMax = 2.0 * lc;

Mesh.Algorithm = 6;
Mesh.Algorithm3D = 2;


// remember to scale the resulting mesh to fit the critical volume!
// $ cp bunny-in-a-box.msh bunny_orig.msh
// $ milonga bunnyscale.mil
// $ mv bunnynew.msh bunny-in-a-box.msh
