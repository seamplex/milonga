//
SetFactory("OpenCASCADE");

// a = 160;
// b = 260;
// lc = 20;

Sphere(1) = {0,0,0, a, 0, Pi/2, Pi/2};
Sphere(2) = {0,0,0, b, 0, Pi/2, Pi/2};
BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; };
Coherence;

Physical Volume("fuel") = {1};
Physical Volume("reflector") = {3};
Physical Surface("internal") = {14, 10, 16, 12, 11, 15};
Physical Surface("external") = {13};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Mesh.Algorithm = 6;
Mesh.Algorithm3D = 2;
Mesh.Optimize = 1;
