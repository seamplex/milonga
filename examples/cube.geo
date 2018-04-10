//
SetFactory("OpenCASCADE");

/*
Rectangle (1) = {0, 0, 0, 1, 1};
Extrude {0,0,1} { Surface{1}; Layers{10}; Recombine; }
*/

// radius of a critical sphere
rc = 7.428998;
// // cube half-length to maintain the same volume as the critical sphere
a = (4*Pi/(8*3))^(1/3) * rc;
lc = a/10;        // element characteristic length
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Box (1) = {-a, -a, -a, 2*a, 2*a, 2*a};

Physical Surface("external") = {1,2,3,4,5,6};
Physical Volume("fuel") = {1};


Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0;

Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 1;
Mesh.Recombine3DAll = 1;
Mesh.Algorithm = 8;
Mesh.Algorithm3D = 2;
Mesh.Recombine3DLevel = 2;
Mesh.Recombine3DConformity = 4;
