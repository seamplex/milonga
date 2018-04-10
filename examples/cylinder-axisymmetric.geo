//
SetFactory("OpenCASCADE");

// radius of a critical sphere
rc = 7.428998;
// cylinder radius to maintain the same volume as the critical sphere
a = (2/3)^(1/3) * rc;

n = 8;
Mesh.CharacteristicLengthMin = a/n;
Mesh.CharacteristicLengthMax = a/n;

Rectangle (1) = {0, -a, 0, a, 2*a};
Extrude{ {0,1,0}, {0,0,0}, Pi }{ Surface{1}; Layers{2*n}; Recombine; }
Extrude{ {0,-1,0}, {0,0,0}, Pi }{ Surface{1}; Layers{2*n}; Recombine; }
Coherence;

Mesh.RecombineAll = 1;
Mesh.Algorithm = 8;

Physical Surface("external") = {2,3,5,6,7,8};
Physical Volume("fuel") = {1, 2};
