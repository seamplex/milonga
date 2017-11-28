SetFactory("OpenCASCADE");

a = 10;
n = 10;
Mesh.CharacteristicLengthMin = a/n;
Mesh.CharacteristicLengthMax = a/n;

Rectangle (1) = {0, -a, 0, a, 2*a};
Extrude{ {0,1,0}, {0,0,0}, Pi }{ Surface{1}; Layers{2*n}; Recombine; }
Extrude{ {0,-1,0}, {0,0,0}, Pi }{ Surface{1}; Layers{2*n}; Recombine; }
Coherence;

Mesh.RecombineAll = 1;
Mesh.Algorithm = 8;

Physical Surface("external") = {4, 8, 7, 3, 2, 6};
Physical Volume("fuel") = {1, 2};
