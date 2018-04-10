//
SetFactory("OpenCASCADE");

// geometry taken from gmsh's t5.geo
rc = 7.428998;    // critical radius according to los alamos report
lc = rc/5;        // element characteristic length

Torus(1) = {0,0,0, 2*rc, 0.75*rc};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Physical Volume("fuel") = {1};
Physical Surface("external") = {1};

Mesh.Algorithm = 6;
Mesh.Algorithm3D = 2;
