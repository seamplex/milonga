//
SetFactory("OpenCASCADE");

rc = 7.428998;    // critical radius according to los alamos report
lc = rc/5;        // element characteristic length

Sphere(1) = {0,0,0, rc};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Physical Volume("fuel") = {1};
Physical Surface("external") = {1};

Mesh.Algorithm = 6;    // (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad)
Mesh.Algorithm3D = 4;  // (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
