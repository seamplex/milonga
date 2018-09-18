//
SetFactory("OpenCASCADE");

// radius of a critical sphere
rc = 7.428998;
// cylinder radius to maintain the same volume as the critical sphere
a = (2/3)^(1/3) * rc;

lc = a/4;        // element characteristic length

Cylinder(1) = {0,0,-a, 0,0,2*a, a};

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Physical Volume("fuel") = {1};
Physical Surface("external") = {1,2,3};

Mesh.Algorithm = 6;    // (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad)
Mesh.Algorithm3D = 1;  // (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
