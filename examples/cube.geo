//
SetFactory("OpenCASCADE");

// radius of a critical sphere
rc = 7.428998;
// // cube half-length to maintain the same volume as the critical sphere
a = (4*Pi/(8*3))^(1/3) * rc;
lc = a/5;        // element characteristic length
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Box (1) = {-a, -a, -a, 2*a, 2*a, 2*a};

Physical Surface("external") = {1,2,3,4,5,6};
Physical Volume("fuel") = {1};


Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0;

Mesh.Algorithm = 6;    // (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad)
Mesh.Algorithm3D = 4;  // (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)

// Mesh.RecombineAll = 1;
// Mesh.RecombinationAlgorithm = 1;
// Mesh.Recombine3DAll = 1;
// Mesh.Recombine3DLevel = 2;
// Mesh.Recombine3DConformity = 4;
