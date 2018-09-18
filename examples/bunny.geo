//
SetFactory("OpenCASCADE");

a() = ShapeFromFile("bunny.brep");

// radius of a critical sphere with the XSs from the problem
rc = 7.428998;
critical_volume = 4/3*Pi*rc^3;

// according to freecad, the original bunny volume is
// App.ActiveDocument.ActiveObject.Shape.Volume = 128833.65653379083
original_volume = 128833.65653379083;
scale = (critical_volume/original_volume)^(1/3);

Dilate {{0, 0, 0}, scale} {Volume{1};}

Mesh.CharacteristicLengthMin = 10;
Mesh.CharacteristicLengthMax = 10;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.Algorithm = 5; //  (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad)

Physical Surface("external") = {1:700};
Physical Volume("fuel") = {1};
