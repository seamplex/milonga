//
SetFactory("OpenCASCADE");
a = 100;

Rectangle(1) = {-a/2, -a/2, 0, a, a};
Physical Surface("fuel", 1) = {1};
Physical Line("external") = {1,2,3,4};

Mesh.CharacteristicLengthMax = a/20;
Mesh.CharacteristicLengthMin = Mesh.CharacteristicLengthMax;

Mesh.Algorithm = 6;
//Mesh.RecombineAll = 1;
//Mesh.ElementOrder = 2;
