// geometry taken from gmsh's t5.geo
rc = 7.428998;    // critical radius according to los alamos report
lc = rc/10;        // element characteristic length

p1 = newp; Point(p1) = {0,   0,   0,   lc} ;
p2 = newp; Point(p2) = {0+rc,0,   0,   lc} ;
p3 = newp; Point(p3) = {0,   0+rc,0,   lc} ;
p4 = newp; Point(p4) = {0,   0,   0+rc,lc} ;
p5 = newp; Point(p5) = {0-rc,0,   0,   lc} ;
p6 = newp; Point(p6) = {0,   0-rc,0,   lc} ;
p7 = newp; Point(p7) = {0,   0,   0-rc,lc} ;

c1 = newreg; Circle(c1) = {p2,p1,p7}; c2 = newreg; Circle(c2) = {p7,p1,p5};
c3 = newreg; Circle(c3) = {p5,p1,p4}; c4 = newreg; Circle(c4) = {p4,p1,p2};
c5 = newreg; Circle(c5) = {p2,p1,p3}; c6 = newreg; Circle(c6) = {p3,p1,p5};
c7 = newreg; Circle(c7) = {p5,p1,p6}; c8 = newreg; Circle(c8) = {p6,p1,p2};
c9 = newreg; Circle(c9) = {p7,p1,p3}; c10 = newreg; Circle(c10) = {p3,p1,p4};
c11 = newreg; Circle(c11) = {p4,p1,p6}; c12 = newreg; Circle(c12) = {p6,p1,p7};

// We need non-plane surfaces to define the spherical holes. Here we
// use ruled surfaces, which can have 3 or 4 sides:

l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1};
l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2};
l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3};
l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4};
l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5};
l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6};
l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7};
l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8};

// We then store the surface loops identification numbers in a list
// for later reference (we will need these to define the final volume):

theloops = newreg ;
Surface Loop(theloops) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

thehole = newreg ;
Volume(thehole) = theloops;

Physical Volume("fuel") = {30};
Physical Surface("external") = {28, 16, 14, 20, 24, 26, 18, 22};

Mesh.Algorithm = 5;
Mesh.RecombineAll = 0;
