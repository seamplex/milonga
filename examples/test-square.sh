#!/bin/bash
# solve a two-dimensional square

. locateruntest.sh

rm -f square-*.dat square.gp square.svg square-*.png square-*.geo square.html

runmilonga square.mil --volumes | tee square-fvm.dat
mv square_phi.msh square-fvm.msh
mv square_phi.vtk square-fvm.vtk

runmilonga square.mil --elements | tee square-fem.dat
mv square_phi.msh square-fem.msh
mv square_phi.vtk square-fem.vtk

cat << EOF > square.gp
set title "Two square solutions a-la-milonga"
set ticslevel 0
splot "square-fvm.dat", "square-fem.dat"
EOF
plot square svg

if [ ! -z "`which gmsh`" ]; then
  for i in fvm fem; do
    cat << EOF > square-$i.geo
Merge "square-$i.msh";
General.SmallAxes = 0;
Print "square-$i.png";
Exit;
EOF
    gmsh square-$i.geo
    convert -trim square-$i.png square-$i.png
  done
fi

cat << EOF > square.md
% FEM vs. FVM

A bare square solved with [milonga](https://www.seamplex.com):

![Flux map computed with finite volumes](square-fvm.png)

![Flux map computed with finite elements](square-fem.png)

![Both solutions as $\phi(x,y)$](square.svg)

# Terminal mimic

Here is the input and the commandline used to run the case:

~~~
$ cat square.mil
`cat square.mil`
$ milonga square.mil
~~~


EOF
if [ ! -z "`which pandoc`" ]; then
  pandoc -s square.md -o square.html
  xdg-open square.html
fi

exit 0
