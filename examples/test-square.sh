#!/bin/bash
# solve a one-dimensional slab with different schemes

. locateruntest.sh

runmilonga square.mil structured   --volumes  > square-struct-vol.dat
mv square_out.msh square-struct-vol.msh
mv square.vtk square-struct-vol.vtk
# callgmsh square-struct-vol.msh

runmilonga square.mil structured   --elements > square-struct-ele.dat
mv square_out.msh square-struct-ele.msh
mv square.vtk square-struct-ele.vtk
# callgmsh square-struct-ele.msh

runmilonga square.mil unstructured --volumes  > square-unstruct-vol.dat
mv square_out.msh square-unstruct-vol.msh
mv square.vtk square-unstruct-vol.vtk
# callgmsh square-unstruct-vol.msh

runmilonga square.mil unstructured --elements > square-unstruct-ele.dat
mv square_out.msh square-unstruct-ele.msh
mv square.vtk square-unstruct-ele.vtk
# callgmsh square-unstruct-ele.msh

# plot "set title 'five square solutions'; \
#       set ticslevel 0; \
#       splot 'square-struct-vol.dat'   w lp pt 2 palette, \
#             'square-struct-ele.dat'   w lp pt 4 palette,\
#             'square-unstruct-vol.dat' w p pt 5 palette,\
#             'square-unstruct-ele.dat' w p pt 6 palette,\
#             (pi/2)**2*sin(pi*x/100)*sin(pi*y/100) w l palette"

exit 0
