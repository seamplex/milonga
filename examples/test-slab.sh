#!/bin/bash
# solve a one-dimensional slab with different schemes

. locateruntest.sh

# we run the four cases with the single-mesh input
runmilonga slab-struct.mil   --volumes  | tee slab-struct-vol.dat
runmilonga slab-struct.mil   --elements | tee slab-struct-ele.dat
runmilonga slab-unstruct.mil --volumes  | tee slab-unstruct-vol.dat
runmilonga slab-unstruct.mil --elements | tee slab-unstruct-ele.dat

# we check that the output coincides with the two-mesh input results
runmilonga slab.mil structured   --volumes  > slab-struct-vol2.dat
runmilonga slab.mil structured   --elements > slab-struct-ele2.dat
runmilonga slab.mil unstructured --volumes  > slab-unstruct-vol2.dat
runmilonga slab.mil unstructured --elements > slab-unstruct-ele2.dat

cat << EOF > slab.gp
set title "Five slab solutions a-la-milonga"
set key above
set grid
a = 2*10.371065
set xrange [0:a]
plot 'slab-struct-vol.dat'   w p pt 2 lt 1,\
     'slab-struct-ele.dat'   w p pt 4 lt 2,\
     'slab-unstruct-vol.dat' w p pt 5 lt 3,\
     'slab-unstruct-ele.dat' w p pt 6 lt 4,\
     pi/2*sin(pi*x/a) w l lt 7
EOF
plot slab pdf show
exit 0
