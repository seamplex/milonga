#!/bin/bash
# run the 2D IAEA PWR benchmark over an unstructured mesh of
# characteristic length lc=3.5cm with finite elements

. locateruntest.sh

checkgmsh
checkm4

# TODO: fix this!
if [ ! -e 2dpwr.geo.m4 ]; then
 cded=1
 cd examples
 milongabin="../milonga"
 testdir="./"
fi

runmilonga 2dpwr.mil structured    volumes
# callgmsh 2dpwr-structured-volumes-3.msh

runmilonga 2dpwr.mil structured    elements
# callgmsh 2dpwr-structured-elements-3.msh

runmilonga 2dpwr.mil unstructured  volumes
# callgmsh 2dpwr-unstructured-volumes-3.msh

runmilonga 2dpwr.mil unstructured  elements
# callgmsh 2dpwr-unstructured-elements-3.msh

plot "set title 'five 2D PWR solutions'; \
      set ticslevel 0; \
      splot '2dpwr-structured-volumes-3.dat'   u 1:2:5 w lp pt 2 palette, \
            '2dpwr-structured-elements-3.dat'  u 1:2:5 w lp pt 3 palette,\
            '2dpwr-unstructured-volumes-3.dat' u 1:2:5 w p pt 5 palette,\
            '2dpwr-unstructured-elements-3.dat'u 1:2:5 w p pt 6 palette"

if [ ! -z "$cded" ]; then
 cd ..
fi
            
exit 0
