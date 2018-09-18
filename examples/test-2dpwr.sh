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

cat << EOF > 2dpwr-table.md
| Grid | Scheme | Matrix size | \$k_\\text{eff}\$ | \$\\max \\phi_2\$ | \$x\$ | \$y\$ | CPU [s] | Memory [Mb] |
| -----| ------ | ----------- | ----------------- | ----------------- | ----- | ----- | ------- | ----------- |
EOF

cat 2dpwr-table.md
i=0
for m in structured unstructured; do
  for s in volumes elements; do
    i=$(($i+1))
    runmilonga 2dpwr.mil $m $s | tee -a 2dpwr-table.md

    cat << EOF > 2dpwr-$m-$s.gp
set ticslevel 0
set view 30,70-10*$i
unset key
splot "2dpwr-$m-$s-3.dat" u 1:2:3 w p pt 37+2*$i ps 0.5 palette
EOF
    plot 2dpwr-$m-$s svg
  done
done

cat << EOF > 2dpwr.md
% IAEA 2D PWR Benchmark

# Multiplication factor and maximum flux

EOF
cat 2dpwr-table.md >> 2dpwr.md

for m in structured unstructured; do
  for s in volumes elements; do
  cat << EOF >> 2dpwr.md
  

## Power in $m $s

![](2dpwr-$m-$s.svg)

EOF
  done
done

if [ ! -z "`which pandoc`" ]; then
  pandoc --number-sections -s 2dpwr.md -o 2dpwr.html
  xdg-open 2dpwr.html
fi

if [ ! -z "$cded" ]; then
 cd ..
fi
            
exit 0
