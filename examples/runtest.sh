#!/bin/sh
# this file should not be called directly from the shell,
# it ought to be included from within run.sh for each case

if [[ -x ./milonga ]] && [[ ! -h ./milonga ]]; then
 milongabin=./milonga
 testdir="examples/"
elif [[ -x ../milonga ]] && [[ ! -d ../milonga ]]; then
 milongabin="../milonga"
 testdir="./"
else
 echo "do not know how to run milonga :("
 exit 1
fi

# runs milonga
function runmilonga {
 $milongabin ${testdir}${1} ${2} ${3} ${4} ${5} ${6}
 outcome=$?
 if [ $outcome -ne 0 ]; then
   exit $outcome
 fi
}

# calls gnuplot with the provided command if it is installed
function plot {
 if [ ! -z "$2" ]; then
   format=$2
 else
   format=svg
 fi
 rm -f $1.$format
#  if [ "x`which pyxplot`" != "x" ]; then
#    echo "set terminal pdf; set output \"$1.$format\"" | pyxplot - $1.gp
 if [ "x`which gnuplot`" != "x" ]; then
   gnuplot -e "set terminal $format; set output \"$1.$format\"" $1.gp
 fi
 if [[ ! -z "$3" ]] && [[ -e $1.$format ]]; then
   xdg-open $1.$format
 fi
}

# checks if gmsh is installed
function checkgmsh {
 if [ -z "`which gmsh`" ]; then
  echo "gmsh is not installed, skipping test"
  exit 77
 fi
}

# checks if m4 is installed
function checkm4 {
 if [ -z "`which m4`" ]; then
  echo "m4 is not installed, skipping test"
  exit 77
 fi
}

# calls gmsh (if exists)
function callgmsh {
 if [ ! -z "`which gmsh`" ]; then
  gmsh $1 &
 fi
}


# calls pandoc (or markdown) and shows the result
function callpandoc {
 if [ ! -z "`which pandoc`" ]; then
   pandoc $1.txt -o $1.html -s
   pandoc $1.txt -o $1.pdf
 elif [ ! -z "`which markdown`" ]; then
   markdown $1.txt > $1.html
 fi

 if [ -f $1.html ]; then
   if [ ! -z "`which x-www-browser`" ]; then
     x-www-browser $1.html &
   elif [ ! -z "`which xdg-open`" ]; then
     xdg-open $1.html &
   fi
 fi
}
