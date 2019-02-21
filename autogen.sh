#!/bin/sh
# 
# Execute this script to generate a configure script
#
# This file is free software: you are free to change and redistribute it.
# There is NO WARRANTY, to the extent permitted by law.
#
plugin=milonga

# check for needed tools (cannot put this into an m4 macro
# because we do not even know if we have m4 available)
for i in git m4 autoconf xargs; do
 if [ -z "`which $i`" ]; then
  echo "error: $i not installed"
  exit 1
 fi
done

./autoclean.sh

if [ -e ../wasora/.git ]; then
  WASORA_REPO=../wasora/.git
else 
  WASORA_REPO=https://bitbucket.org/seamplex/wasora.git
fi

rm -rf wasora
git clone ${WASORA_REPO} || exit 1

m4 wasora/m4/bootstrap.m4 - << EOF | sh -s $1 || exit 1
plugin=${plugin}
WASORA_CHECK_VCS
PLUGIN_VERSION_VCS
WASORA_README_INSTALL
EOF

# build makefile.am
am="src/Makefile.am"
echo -n "building $am... "
cat << EOF > $am
include \$(SLEPC_DIR)/lib/slepc/conf/slepc_variables
EOF
# this variable is set by petsc somewhere and we need it empty to make install
# apparently undefine DESTDIR is supported by GNU Make >= 3.82
makever=`make -v | head -n1 | cut -d' ' -f3`
makemajor=`echo $makever | cut -d. -f1`
makeminor=`echo $makever | cut -d. -f2`
if [ $makemajor -ge 4 -o $makeminor -ge 82 ]; then
  echo unexport DESTDIR >> $am
fi

cat << EOF >> $am
AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = \$(ACLOCAL_FLAGS)

bin_PROGRAMS = ${plugin}

${plugin}_CFLAGS = -I../wasora/src \$(SLEPC_CC_INCLUDES) \$(PETSC_CC_INCLUDES) \$(CC_INCLUDES) \$(all_includes) -DHARDCODEDPLUGIN 
${plugin}_LDADD = \$(STANDALONELIBS) \$(SLEPC_LIB) \$(all_libraries)
${plugin}_LDFLAGS = -rdynamic

${plugin}_SOURCES = \\
EOF

cd src
find . -maxdepth 2 \( -name "*.c" -o -name "*.h" \) | xargs echo -n >> ../$am
echo "\\" >> ../$am
find ../wasora/src -maxdepth 2 \( -name "*.c" -o -name "*.h" \) | xargs echo >> ../$am
cd ..

cat << EOF >> $am

version.\$(OBJEXT): version.h
version.h: Makefile
	./version.sh
EOF
echo "done"

echo "calling autoreconf... "
autoreconf -i
echo "done"

