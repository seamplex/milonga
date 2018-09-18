plugin = milonga

SUBDIRS = src
DIST_SUBDIRS = src wasora

ACLOCAL_AMFLAGS = -I wasora/m4
dist_doc_DATA = AUTHORS ChangeLog README TODO 
EXTRA_DIST = src/version.sh locateruntest.sh examples \
             doc/README.md \
             doc/README.html  doc/README.pdf \
             doc/INSTALL.html doc/INSTALL.pdf

TESTS = examples/test-slab.sh \
        examples/test-square.sh \
        examples/test-2dpwr.sh
#         examples/test-3dshapes.sh
#         examples/test-3dpwr-struct.sh

all-local:
	if [ -e src/.libs/lib${plugin}.so.0.0.0 ]; then cp -u  src/.libs/lib${plugin}.so.0.0.0 ./${plugin}.so; fi
	if [ -e src/${plugin} ]; then cp -u  src/${plugin} .; fi
	if [ -e ${plugin} ]; then cd examples; ln -sf ../${plugin}; cd;  fi
	rm -f src/version.$(OBJEXT) src/version.h src/version.lo src/.lib/version.$(OBJEXT)

clean-local:
	rm -f src/wasora/${plugin} examples/${plugin}
	rm -f ${plugin} ${plugin}.so
	rm -f *~ examples/*~
	rm -f examples/*.dat examples/*.pos examples/*.vtk
	rm -f examples/*pwr*.txt examples/*pwr*.msh
	rm -f examples/*pwr*.pdf examples/../*pwr*.pdf 
	rm -f examples/*pwr*.html
	rm -f *.log *.trs
check-local:
	@echo
	@echo See the README and the directory examples for details about the test suite.
	@echo
