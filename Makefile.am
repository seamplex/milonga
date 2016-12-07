plugin = milonga

SUBDIRS = src $(DIR_STANDALONE)
if WANT_STANDALONE
  SUBDIRS += src/wasora
endif
DIST_SUBDIRS = src src/wasora

ACLOCAL_AMFLAGS = -I m4
dist_doc_DATA = AUTHORS ChangeLog README PLUGINS TODO 
EXTRA_DIST = src/version.sh locateruntest.sh examples \
             doc/README.md \
             doc/README.html  doc/README.pdf \
             doc/INSTALL.html doc/INSTALL.pdf \
             doc/PLUGINS.html doc/PLUGINS.pdf

TESTS = examples/test-slab.sh \
        examples/test-square.sh \
        examples/test-2dpwr.sh \
        examples/test-3dshapes.sh
#         examples/test-3dpwr-struct.sh

all-local:
	if [ -e src/.libs/lib${plugin}.so.0.0.0 ]; then cp -u  src/.libs/lib${plugin}.so.0.0.0 ./${plugin}.so; fi
	if [ -e src/wasora/${plugin}$(EXEEXT) ]; then cp -u  src/wasora/${plugin}$(EXEEXT) .; fi
	if [ -e ${plugin}$(EXEEXT) ]; then cd examples; ln -sf ../${plugin}$(EXEEXT); cd;  fi
	rm -f src/version.$(OBJEXT) src/version.h src/version.lo src/.lib/version.$(OBJEXT)

clean-local:
	rm -f src/wasora/${plugin}$(EXEEXT) examples/${plugin}$(EXEEXT)
	rm -f ${plugin}$(EXEEXT) ${plugin}.so
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
