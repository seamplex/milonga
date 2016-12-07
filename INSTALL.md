Installing milonga
==================

This file contains brief instructions to download, and/or compile and/or install [milonga](http://www.talador.com.ar/jeremy/wasora/milonga/). The detailed discussed is deferred to the full documentation (see directory `doc`).

Milonga works on top of the framework provided by [wasora](http://www.talador.com.ar/jeremy/wasora). Actually, milonga is a plugin for wasora that can be dynamically loaded at run-time (the set of wasora plus one its plugins is referred to as the _wasora suite_). Nevertheless, milonga can be compiled and executed as a stand-alone binary. In any case, milonga follows the same design principles that are built into wasora. See the wasora [home page](http://www.talador.com.ar/jeremy/wasora) and [source repository](https://bitbucket.org/gtheler/wasora) for further details.

Getting milonga
---------------

The wasora code can be obtained essentially in either source or binary form. Sinces the v0.3.x series, milonga uses [Mercurial](http://mercurial.selenic.com/) as the version control system and the repository is hosted at [Bitbucket](https://bitbucket.org/gtheler/milonga) (v0.2.x series used [Bazaar](http://bazaar.canonical.com/en/) and [Launchpad](https://launchpad.net/)]). Tarballs containing either sources or binaries are periodically prepared from the repository sources and may be downloaded from <https://bitbucket.org/gtheler/milonga/downloads>.

In order of increasing level of user expertise, to get milonga either

 * Download the latest binary tarball for your architecture as listed in milonga's [home page](http://www.talador.com.ar/jeremy/wasora/milonga). Currently, the available options are

    - linux-amd64: GNU/Linux 64-bit binary executable

    - linux-i686: GNU/Linux 32-bit binary executable

    - cygwin-i686: Cygwin-based 32-bit binary Windows executable

   and proceed to the [Installing] section.

   The provided binaries are statically linked to the required libraries to avoid having a user that expects to run milonga out of the box dealing with unresolved library dependences. However, there may be some libraries that are not available for certain configurations. A full source compilation is recommended.

   > Note that using wasora in either non-free and/or non-UNIX operating systems is highly discouraged.  Binaries and support for other platforms are provided only as a polite way to allow potential users to be introduced to milonga. Real usage of scientific and engineering computations codes _ought_ to be performed in UNIX-based operating systems. Please switch to [GNU/Linux](http://www.debian.org/).

 * Download the latest source tarball listed in milonga's [home page](http://www.talador.com.ar/jeremy/wasora/milonga). A list of files and versions available to be downloaded can be found by browsing

     <http://www.talador.com.ar/jeremy/wasora/download>

   Then proceed to the [Compiling] section.

 * Clone the [Mercurial](http://mercurial.selenic.com/) [Bitbucket](https://bitbucket.org/gtheler/wasora) repository:

        $ hg clone https://bitbucket.org/gtheler/milonga

   and proceed to the [Bootstrapping] section. You may keep your tree up-to-date by pulling incremental changes:

        $ cd milonga
        $ hg pull
        pulling from https://bitbucket.org/gtheler/milonga
        searching for changes
        no changes found
        $ 


Bootstrapping
-------------

Skip this section if you did not clone the [repository](https://bitbucket.org/gtheler/wasora).

The repository development tree has to be bootstrapped by [autoconf & friends](http://en.wikipedia.org/wiki/GNU_build_system) to be able to configure and make the code. As with wasora, the script `autogen.sh` generates the files that `autoconf` needs to produce a working `configure` script (see wasora's `INSTALL` file). However, being a plugin, milonga needs access to the wasora source tree either to compile a standalone binary or the access the required header files to build a loadable shared-object file. Note that `libtool` is also needed, but `autoreconf` will fail with a non-clear error message if it not installed. It the script `autogen.sh` complains with a cryptic message, be sure to check if `libtool` is properly installed.

As this tree generates a plugin for [wasora](http:/www.talador.com.ar/jeremy/wasora), access its source tree is needed. By default, `autogen.sh` tries to find a valid wasora source tree (either an `hg`-cloned or an uncompressed source tarball) in the following relative locations:

 * `../wasora`
 * `../../wasora`
 * `../../../wasora`

A different location can be provided either

 1. by passing the path as an argument of `autogen.sh`

        $ ./autogen.sh $HOME/wasora

 2. by setting the environment variable `WASORA_DIR`

        $ export WASORA_DIR=$HOME/wasora
        $ ./autogen.sh

The wasora source directory `src` is copied into the plugins tree under `src/wasora`. The associated script `autoclean.sh` should remove any automatically-generated files (including `src/wasora`) and revert the tree to a freshly `hg`-cloned status.

$ ./autogen.sh 
cleaning... done
getting hg revision id... done
building changelog... done
formatting readme & install... done
building configure.ac... done
building src/Makefile.am... done
calling autoreconf... 
configure.ac:21: installing './compile'
configure.ac:16: installing './config.guess'
configure.ac:16: installing './config.sub'
configure.ac:18: installing './install-sh'
configure.ac:18: installing './missing'
parallel-tests: installing './test-driver'
src/Makefile.am: installing './depcomp'
done
$ 

At this point, a tree similar to the source distribution tarball is obtained, which can be configured and compiled as described in the section [Compiling] below. The bootstrapped files are listed in `.hgignore` so Mercurial should not report any changes in the status of the working tree after executing `autogen.sh`:

    $ hg status
    $

To clean the working tree and revert it to a fresh-clone status, the `autoclean.sh` script should be executed:

    $ ./autoclean.sh

Note that `autogen.sh` calls `autoclean.sh` first.



After a successful execution of `autogen.sh`, the `configure` script should be ready to be executed. 


Compiling
---------

Skip this section if you downloaded a binary tarball.

wasora (and thus milonga) follows the standard GNU `./configure && make` procedure. So uncompress the downloaded tarball into a proper location within your home directory:


    $ tar xvzf milonga-0.4.14.tar.gz
    $ cd wasora

### PETSc and SLEPc

milonga relies on the [SLEPc](http://www.grycap.upv.es/slepc/) library to solve the eigenvalue problem, which in turn depends on the [PETSc](http://www.mcs.anl.gov/petsc/) library. Therefore, they should be correctly installed and available before compiling milonga. At least version 3.5.0 is needed to compile milonga, so chances are that they have to be downloaded and compiled from scratch (i.e. distribution repositories are known to have previous versions).

If you already have installed them and are familiar with PETSc and SLEPc, be sure to correctly set the environment variables `PETSC_DIR`, `PETSC_ARCH` and `SLEPC_DIR` and then proceed to the [milonga] subsection below.  
If you are not familiar with PETSc and SLEPc, first check this URL

<http://www.mcs.anl.gov/petsc/documentation/installation.html>

and then follow this quick-start steps:

    $ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.6.2.tar.gz
    $ tar xvzf petsc-3.6.2.tar.gz
    $ cd petsc-3.6.2
    $ export PETSC_DIR=$PWD
    $ export PETSC_ARCH=arch-linux2-c-opt
    $ ./configure --with-fc=0 --with-mpi=0 --with-debugging=0 --with-x=0 --download-cblaslapack
    $ make all
    $ cd ..
    $ wget http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.6.2.tar.gz -O slepc-3.6.2.tar.gz
    $ tar xvzf slepc-3.6.2.tar.gz
    $ cd slepc-3.6.2
    $ export SLEPC_DIR=$PWD
    $ ./configure
    $ make
    $ cd ..
    $

Remember that exported environment variables last as long as the session in which they were exported is open. To have them always defined, add the lines 

    export PETSC_DIR=$HOME/petsc-3.6.2
    export SLEPC_DIR=$HOME/slepc-3.6.2
    export PETSC_ARCH=arch-linux2-c-opt

to your `~/.bashrc` file. See PETSc's and SLEPc's full documentation for further details.


### milonga

Execute the `configure` script so it can check the environment is able to build milonga:

```
$ ./configure
checking build system type... x86_64-unknown-linux-gnu
checking host system type... x86_64-unknown-linux-gnu
checking for a BSD-compatible install... /usr/bin/install -c
checking whether build environment is sane... yes
[...]
config.status: creating src/Makefile
config.status: executing depfiles commands

## --------------------- ##
## Configuration summary ##
## --------------------- ##
  GSL library (required): yes, version 1.16

  IDA library (optional): yes, version unknown
    differential-algebraic systems will be solved

  Readline library (opt): yes, version 6.3
    run-time debugging-like capabilities will be provided

Now proceed to compile with 'make'

$
```


By default, both a shared-object file with the dynamically-loadable plugin (`milonga.so`) and a standalone executable of wasora with the plugin statically linked into it (`milonga`). To disable the generation of either one, `configure` provides the options `--disable-plugin` or `--disable-standalone`.


Once the `configure` step successfully tells us to proceed to compile with `make`, that is what we do:

    $ make

In some cases, if the `dash` shell is used there may be some problems with the complex compilation lines generated by libtool. The `bash` shell is smarter and ought to be used:

    $ make SHELL=/bin/bash

The executable binary will be located in the current directory and called `milonga`. At this point you may want to check if milonga actually works by executing the test suite:

    $ make check

If [gnuplot](http://www.gnuplot.info/) and/or [Gmsh](http://geuz.org/gmsh/) are installed, some graphical windows should pop up. See the `README` for a full list of the test involved.

The usual way to finish the compilation of a program that follows the GNU standard is to perform a system-wide installation by executing (as root):

    # make install

This would leave the executable `milonga` available to be executed by any user of the system. However, other workflows may be used to run wasora. See the section [Installing] for details.


Required libraries
------------------

The code depends on a few libraries (and its development headers). Some of them are mandatory and some of the are optional. The name in parenthesis refers to the Debian-based package name

 * Mandatory for compilation:
     - [GNU Scientific Library](http://www.gnu.org/software/gsl/) (`libgsl0-dev`)

 * Needed if DAE systems are to be solved:
     - [IDA SUNDIALS Library](http://computation.llnl.gov/casc/sundials/main.html) (`libsundials-serial-dev`)

 * Needed if debug-mode is desired:
     - [GNU Readline](http://cnswww.cns.cwru.edu/php/chet/readline/rltop.html) (`libreadline-dev`)

Therefore, in Debian-based GNU/Linux boxes, one would do

    # apt-get install libgsl0-dev libsundials-serial-dev libreadline-dev

and all the required libraries (and development headers) should be detected by configure. Note that some wasora plugins (such as [milonga](http://www.talador.com.ar/jeremy/wasora/milonga)) may need further additional libraries (for instance [PETSc](http://www.mcs.anl.gov/petsc/) and [SLEPc](http://www.grycap.upv.es/slepc/)).

If `configure` is still unable to detect the GSL, it can be instructed to download and compile it in a local subdirectory using the `--enable-download-gsl` option:

    $ ./configure --enable-download-gsl

If no Internet connection is available, the file `gsl-1.16.tar.gz` may be separately obtained (for example from <http://ftpmirror.gnu.org/gsl/>) and copied into the wasora directory.

When giving the `--enable-download-gsl` option, the generated binary will be statically linked against the downloaded library. 


By default, `configure` checks for the optional libraries. However, they can be explicitly disabled by using the `--without-ida` and `--without-readline` options in `configure`:

    $ ./configure --without-ida --without-readline 
    [...]
    config.status: creating src/Makefile
    config.status: executing depfiles commands

    ## --------------------- ##
    ## Configuration summary ##
    ## --------------------- ##
      GSL library (required): yes, version 1.16

      IDA library (optional): no
        differential-algebraic systems will NOT be solved

      Readline library (opt): no
        run-time debugging-like capabilities will NOT be provided

    WARNING: there is at least one optional library missing.
    If this was not the desired result, check config.log for clues.

    Now proceed to compile with 'make'

    $

Call `./configure --help` to see all the available options. Which libraries wasora was finally linked against to can be checked by executing it with the `-v` option:

    $ ./milonga -v
    milonga 0.4.14  (10996c271131 2015-11-18 08:02 -0300)
    free nuclear reactor core analysis code
    
     rev hash 10996c271131758544a4b02123420fdb51ed1ff2
     last commit on 2015-11-18 08:02 -0300 (rev 235)
     compiled on 2015-11-19 11:35:38 by gtheler@frink (linux-gnu x86_64)
     with gcc (Debian 4.9.2-10) 4.9.2 using -O2 linked against
      SLEPc Release Version 3.6.2, Nov 03, 2015
      Petsc Release Version 3.6.2, Oct, 02, 2015  arch-linux2-c-opt
     running on Linux 3.16.0-4-amd64 #1 SMP Debian 3.16.7-ckt11-1+deb8u6 (2015-11-09) x86_64
     8  Intel(R) Core(TM) i7-3770K CPU @ 3.50GHz
    
    
     milonga is copyright (c) 2010-2015 jeremy theler
     licensed under GNU GPL version 3 or later.
     milonga is free software: you are free to change and redistribute it.
     There is NO WARRANTY, to the extent permitted by law.
    
    
    ---------------             ----------          ---------
    wasora 0.4.19  (040eac0f20eb 2015-11-14 11:46 -0300)
    wasora's an advanced suite for optimization & reactor analysis
    
     rev hash 040eac0f20ebbac329622568089b98582812bbd2
     last commit on 2015-11-14 11:46 -0300 (rev 173)
    
     compiled on 2015-11-17 08:48:31 by gtheler@frink (linux-gnu x86_64)
     with gcc (Debian 4.9.2-10) 4.9.2 using -O2 and linked against
      GNU Scientific Library version 1.16
      GNU Readline version 6.3
      SUNDIALs Library version 2.5.0
    
    
     wasora is copyright (C) 2009-2015 jeremy theler
     licensed under GNU GPL version 3 or later.
     wasora is free software: you are free to change and redistribute it.
     There is NO WARRANTY, to the extent permitted by law.
    $ 
 

If something goes wrong and the compilation fails, please feel free to ask for help at the wasora mailing list at <wasora@talador.com.ar>.
 


Further information
-------------------

Home page: <http://www.talador.com.ar/jeremy/wasora/milonga>  
Mailing list and bug reports: <wasora@talador.com.ar>  


milonga is copyright (C) 2009--2015 jeremy theler  
milonga is licensed under [GNU GPL version 3](http://www.gnu.org/copyleft/gpl.html) or (at your option) any later version.  
milonga is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  
See the file `COPYING` for copying conditions.  

The text below the cutting line corresponds to the original FSF instructions for installing software (as wasora) that follows the GNU configure & make convention.




-----------------------------------------------------------------------


Installation Instructions
=========================

Copyright (C) 1994-1996, 1999-2002, 2004-2012 Free Software Foundation,
Inc.

   Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without warranty of any kind.

Basic Installation
------------------

   Briefly, the shell commands `./configure; make; make install` should
configure, build, and install this package.  The following
more-detailed instructions are generic; see the `README` file for
instructions specific to this package.  Some packages provide this
`INSTALL` file but do not implement all of the features documented
below.  The lack of an optional feature in a given package is not
necessarily a bug.  More recommendations for GNU packages can be found
in *note Makefile Conventions: (standards)Makefile Conventions.

   The `configure` shell script attempts to guess correct values for
various system-dependent variables used during compilation.  It uses
those values to create a `Makefile` in each directory of the package.
It may also create one or more `.h` files containing system-dependent
definitions.  Finally, it creates a shell script `config.status` that
you can run in the future to recreate the current configuration, and a
file `config.log` containing compiler output (useful mainly for
debugging `configure`).

   It can also use an optional file (typically called `config.cache`
and enabled with `--cache-file=config.cache` or simply `-C`) that saves
the results of its tests to speed up reconfiguring.  Caching is
disabled by default to prevent problems with accidental use of stale
cache files.

   If you need to do unusual things to compile the package, please try
to figure out how `configure` could check whether to do them, and mail
diffs or instructions to the address given in the `README` so they can
be considered for the next release.  If you are using the cache, and at
some point `config.cache` contains results you don`t want to keep, you
may remove or edit it.

   The file `configure.ac` (or `configure.in`) is used to create
`configure` by a program called `autoconf`.  You need `configure.ac` if
you want to change it or regenerate `configure` using a newer version
of `autoconf`.

   The simplest way to compile this package is:

  1. `cd` to the directory containing the package's source code and type
     `./configure` to configure the package for your system.

     Running `configure` might take a while.  While running, it prints
     some messages telling which features it is checking for.

  2. Type `make` to compile the package.

  3. Optionally, type `make check` to run any self-tests that come with
     the package, generally using the just-built uninstalled binaries.

  4. Type `make install` to install the programs and any data files and
     documentation.  When installing into a prefix owned by root, it is
     recommended that the package be configured and built as a regular
     user, and only the `make install` phase executed with root
     privileges.

  5. Optionally, type `make installcheck` to repeat any self-tests, but
     this time using the binaries in their final installed location.
     This target does not install anything.  Running this target as a
     regular user, particularly if the prior `make install` required
     root privileges, verifies that the installation completed
     correctly.

  6. You can remove the program binaries and object files from the
     source code directory by typing `make clean`.  To also remove the
     files that `configure` created (so you can compile the package for
     a different kind of computer), type `make distclean`.  There is
     also a `make maintainer-clean` target, but that is intended mainly
     for the package`s developers.  If you use it, you may have to get
     all sorts of other programs in order to regenerate files that came
     with the distribution.

  7. Often, you can also type `make uninstall` to remove the installed
     files again.  In practice, not all packages have tested that
     uninstallation works correctly, even though it is required by the
     GNU Coding Standards.

  8. Some packages, particularly those that use Automake, provide `make
     distcheck`, which can by used by developers to test that all other
     targets like `make install` and `make uninstall` work correctly.
     This target is generally not run by end users.

Compilers and Options
---------------------

   Some systems require unusual options for compilation or linking that
the `configure` script does not know about.  Run `./configure --help`
for details on some of the pertinent environment variables.

   You can give `configure` initial values for configuration parameters
by setting variables in the command line or in the environment.  Here
is an example:

     ./configure CC=c99 CFLAGS=-g LIBS=-lposix

   *Note Defining Variables::, for more details.

Compiling For Multiple Architectures
------------------------------------

   You can compile the package for more than one kind of computer at the
same time, by placing the object files for each architecture in their
own directory.  To do this, you can use GNU `make`.  `cd` to the
directory where you want the object files and executables to go and run
the `configure` script.  `configure` automatically checks for the
source code in the directory that `configure` is in and in `..`.  This
is known as a "VPATH" build.

   With a non-GNU `make`, it is safer to compile the package for one
architecture at a time in the source code directory.  After you have
installed the package for one architecture, use `make distclean` before
reconfiguring for another architecture.

   On MacOS X 10.5 and later systems, you can create libraries and
executables that work on multiple system types--known as "fat" or
"universal" binaries--by specifying multiple `-arch` options to the
compiler but only a single `-arch` option to the preprocessor.  Like
this:

     ./configure CC="gcc -arch i386 -arch x86_64 -arch ppc -arch ppc64" \
                 CXX="g++ -arch i386 -arch x86_64 -arch ppc -arch ppc64" \
                 CPP="gcc -E" CXXCPP="g++ -E"

   This is not guaranteed to produce working output in all cases, you
may have to build one architecture at a time and combine the results
using the `lipo` tool if you have problems.

Installation Names
------------------

   By default, `make install` installs the package`s commands under
`/usr/local/bin`, include files under `/usr/local/include`, etc.  You
can specify an installation prefix other than `/usr/local` by giving
`configure` the option `--prefix=PREFIX`, where PREFIX must be an
absolute file name.

   You can specify separate installation prefixes for
architecture-specific files and architecture-independent files.  If you
pass the option `--exec-prefix=PREFIX` to `configure`, the package uses
PREFIX as the prefix for installing programs and libraries.
Documentation and other data files still use the regular prefix.

   In addition, if you use an unusual directory layout you can give
options like `--bindir=DIR` to specify different values for particular
kinds of files.  Run `configure --help` for a list of the directories
you can set and what kinds of files go in them.  In general, the
default for these options is expressed in terms of `${prefix}`, so that
specifying just `--prefix` will affect all of the other directory
specifications that were not explicitly provided.

   The most portable way to affect installation locations is to pass the
correct locations to `configure`; however, many packages provide one or
both of the following shortcuts of passing variable assignments to the
`make install` command line to change installation locations without
having to reconfigure or recompile.

   The first method involves providing an override variable for each
affected directory.  For example, `make install
prefix=/alternate/directory` will choose an alternate location for all
directory configuration variables that were expressed in terms of
`${prefix}`.  Any directories that were specified during `configure`,
but not in terms of `${prefix}`, must each be overridden at install
time for the entire installation to be relocated.  The approach of
makefile variable overrides for each directory variable is required by
the GNU Coding Standards, and ideally causes no recompilation.
However, some platforms have known limitations with the semantics of
shared libraries that end up requiring recompilation when using this
method, particularly noticeable in packages that use GNU Libtool.

   The second method involves providing the `DESTDIR` variable.  For
example, `make install DESTDIR=/alternate/directory` will prepend
`/alternate/directory` before all installation names.  The approach of
`DESTDIR` overrides is not required by the GNU Coding Standards, and
does not work on platforms that have drive letters.  On the other hand,
it does better at avoiding recompilation issues, and works well even
when some directory options were not specified in terms of `${prefix}`
at `configure` time.

Optional Features
-----------------

   If the package supports it, you can cause programs to be installed
with an extra prefix or suffix on their names by giving `configure` the
option `--program-prefix=PREFIX` or `--program-suffix=SUFFIX`.

   Some packages pay attention to `--enable-FEATURE` options to
`configure`, where FEATURE indicates an optional part of the package.
They may also pay attention to `--with-PACKAGE` options, where PACKAGE
is something like `gnu-as` or `x` (for the X Window System).  The
`README` should mention any `--enable-` and `--with-` options that the
package recognizes.

   For packages that use the X Window System, `configure` can usually
find the X include and library files automatically, but if it doesn`t,
you can use the `configure` options `--x-includes=DIR` and
`--x-libraries=DIR` to specify their locations.

   Some packages offer the ability to configure how verbose the
execution of `make` will be.  For these packages, running `./configure
--enable-silent-rules` sets the default to minimal output, which can be
overridden with `make V=1`; while running `./configure
--disable-silent-rules` sets the default to verbose, which can be
overridden with `make V=0`.

Particular systems
------------------

   On HP-UX, the default C compiler is not ANSI C compatible.  If GNU
CC is not installed, it is recommended to use the following options in
order to use an ANSI C compiler:

     ./configure CC="cc -Ae -D_XOPEN_SOURCE=500"

and if that doesn`t work, install pre-built binaries of GCC for HP-UX.

   HP-UX `make` updates targets which have the same time stamps as
their prerequisites, which makes it generally unusable when shipped
generated files such as `configure` are involved.  Use GNU `make`
instead.

   On OSF/1 a.k.a. Tru64, some versions of the default C compiler cannot
parse its `<wchar.h>` header file.  The option `-nodtk` can be used as
a workaround.  If GNU CC is not installed, it is therefore recommended
to try

     ./configure CC="cc"

and if that doesn`t work, try

     ./configure CC="cc -nodtk"

   On Solaris, don`t put `/usr/ucb` early in your `PATH`.  This
directory contains several dysfunctional programs; working variants of
these programs are available in `/usr/bin`.  So, if you need `/usr/ucb`
in your `PATH`, put it _after_ `/usr/bin`.

   On Haiku, software installed for all users goes in `/boot/common`,
not `/usr/local`.  It is recommended to use the following options:

     ./configure --prefix=/boot/common

Specifying the System Type
--------------------------

   There may be some features `configure` cannot figure out
automatically, but needs to determine by the type of machine the package
will run on.  Usually, assuming the package is built to be run on the
_same_ architectures, `configure` can figure that out, but if it prints
a message saying it cannot guess the machine type, give it the
`--build=TYPE` option.  TYPE can either be a short name for the system
type, such as `sun4`, or a canonical name which has the form:

     CPU-COMPANY-SYSTEM

where SYSTEM can have one of these forms:

     OS
     KERNEL-OS

   See the file `config.sub` for the possible values of each field.  If
`config.sub` isn`t included in this package, then this package doesn't
need to know the machine type.

   If you are _building_ compiler tools for cross-compiling, you should
use the option `--target=TYPE` to select the type of system they will
produce code for.

   If you want to _use_ a cross compiler, that generates code for a
platform different from the build platform, you should specify the
"host" platform (i.e., that on which the generated programs will
eventually be run) with `--host=TYPE`.

Sharing Defaults
----------------

   If you want to set default values for `configure` scripts to share,
you can create a site shell script called `config.site` that gives
default values for variables like `CC`, `cache_file`, and `prefix`.
`configure` looks for `PREFIX/share/config.site` if it exists, then
`PREFIX/etc/config.site` if it exists.  Or, you can set the
`CONFIG_SITE` environment variable to the location of the site script.
A warning: not all `configure` scripts look for a site script.

Defining Variables
----------------

   Variables not defined in a site shell script can be set in the
environment passed to `configure`.  However, some packages may run
configure again during the build, and the customized values of these
variables may be lost.  In order to avoid this problem, you should set
them in the `configure` command line, using `VAR=value`.  For example:

     ./configure CC=/usr/local2/bin/gcc

causes the specified `gcc` to be used as the C compiler (unless it is
overridden in the site shell script).

Unfortunately, this technique does not work for `CONFIG_SHELL` due to
an Autoconf limitation.  Until the limitation is lifted, you can use
this workaround:

     CONFIG_SHELL=/bin/bash ./configure CONFIG_SHELL=/bin/bash

`configure` Invocation
----------------

   `configure` recognizes the following options to control how it
operates.

`--help`
`-h`
     Print a summary of all of the options to `configure`, and exit.

`--help=short`
`--help=recursive`
     Print a summary of the options unique to this package`s
     `configure`, and exit.  The `short` variant lists options used
     only in the top level, while the `recursive` variant lists options
     also present in any nested packages.

`--version`
`-V`
     Print the version of Autoconf used to generate the `configure`
     script, and exit.

`--cache-file=FILE`
     Enable the cache: use and save the results of the tests in FILE,
     traditionally `config.cache`.  FILE defaults to `/dev/null` to
     disable caching.

`--config-cache`
`-C`
     Alias for `--cache-file=config.cache`.

`--quiet`
`--silent`
`-q`
     Do not print messages saying which checks are being made.  To
     suppress all normal output, redirect it to `/dev/null` (any error
     messages will still be shown).

`--srcdir=DIR`
     Look for the package`s source code in directory DIR.  Usually
     `configure` can determine that directory automatically.

`--prefix=DIR`
     Use DIR as the installation prefix.  *note Installation Names::
     for more details, including other options available for fine-tuning
     the installation locations.

`--no-create`
`-n`
     Run the configure checks, but stop before creating any output
     files.

`configure` also accepts some other, not widely useful, options.  Run
`configure --help` for more details.

