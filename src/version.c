/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's version banner
 *
 *  Copyright (C) 2010--2018 jeremy theler
 *
 *  This file is part of milonga.
 *
 *  milonga is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  milonga is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with wasora.  If not, see <http://www.gnu.org/licenses/>.
 *------------------- ------------  ----    --------  --     -       -         -
 */

#include <stdio.h>
#include <sys/utsname.h>

#include "milonga.h"
#include "version.h"

// global static so the compiler locates these strings in the text section
// so when the plugin_* functions return pointers to the strings there is
// no need to free them afterward
const char milonganame[] = "milonga";
const char milongadescription[] = "free nuclear reactor core analysis code";
char milongashortversion[128];
char milongalongversion[2048];
const char milongausage[] = "\
  --diffusion           use the diffusion formulation\n\
  --s2                  use the discrete ordinates S2 formulation\n\
  --s4                  use the discrete ordinates S4 formulation\n\
  --s6                  use the discrete ordinates S6 formulation\n\
  --s8                  use the discrete ordinates S8 formulation\n\
                        (they override the MILONGA_PROBLEM FORMULATION keyword))\n\
  --elements            use a finite-elements-based spatial discretization scheme\n\
  --volumes             use a finite-volumes-based spatial discretization scheme\n\
                        (both override the MILONGA_PROBLEM SCHEME keyword))\n\
  --largest             solve R phi = 1/k F phi (largest eigenvalue spectrum)\n\
  --smallest            solve F phi =  k  R phi (smallest eigenvalue spectrum)\n\
                        (both override the MILONGA_SOLVER SPECTRUM keyword))\n\
  --petsc_opt <option[=argument]>\n\
                        pass \"-option argument\" directly to PETSc/SLEPc, e.g.\n\
    $ milonga slab.was --petsc_opt malloc_dump --petsc_opt eps_type=gd\n\
\n\
  command line options override keywords given in the input, e.g.\n\
    $ milonga slab.was --elements --largest --slecp_opt st_type=sinvert\n\
  uses finite elements, largest eigenvalue formulation and shift & invert transform\n\
  even if the input contains\n\
    MILONGA_PROBLEM SCHEME volumes\n\
    MILONGA_SOLVER SPECTRUM smallest_eigenvalue ST_TYPE shift\n";

const char milongacopyright[] = "\
 milonga is copyright (c) 2010-2018 jeremy theler\n\
 licensed under GNU GPL version 3 or later.\n\
 milonga is free software: you are free to change and redistribute it.\n\
 There is NO WARRANTY, to the extent permitted by law.";

const char milongahmd5[] = PLUGIN_HEADERMD5;


const char *plugin_name(void) {
  return milonganame;
}

const char *plugin_longversion(void) {
#ifdef PLUGIN_VCS_BRANCH
  char slepcversion[BUFFER_SIZE];
  char petscversion[BUFFER_SIZE];
  char petscarch[BUFFER_SIZE];
  struct utsname computer;
  
  PetscGetVersion(petscversion, BUFFER_SIZE);
  SlepcGetVersion(slepcversion, BUFFER_SIZE);  
  PetscGetArchType(petscarch, BUFFER_SIZE);

  uname(&computer);
  if (milonga_debug_n_processors() != 0) {
    milonga_debug_cpu_info();
  }
  
  sprintf(milongalongversion, "\n\
 last commit on %s\n\
 compiled on %s by %s@%s (%s)\n\
 with %s using %s linked against\n\
  %s\n\
  %s %s\n\
 running on %s %s %s %s\n\
 %d %s\n",
   PLUGIN_VCS_DATE,
   COMPILATION_DATE,
   COMPILATION_USERNAME,
   COMPILATION_HOSTNAME,
   COMPILATION_ARCH,
   CCOMPILER_VERSION,
   CCOMPILER_FLAGS,
   slepcversion, petscversion, petscarch,
   computer.sysname, computer.release, computer.version, computer.machine,
   cpuinfo.n_processors,
   (cpuinfo.model_name != NULL)?cpuinfo.model_name[0]:"?");

   milonga_debug_cpu_info_free();

#endif  
  
  
  
  return milongalongversion;
}

const char *plugin_wasorahmd5(void) {
  return milongahmd5;
}
const char *plugin_copyright(void) {
  return milongacopyright;
}


const char *plugin_version(void) {
#ifdef PLUGIN_VCS_BRANCH
  sprintf(milongashortversion, "%s%s %s", PLUGIN_VCS_VERSION,
                                       (PLUGIN_VCS_CLEAN==0)?"":"+Δ",
                                       strcmp(PLUGIN_VCS_BRANCH, "master")?PLUGIN_VCS_BRANCH:"");
#else
  sprintf(milongashortversion, "%s", PACKAGE_VERSION);
#endif  
    
  return milongashortversion;
}

const char *plugin_description(void) {
  return milongadescription;
}

const char *plugin_usage(void) {
  return milongausage;
}
