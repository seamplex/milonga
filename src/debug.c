/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga debugging and benchmarking routines
 *
 *  Copyright (C) 2010--2015 jeremy theler
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
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <sys/resource.h>



#include <gsl/gsl_errno.h>

#include <petsc.h>
#include <slepcsys.h>
#include <petscdraw.h>

#include "milonga.h"
#include "mat2sng.h"
#include "version.h"

#define milonga_debug_insert_spaces(n)     for (j = 0; j < (n); j++) petsc_call(PetscViewerASCIIPrintf(debug->viewer, " "));
#define wasora_null_free(p) if (p!=NULL) {free(p); p=NULL;}

#undef  __FUNCT__
#define __FUNCT__ "milonga_debug_n_processors"
int milonga_debug_n_processors(void) {

  FILE *proccpuinfo;
  int n_processors = 0;
  char *buffer;

  assert(proccpuinfo = fopen("/proc/cpuinfo", "r"));
  buffer = malloc(BUFFER_SIZE);

  while (fscanf(proccpuinfo, "%s", buffer) != EOF) {
    if (strcmp("processor", buffer) == 0) {
      n_processors++;
    }
  }

  free(buffer);
  fclose(proccpuinfo);

  return n_processors;

}

#undef  __FUNCT__
#define __FUNCT__ "milonga_debug_cpu_info"
int milonga_debug_cpu_info(void) {

  FILE *proccpuinfo;
  int n;
  int field_width;

  char *line;
  char *value;
  char *dummy;


  line = malloc(BUFFER_SIZE);   // no son tan largas como para usar MAX_LINE
  value = malloc(BUFFER_SIZE);

  cpuinfo.n_processors = milonga_debug_n_processors();

  cpuinfo.model_name = calloc(cpuinfo.n_processors, sizeof(char *));
  cpuinfo.mhz = calloc(cpuinfo.n_processors, sizeof(char *));
  cpuinfo.bogomips = calloc(cpuinfo.n_processors, sizeof(char *));
  cpuinfo.cache_size = calloc(cpuinfo.n_processors, sizeof(char *));

  if ((proccpuinfo = fopen("/proc/cpuinfo", "r")) == NULL) {
    wasora_push_error_message("cannot open '/proc/cpuinfo/'");
    return WASORA_RUNTIME_ERROR;
  }

  n = -1;
  while (fgets(line, BUFFER_SIZE-1, proccpuinfo)) {

    if (strncmp(line, "processor", 9) == 0) {
      n++;
    }

    if ((dummy = strchr(line, ':')) != NULL) {
      field_width = dummy-line;

      strcpy(value, line+field_width+1);

      if (strncmp(line, "model name", 10) == 0) {
        cpuinfo.model_name[n] = strdup(value);
        cpuinfo.model_name[n][strlen(cpuinfo.model_name[n])-1] = '\0';
      } else if (strncmp(line, "cpu MHz", 7) == 0) {
        cpuinfo.mhz[n] = strdup(value);
        cpuinfo.mhz[n][strlen(cpuinfo.mhz[n])-1] = '\0';
      } else if (strncmp(line, "cache size", 10) == 0) {
        cpuinfo.cache_size[n] = strdup(value);
        cpuinfo.cache_size[n][strlen(cpuinfo.cache_size[n])-1] = '\0';
      } else if (strncmp(line, "bogomips", 8) == 0) {
        cpuinfo.bogomips[n] = strdup(value);
        cpuinfo.bogomips[n][strlen(cpuinfo.bogomips[n])-1] = '\0';
      }

    }
  }

  for (n = 0; n < cpuinfo.n_processors; n++) {
    if (cpuinfo.model_name[n] == NULL) {
      cpuinfo.model_name[n] = strdup("unknown");
    }
    if (cpuinfo.mhz[n] == NULL) {
      cpuinfo.mhz[n] = strdup("unknown");
    }
    if (cpuinfo.cache_size[n] == NULL) {
      cpuinfo.cache_size[n] = strdup("unknown");
    }
    if (cpuinfo.bogomips[n] == NULL) {
      cpuinfo.bogomips[n] = strdup("unknown");
    }
  }

  free(line);
  free(value);
  fclose(proccpuinfo);

  return WASORA_RUNTIME_OK;


}

#undef  __FUNCT__
#define __FUNCT__ "milonga_debug_cpu_info_free"
int milonga_debug_cpu_info_free(void) {
  int n;
  
  if (cpuinfo.model_name != NULL) {
    for (n = 0; n < cpuinfo.n_processors; n++) {
      wasora_null_free(cpuinfo.model_name[n]);
      wasora_null_free(cpuinfo.mhz[n]);
      wasora_null_free(cpuinfo.bogomips[n]);
      wasora_null_free(cpuinfo.cache_size[n]);
    }
    wasora_null_free(cpuinfo.model_name);
    wasora_null_free(cpuinfo.mhz);
    wasora_null_free(cpuinfo.bogomips);
    wasora_null_free(cpuinfo.cache_size);
  }
  
  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "milonga_debug_open"
int milonga_debug_open(struct debug_t *debug) {
  int i, j;
  time_t tm;
  struct utsname computer;
  char *newfilepath;
  char *markdownpath;

  char libversion[BUFFER_SIZE];

  uname(&computer);
  time(&tm);

  if (debug->file == NULL) {
    wasora_push_error_message("no FILE given to MILONGA_DEBUG");
    return WASORA_RUNTIME_ERROR;
  }

  // hacemos como en openfile pero petsc necesita su call
  if ((newfilepath = wasora_evaluate_string(debug->file->format, debug->file->n_args, debug->file->arg)) == NULL) {
    return WASORA_RUNTIME_ERROR;
  }
  debug->file->path = strdup(newfilepath);
  markdownpath = malloc(strlen(newfilepath) + 8);
  snprintf(markdownpath, strlen(newfilepath) + 6, "%s.md", newfilepath);
  
  debug->file_opened = 1;
  petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, markdownpath, &debug->viewer));
  
  free(newfilepath);
  free(markdownpath);

  

  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%% milonga debugging and benchmarking output\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%% %s\n", getlogin()));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%% %s\n", ctime(&tm)));
  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "milonga debugging and benchmarking output\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "=========================================\n\n"));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "code invocation\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "---------------\n\n"));
   
// PetscGetUserName()
// PetscGetHostName()  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "    %s@%s$", getlogin(), computer.nodename));
  for (i = 0; i < wasora.argc; i++) {
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, " %s", wasora.argv[i]));
  }
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\non %s\n\n", ctime(&tm)));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "code version\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "------------\n\n"));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "[%s](http://www.talador.com.ar/jeremy/wasora/milonga) %s  \n", plugin_name(), plugin_version()));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s  \n\n", plugin_description()));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s\n", plugin_longversion()));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s\n", plugin_copyright()));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n"));
  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "milonga was compiled on %s  \n", COMPILATION_DATE));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "by %s on %s (%s) using %s  \n", COMPILATION_USERNAME, COMPILATION_HOSTNAME, COMPILATION_ARCH, CCOMPILER_VERSION));
  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "and linked against  \n"));
  petsc_call(PetscGetVersion(libversion, BUFFER_SIZE));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s  \n", libversion));
  petsc_call(SlepcGetVersion(libversion, BUFFER_SIZE));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s  \n", libversion));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n"));


  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));
  
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "computer information\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "--------------------\n\n"));

  milonga_debug_cpu_info();

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "------------                     ------------\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "architecture                     %s\n", computer.machine));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "CPU model name                   %s\n", cpuinfo.model_name[0]));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of processors present     %d\n", cpuinfo.n_processors));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of online processors      %ld\n", sysconf(_SC_NPROCESSORS_ONLN)));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "total physical memory            %ld Gb\n", sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGESIZE)/(1024*1024*1024)));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "page size                        %ld bytes\n", sysconf(_SC_PAGESIZE)));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "system kernel name               %s\n", computer.sysname));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "kernel version                   %s\n", computer.version));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "kernel release                   %s\n", computer.release));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "------------                     ------------\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));
  

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "------     -----      -----            -----     \n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "processor  bogomips   MHz              cache\n"));
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "------     -----      -----            -----     \n"));
  for (i = 0; i < cpuinfo.n_processors; i++) {
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%d", i));
    milonga_debug_insert_spaces(9);
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s", cpuinfo.bogomips[i]));
    milonga_debug_insert_spaces(11-strlen(cpuinfo.bogomips[i]));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s", cpuinfo.mhz[i]));
    milonga_debug_insert_spaces(17-strlen(cpuinfo.mhz[i]));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "%s\n\n", cpuinfo.cache_size[i]));
  }
  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "------     -----      -----            -----     \n"));

  petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));
   


  return WASORA_RUNTIME_OK;
}


#undef  __FUNCT__
#define __FUNCT__ "milonga_instruction_debug"
int milonga_instruction_debug(void *arg) {
  
  debug_t *debug = (debug_t *)arg;

  struct rusage resource_usage;

  PetscInt size;
  PetscViewer viewer, viewer2;
  char *filename;

  wasora_value(milonga.vars.available_memory) = sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGESIZE);
  getrusage(RUSAGE_SELF, &resource_usage);
  wasora_value(milonga.vars.memory_usage_global) = (double)(1024.0*resource_usage.ru_maxrss);
  petsc_call(PetscMemoryGetMaximumUsage(wasora_value_ptr(milonga.vars.memory_usage_petsc)));
  petsc_call(PetscGetFlops(wasora_value_ptr(milonga.vars.flops_petsc)));
  
  if (debug->matrices_size.n_tokens != 0) {
    size = (int)(wasora_evaluate_expression(&debug->matrices_size));
  } else {
    size = DEFAULT_MATRICES_SIZE;
  }
  

  
  if (debug->file != NULL) {
    if (debug->file_opened == 0) {
      wasora_call(milonga_debug_open(debug));
    }
    filename = malloc(strlen(debug->file->path)+32);

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "problem static step %d\n", (int)(wasora_var(wasora_special_var(step_static)))));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-----------------------\n", (int)(wasora_var(wasora_special_var(step_static)))));


    // eigencurrent solution
    if (milonga.eps != NULL) {
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n### SLEPc's EPSView output\n\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
      petsc_call(EPSView(milonga.eps, debug->viewer));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n\n"));
        
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "### Eigenvalue problem result\n\n"));

      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "---------------------- ---------------                                              ------------------\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "requested eigenvalue   $\\lambda$                                                    %.10f\n", wasora_value(milonga.vars.keff)));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "residual norm          $\\|R \\phi - \\lambda F \\phi\\|_2$                              %.4g\n", wasora_value(milonga.vars.residual_norm)));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "relative error         $\\|R \\phi - \\lambda F \\phi\\|_2 / \\| \\lambda F \\phi \\|_2$     %.4g\n", wasora_value(milonga.vars.rel_error)));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "error estimate         $\\|\\lambda - \\lambda_\\text{real}\\|$                          %.4g\n", wasora_value(milonga.vars.error_estimate)));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "---------------------- ---------------                                              ------------------\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));
        
    } else if (milonga.ksp != NULL) {
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n### PETSc's KSPView output\n\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
      petsc_call(KSPView(milonga.ksp, debug->viewer));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n\n"));
    }

    // system resources
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "### System resource usage\n\n"));

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "user CPU time                          %.4g seconds\n", resource_usage.ru_utime.tv_sec + 1e-6*resource_usage.ru_utime.tv_usec));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "system CPU time                        %.4g seconds\n", resource_usage.ru_stime.tv_sec + 1e-6*resource_usage.ru_stime.tv_usec));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "total                                  %.4g seconds\n", resource_usage.ru_utime.tv_sec + resource_usage.ru_stime.tv_sec + 1e-6*(resource_usage.ru_utime.tv_usec + resource_usage.ru_stime.tv_usec)));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------ ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "time insumed to                        cpu [secs]         wall [secs]\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------ ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "build the matrices                     %.4g               %.4g\n", wasora_value(milonga.vars.time_cpu_build), wasora_value(milonga.vars.time_wall_build)));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "solve the eigencurrent problem         %.4g               %.4g\n", wasora_value(milonga.vars.time_cpu_solve), wasora_value(milonga.vars.time_wall_solve)));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "total                                  %.4g               %.4g\n", wasora_value(milonga.vars.time_cpu_build)+wasora_value(milonga.vars.time_cpu_solve), wasora_value(milonga.vars.time_wall_build)+wasora_value(milonga.vars.time_wall_solve)));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------ ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "maximum memory resident set size       %ld kb\n", resource_usage.ru_maxrss));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "PETSc's maximum memory set size        %.0f kb\n", wasora_value(milonga.vars.memory_usage_petsc)/1024.0));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of soft page faults             %ld\n", resource_usage.ru_minflt));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of hard page faults             %ld\n", resource_usage.ru_majflt));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of swaps                        %ld\n", resource_usage.ru_majflt));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of block input operations       %ld\n", resource_usage.ru_inblock));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "number of block output operations      %ld\n", resource_usage.ru_oublock));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "-------------------------------------- ------------------\n"));

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n"));
    

    if (debug->matrices & DEBUG_MATRICES_ASCII) {
      PetscViewer ascii_file;

      sprintf(filename, "%s-R.txt", debug->file->path);
      petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
      milonga_print_petsc_matrix(milonga.R, ascii_file);
      petsc_call(PetscViewerDestroy(&ascii_file));
      
      if (milonga.has_fission && !milonga.has_sources) {
        sprintf(filename, "%s-F.txt", debug->file->path);
        petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
        milonga_print_petsc_matrix(milonga.F, ascii_file);
        petsc_call(PetscViewerDestroy(&ascii_file));
      } else if (milonga.has_sources) {
        sprintf(filename, "%s-S.txt", debug->file->path);
        petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
        milonga_print_petsc_vector(milonga.S, ascii_file);
        petsc_call(PetscViewerDestroy(&ascii_file));
      }

    }

    if (debug->matrices & DEBUG_MATRICES_ASCII_STRUCT) {
      PetscViewer ascii_file;

      sprintf(filename, "%s-R.str", debug->file->path);
      petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
      milonga_print_petsc_matrix_struct(milonga.R, ascii_file);
      petsc_call(PetscViewerDestroy(&ascii_file));

      if (milonga.has_fission && !milonga.has_sources) {
        sprintf(filename, "%s-F.str", debug->file->path);
        petsc_call(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &ascii_file));
        milonga_print_petsc_matrix_struct(milonga.F, ascii_file);
        petsc_call(PetscViewerDestroy(&ascii_file));
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_BINARY) {
      sprintf(filename, "%s-R.bin", debug->file->path);
      petsc_call(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer));
      petsc_call(MatView(milonga.R, viewer));
      petsc_call(PetscViewerDestroy(&viewer));
      
      if (milonga.has_fission && !milonga.has_sources) {
        sprintf(filename, "%s-F.bin", debug->file->path);
        petsc_call(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer));
        MatView(milonga.F, viewer);
        petsc_call(PetscViewerDestroy(&viewer));
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_COMPRESSED_BINARY) {
      sprintf(filename, "%s-R.gz", debug->file->path);
      petsc_call(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer));
      MatView(milonga.R, viewer);
      petsc_call(PetscViewerDestroy(&viewer));
      
      if (milonga.has_fission && !milonga.has_sources) {
        sprintf(filename, "%s-F.gz", debug->file->path);
        petsc_call(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer));
        MatView(milonga.F, viewer);
        petsc_call(PetscViewerDestroy(&viewer));
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_ASCII) {
      sprintf(filename, "%s-R.asc", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT);
      MatView(milonga.R, viewer);
      petsc_call(PetscViewerDestroy(&viewer));

      if (milonga.has_fission && !milonga.has_sources) {
        sprintf(filename, "%s-F.asc", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT);
        MatView(milonga.F, viewer);
        petsc_call(PetscViewerDestroy(&viewer));
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_OCTAVE) {
      sprintf(filename, "%s-R.m", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
      MatView(milonga.R, viewer);
      petsc_call(PetscViewerDestroy(&viewer));

      if (milonga.has_fission && !milonga.has_sources) {
        sprintf(filename, "%s-F.m", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
        MatView(milonga.F, viewer);
        petsc_call(PetscViewerDestroy(&viewer));
      }
    }

    if (debug->matrices & DEBUG_MATRICES_PETSC_DENSE) {
      sprintf(filename, "%s-R.den", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
      MatView(milonga.R, viewer);
      petsc_call(PetscViewerDestroy(&viewer));

      if (milonga.has_fission) {
        sprintf(filename, "%s-F.den", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
        MatView(milonga.F, viewer);
        petsc_call(PetscViewerDestroy(&viewer));
      }
    }
    
    if (debug->matrices & DEBUG_MATRICES_SNG) {
      sprintf(filename, "%s-R.sng", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT);
      mat2sng(milonga.R, size, (PetscInt)(wasora_evaluate_expression(&debug->matrices_stride)), 0, viewer);
      petsc_call(PetscViewerDestroy(&viewer));

      if (milonga.has_fission && !milonga.has_sources) {
        sprintf(filename, "%s-F.sng", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT);
        mat2sng(milonga.F, size, (PetscInt)(wasora_evaluate_expression(&debug->matrices_stride)), 0, viewer);
        petsc_call(PetscViewerDestroy(&viewer));
      }
    }    

    if (debug->matrices & DEBUG_MATRICES_SNG_STRUCT) {
      sprintf(filename, "%s-str-R.sng", debug->file->path);
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
      PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT);
      mat2sng(milonga.R, size, (PetscInt)(wasora_evaluate_expression(&debug->matrices_stride)), 1, viewer);
      petsc_call(PetscViewerDestroy(&viewer));

      if (milonga.has_fission && !milonga.has_sources) {
        sprintf(filename, "%s-str-F.sng", debug->file->path);
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
        PetscViewerPushFormat(viewer, PETSC_VIEWER_DEFAULT);
        mat2sng(milonga.F, size, (PetscInt)(wasora_evaluate_expression(&debug->matrices_stride)), 1, viewer);
        petsc_call(PetscViewerDestroy(&viewer));
      }
    }    
    
    free(filename);

  }

  if (debug->matrices & DEBUG_MATRICES_X) {

    PetscDraw draw;
    PetscDraw draw2;

    PetscViewerDrawOpen(PETSC_COMM_WORLD, PETSC_NULL, "R", 100, 100, size, size, &viewer);
    MatView(milonga.R, viewer);
    PetscViewerDrawGetDraw(viewer, 0, &draw);
    PetscDrawSetPause(draw, -1);

    if (milonga.has_fission && !milonga.has_sources) {
      PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL, "F", size+100, size+100, size, size, &viewer2);
      MatView(milonga.F, viewer2);
      PetscViewerDrawGetDraw(viewer2, 0, &draw2);
      PetscDrawSetPause(draw2, -1);
//      PetscDrawPause(draw2);
    }
    
    PetscDrawPause(draw);

    petsc_call(PetscViewerDestroy(&viewer));
    if (milonga.has_fission && !milonga.has_sources) {
      petsc_call(PetscViewerDestroy(&viewer2));
    }
  }
  
  if ((int)(wasora_var(wasora_special_var(static_steps))) == 1 || (int)(wasora_var(wasora_special_var(done))) == 1) {
    milonga_debug_close(debug);
  }

  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "milonga_debug_close"
int milonga_debug_close(debug_t *debug) {

  int c;
  FILE *finput;

  if (debug->file != NULL) {
  
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "PETSc's LogView output\n"));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "----------------------\n"));
  
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
    petsc_call(PetscLogView(debug->viewer));
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n\n"));
  
    if (debug->include_input) {
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "transcription of input file\n"));
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "---------------------------\n\n"));

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n"));
    finput = fopen(wasora.argv[1], "r");
      while ((c = fgetc(finput)) != EOF) {
        fputc(c, debug->file->pointer);
      }
      fclose(finput);
      petsc_call(PetscViewerASCIIPrintf(debug->viewer, "~~~~\n\n\n"));
    }

    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\n\n\
*  *  *  *\n\n\
                           .,,.\n\
                       ,dKWMMMMW0d:.\n\
                     .XMMMMMMMMMMMMMWXkdlllool:.\n\
                    .NMMMMMMMMMMMMMMMMMMMMMMMMMMX,\n\
                    dMMMMMMMMMMMMMMMMMMMMMMMMMMMMMc\n\
                    KMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM:\n\
                  .oMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM:\n\
        ...',;coxKMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMd.\n\
      xkkkkkkxxddoolllllllodxO0KKKKKKKKKKKKKXXXXXNNNNNNNXK0000000000000OOx;\n\n"));
  
  
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "[milonga"));
#ifdef PLUGIN_VCS_BRANCH
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, " %s", PLUGIN_VCS_VERSION));
#else
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, " %s", PACKAGE_VERSION));
#endif
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "](http://www.seamplex.com/milonga)\n\n"));

  /*
    petsc_call(PetscViewerASCIIPrintf(debug->viewer, "\
                       .''.\n\
                   ....    .....\n\
                 ..             ............\n\
                .'                          ..\n\
                ;                             '\n\
                ;                             .'\n\
               '                               .'\n\
    ..........                                   ..\n\
  ,................................................'''''''''''''''''''.\n\
                                                                     .\n\n"));
*/
    

//    PetscViewerDestroy(&debug->viewer);
  }

  milonga_debug_cpu_info_free();
  return WASORA_RUNTIME_OK;

}


#undef  __FUNCT__
#define __FUNCT__ "milonga_print_petsc_vector"
int milonga_print_petsc_vector(Vec b, PetscViewer viewer) {

  double xi;
  int i;
  int m;

  VecGetSize(b, &m);

  for (i = 0; i < m; i++) {
    VecGetValues(b, 1, &i, &xi);
    if (xi != 0) {
      petsc_call(PetscViewerASCIIPrintf(viewer, "% .1e ", xi));
    } else {
      petsc_call(PetscViewerASCIIPrintf(viewer, "    0    "));
    }
    petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
  }
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "milonga_print_petsc_matrix"
int milonga_print_petsc_matrix(Mat A, PetscViewer viewer) {

  double xi;
  int i, j;
  int m, n;

  MatGetSize(A, &m, &n);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      MatGetValues(A, 1, &i, 1, &j, &xi);
      if (xi != 0) {
        petsc_call(PetscViewerASCIIPrintf(viewer, "% .1e ", xi));
      } else {
        petsc_call(PetscViewerASCIIPrintf(viewer, "    0    "));
      }
    }
    petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
  }
  
  return WASORA_RUNTIME_OK;

}

#undef  __FUNCT__
#define __FUNCT__ "milonga_print_petsc_matrix_struct"
int milonga_print_petsc_matrix_struct(Mat A, PetscViewer viewer) {

  double xi;
  int i, j;
  int m, n;

  MatGetSize(A, &m, &n);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      MatGetValues(A, 1, &i, 1, &j, &xi);
      if (xi != 0) {
         petsc_call(PetscViewerASCIIPrintf(viewer, "#"));
      } else {
        petsc_call(PetscViewerASCIIPrintf(viewer, " "));
      }
    }
    petsc_call(PetscViewerASCIIPrintf(viewer, "\n"));
  }

  return WASORA_RUNTIME_OK;
  
}
