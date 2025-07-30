#include <stdio.h>
#include "initialization.h"

/**
 * @brief   Prints usage of SPARC through command line.
 */
void print_usage() {
  printf("\n");
  printf("USAGE:\n");
  printf("    mpirun -np <nproc> {SPARCROOT}/lib/sparc -name <filename>\n");
  printf("\n");
  printf("    {SPARCROOT} is the location of the SPARC folder\n");
  printf("\n");
  printf("REQUIRED ARGUMENT:\n");
  printf("    -name <filename>\n");
  printf("           The filename shared by .inpt file and .ion\n");
  printf("           file (without extension)\n");
  printf("\n");
  printf("OPTIONS: \n");
  printf("    -h, --help\n");
  printf("           Display help (from command line).\n");
  printf("    -n <number of Nodes>\n");
  printf("    -c <number of CPUs per node>\n");
  printf("    -a <number of Accelerators (e.g., GPUs) per node>\n");
  printf("\n");
  printf("EXAMPLE:\n");
  printf("\n");
  printf("    mpirun -np 8 {SPARCROOT}/lib/sparc -name test\n");
  printf("\n");
  printf("    The example command runs sparc with 8 cores, with input file "
         "named\n");
  printf("    test.inpt, and ion file named test.ion.\n");
  printf("\n");
  printf("NOTE: \n");
  printf("    This is a short description of the usage of SPARC. For a "
         "detailed \n");
  printf("    discription, refer to the manual online at\n");
  printf("\n");
  printf("        https://github.com/SPARC-X/SPARC/tree/master/doc \n");
  printf("\n");
}

