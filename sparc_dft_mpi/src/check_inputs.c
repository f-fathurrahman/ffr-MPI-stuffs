#include <stdio.h>
#include <stdlib.h>

#include "initialization.h"

/**
 * @brief   Check input arguments and read filename.
 */
void check_inputs(SPARC_INPUT_OBJ *pSPARC_Input, int argc, char *argv[]) {
#ifdef DEBUG
  printf("Checking input arguments parsed by command line.\n");
#endif
  int i;
  char *pch, name_flag = 'N';

  // extract SPARCROOT from the executable's name
  char *libbase = argv[0];
  pch = strrchr(libbase, '/'); // find last occurrence of '/'
  // printf (RED "Last occurence of '/' found at %ld \n" RESET,pch-libbase+1);
  if (pch == NULL) {
    // in fact, unless SPARCROOT/lib is added to $(PATH), ./sparc is needed
    // strcpy(pSPARC_Input->SPARCROOT,".."); // in case '/' is not found
    snprintf(pSPARC_Input->SPARCROOT, L_STRING, "..");
  } else {
    memcpy(pSPARC_Input->SPARCROOT, libbase, pch - libbase);
    pSPARC_Input->SPARCROOT[pch - libbase] = '\0';
    if (strcmp(pSPARC_Input->SPARCROOT + (int)(pch - libbase) - 4, "/lib") == 0) {
      pSPARC_Input->SPARCROOT[pch - libbase - 4] = '\0'; // directly truncate string
    } else {
      strcat(pSPARC_Input->SPARCROOT, "/..");
    }
  }

  // save input filename
  memset(pSPARC_Input->filename, '\0', sizeof(pSPARC_Input->filename));
  pSPARC_Input->num_node = 0;
  pSPARC_Input->num_cpu_per_node = 0;
  pSPARC_Input->num_acc_per_node = 0;
  for (i = 1; i < argc - 1; i++) {
    if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      print_usage();
      exit(EXIT_FAILURE);
    }
    if (strcmp(argv[i], "-name") == 0) {
      name_flag = 'Y';
      simplifyPath(argv[i + 1], pSPARC_Input->filename, L_STRING);
      // break;
    }
    if (strcmp(argv[i], "-n") == 0) {
      pSPARC_Input->num_node = atoi(argv[i + 1]);
    }
    if (strcmp(argv[i], "-c") == 0) {
      pSPARC_Input->num_cpu_per_node = atoi(argv[i + 1]);
    }
    if (strcmp(argv[i], "-a") == 0) {
      pSPARC_Input->num_acc_per_node = atoi(argv[i + 1]);
    }
  }

  if (name_flag != 'Y') {
    print_usage();
    exit(EXIT_FAILURE);
  }
}
