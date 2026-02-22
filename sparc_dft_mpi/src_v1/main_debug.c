
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <mpi.h>

#include "isddft.h"

// ... declare debugging functions here ...
//void my_Initialize(SPARC_OBJ *pSPARC, int argc, char *argv[]);

void investigate_sparc_obj(SPARC_OBJ *pSPARC) {
  printf("Nx_d = %d\n", pSPARC->Nx_d);
  return;
}


int main(int argc, char* argv[])
{  
  // get communicator size and my rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  int print_restart_typ = 0;
  double t_init, t_acc, *avgvel, *maxvel, *mindis;

  // set up MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // This is our central object
  SPARC_OBJ SPARC;
  
  double t1, t2;
  // why this is needed ?
  MPI_Barrier(MPI_COMM_WORLD);
  
  // start timer
  t1 = MPI_Wtime();
  SPARC.time_start = t1;

  // Read files and initialize
  my_Initialize(&SPARC, argc, argv);

  investigate_sparc_obj(&SPARC);

  t_init = MPI_Wtime();
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  avgvel = (double *)malloc(SPARC.Ntypes * sizeof(double) );
  maxvel = (double *)malloc(SPARC.Ntypes * sizeof(double) );
  mindis = (double *)malloc(SPARC.Ntypes*(SPARC.Ntypes+1)/2 * sizeof(double) );

  // Check whether the restart has to be performed
  if(SPARC.RestartFlag != 0){
    // Check if .restart file present
    if(rank == 0){
        FILE *rst_fp = NULL;
        if( access(SPARC.restart_Filename, F_OK ) != -1 )
        rst_fp = fopen(SPARC.restart_Filename,"r");
      else if( access(SPARC.restartC_Filename, F_OK ) != -1 )
        rst_fp = fopen(SPARC.restartC_Filename,"r");
      else
        rst_fp = fopen(SPARC.restartP_Filename,"r");
      
        if(rst_fp == NULL)
            SPARC.RestartFlag = 0;
      }
      MPI_Bcast(&(SPARC.RestartFlag), 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  // Initialize the MD for the very first step only
  Calculate_electronicGroundState(&SPARC);

  MPI_Barrier(MPI_COMM_WORLD);
  // end timer
  t2 = MPI_Wtime();
  if (rank == 0) {
    printf("The program took %.3f s.\n", t2 - t1);
  }

  // ensure stdout flushed to prevent Finalize hang
  fflush(stdout);

  // finalize MPI
  MPI_Finalize();
  return 0;

}