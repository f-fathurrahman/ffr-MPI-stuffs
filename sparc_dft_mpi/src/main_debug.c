#include <mpi.h>
#include <stdio.h>

#include "isddft.h"

void my_Initialize(SPARC_OBJ *pSPARC, int argc, char *argv[]);

int main(int argc, char* argv[])
{
  // set up MPI
  MPI_Init(&argc, &argv);
  
  // get communicator size and my rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
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