#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  int ierr, Nprocs, my_rank;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  printf("Hello, I am proc: %d from %d processes\n", my_rank, Nprocs);
  ierr = MPI_Finalize();
  
  return 0;
}

