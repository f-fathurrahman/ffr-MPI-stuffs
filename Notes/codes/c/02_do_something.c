#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  int ierr, my_rank, Nprocs;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if( my_rank == 0) {
    printf("Doing something for my_rank == 0\n");
  }
  else if( my_rank == 1 ) {
    printf("Doing something for my_rank == 1\n");
  }
  else {
    printf("Doing something for my_rank == %d\n", my_rank);
  }

  ierr = MPI_Finalize();

  printf("This is done outside of MPI\n");

  return 0;
}

