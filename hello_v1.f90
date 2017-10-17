! ffr, June 2016

PROGRAM hello

  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: my_rank
  INTEGER :: ierr
  INTEGER :: nprocs

  CALL MPI_Init( ierr )

  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )

  WRITE(*,*) 'Hello world MPI from proc:', my_rank

  CALL MPI_Finalize( ierr )

END PROGRAM

