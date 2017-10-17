! ffr, June 2016

PROGRAM hello

  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: my_rank
  INTEGER :: ierr
  INTEGER :: nprocs

  CALL MPI_Init( ierr )

  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )

  IF( my_rank == 0 ) THEN
    WRITE(*,*) 'Hello world MPI'
    WRITE(*,*) 'Number of procs: ', nprocs
  ENDIF

  CALL MPI_Finalize( ierr )

END PROGRAM

