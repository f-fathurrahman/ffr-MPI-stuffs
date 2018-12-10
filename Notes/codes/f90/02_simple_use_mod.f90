PROGRAM simple
  USE mpi
  IMPLICIT NONE
  INTEGER :: rank, Nprocs, ierr
  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, Nprocs, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  WRITE(*,*) 'Hello, I am proc: ', rank, ' from ', Nprocs, ' total processes'
  CALL MPI_Finalize(ierr)
END PROGRAM

