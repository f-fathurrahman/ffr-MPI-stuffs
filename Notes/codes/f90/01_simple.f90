PROGRAM simple
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: my_rank, Nprocs, ierr
  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, Nprocs, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  WRITE(*,*) 'Hello, I am proc: ', my_rank, ' from ', Nprocs, ' total processes'
  CALL MPI_Finalize(ierr)
END PROGRAM

