PROGRAM ex_gather

  IMPLICIT NONE 
  INCLUDE 'mpif.h'

  INTEGER :: ierr
  INTEGER :: Nprocs, my_rank
  INTEGER :: irecv(3), isend

  ! Nprocs should be equal to 3 because irecv is 3

  CALL MPI_Init( ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, Nprocs, ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )
 
  ! isend = my_rank + 1

  IF( my_rank == 1 ) THEN 
    isend = 77
  ELSEIF( my_rank == 2 ) THEN 
    isend = 4
  ELSE
    isend = 99
  ENDIF 

  CALL MPI_Gather( isend, 1, MPI_INTEGER, &
                   irecv, 1, MPI_INTEGER, &
                   0, MPI_COMM_WORLD, ierr )

  IF( my_rank == 0 ) THEN 
    WRITE(*,*) 'irecv = ', irecv
  ENDIF 

  CALL MPI_Finalize( ierr )

END PROGRAM 
