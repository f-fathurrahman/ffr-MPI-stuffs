PROGRAM ex_bcast
  IMPLICIT NONE 
  INCLUDE 'mpif.h'

  INTEGER :: imsg(4)
  INTEGER :: ierr, i
  INTEGER :: Nprocs
  INTEGER :: my_rank

  CALL MPI_Init( ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, Nprocs, ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )

  IF( my_rank == 0 ) THEN 
    DO i = 1,4
      imsg(i) = i
    ENDDO 
  ELSE 
    DO i = 1,4
      imsg(i) = 0
    ENDDO 
  ENDIF 

  WRITE(*,*) 'Before: ', imsg(:)

  CALL MPI_Bcast( imsg, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

  WRITE(*,*) 'After: ', imsg(:)

  CALL MPI_Finalize( ierr )

END PROGRAM 

