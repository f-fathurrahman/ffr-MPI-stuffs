PROGRAM isend

  IMPLICIT NONE 
  include 'mpif.h'

  INTEGER :: istatus(MPI_STATUS_SIZE)
  INTEGER :: itag
  INTEGER :: isbuf, ireq, irbuf

  INTEGER :: ierr
  INTEGER :: Nprocs
  INTEGER :: my_rank

  CALL MPI_Init( ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, Nprocs, ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )

  itag = 1

  IF( my_rank == 0 ) THEN 
    isbuf = 9
    CALL MPI_Isend( isbuf, 1, MPI_INTEGER, 1, itag, MPI_COMM_WORLD, ireq, ierr )
    CALL MPI_Wait( ireq, istatus, ierr )
  ELSEIF( my_rank == 1 ) THEN 
    CALL MPI_Irecv( irbuf, 1, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, ireq, ierr )
    CALL MPI_Wait( ireq, istatus, ierr )
    WRITE(*,*) 'irbuf = ', irbuf
  ENDIF 

  CALL MPI_Finalize( ierr )

END PROGRAM 

