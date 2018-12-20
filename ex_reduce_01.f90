PROGRAM ex_reduce
  IMPLICIT NONE 
  INCLUDE 'mpif.h'

  INTEGER :: istart, iend
  REAL(8) :: a(9), ss, tmp

  INTEGER :: ierr, i
  INTEGER :: Nprocs
  INTEGER :: my_rank

  CALL MPI_Init( ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, Nprocs, ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )

  istart = my_rank*3 + 1
  iend   = istart + 2

  DO i = istart,iend
    a(i) = i
  ENDDO 

  IF( my_rank == 0 ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'my_rank = ', my_rank
    DO i = 1,9
      WRITE(*,'(1x,I4,A,F18.10)') i, ' ', a(i)
    ENDDO 
  ELSEIF( my_rank == 1 ) THEN 
    WRITE(*,*)
    WRITE(*,*) 'my_rank = ', my_rank
    DO i = 1,9
      WRITE(*,'(1x,I4,A,F18.10)') i, ' ', a(i)
    ENDDO 
  ELSE 
    WRITE(*,*)
    WRITE(*,*) 'my_rank = ', my_rank
    DO i = 1,9
      WRITE(*,'(1x,I4,A,F18.10)') i, ' ', a(i)
    ENDDO 
  ENDIF 

  ss = 0.d0
  DO i = istart,iend
    ss = ss + a(i)
  ENDDO 

  CALL MPI_Reduce( ss, tmp, 1, MPI_REAL8, MPI_SUM, 0, &
                   MPI_COMM_WORLD, ierr )

  ss = tmp
  WRITE(*,'(1x,A,I4,A,F18.10)') 'my_rank = ', my_rank, ' ss = ', ss

  IF( my_rank == 0 ) THEN 
    WRITE(*,*) 'sum = ', ss
  ENDIF 

  CALL MPI_Finalize( ierr )

END PROGRAM 
