! ffr, June 2016
!
! calculate dot product between two vectors

PROGRAM dot
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: Ndata, Ndata_loc
  REAL(8), ALLOCATABLE :: vec1_loc(:), vec2_loc(:)
  INTEGER :: i, ierr, nprocs, my_rank, istart, iend
  REAL(8) :: res_loc, res

  Ndata = 10

  CALL MPI_Init( ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )

  Ndata_loc = Ndata/nprocs ! assume that nprocs divides Ndata evenly

  IF( my_rank == 0 ) THEN
    WRITE(*,*) 'Ndata_loc = ', Ndata_loc
  ENDIF

  ALLOCATE( vec1_loc(Ndata_loc) )
  ALLOCATE( vec2_loc(Ndata_loc) )
  
  istart = my_rank*Ndata_loc + 1
  iend   = (my_rank+1)*Ndata_loc

  WRITE(*,*) 'my_rank, istart, iend = ', my_rank, istart, iend

  DO i = 1, Ndata_loc
    vec1_loc(i) = 0.2d0*istart/Ndata
    vec2_loc(i) = -0.3d0*istart/Ndata
    WRITE(*,'(1x,A,2I4,2F18.10)') 'my_rank, i, vec1, vec2: ', my_rank, i, vec1_loc(i), vec2_loc(i)
    istart = istart + 1
  ENDDO

  res_loc = serial_dot( Ndata_loc, vec1_loc, vec2_loc )

  WRITE(*,*) 'my_rank, result = ', my_rank, res_loc

  CALL MPI_Reduce( res_loc, res, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

  IF( my_rank == 0 ) THEN
    WRITE(*,*) 'res = ', res
  ENDIF
  
  ! is this safe ?
  DEALLOCATE( vec1_loc )
  DEALLOCATE( vec2_loc )

  CALL MPI_Finalize( ierr )

CONTAINS

  REAL(8) FUNCTION serial_dot( N, x, y )
    IMPLICIT NONE
    INTEGER :: N
    REAL(8) :: x(N), y(N)
    REAL(8) :: ss
    ss = 0.d0
    DO i = 1, N
      ss = ss + x(i)*y(i)
    ENDDO
    serial_dot = ss
  END FUNCTION

END PROGRAM


