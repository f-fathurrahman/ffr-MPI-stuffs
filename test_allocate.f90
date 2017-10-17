!

PROGRAM test_allocate
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: ierr
  REAL(8), ALLOCATABLE :: vec1_loc(:), vec2_loc(:)
  INTEGER :: nprocs, my_rank
  INTEGER :: Ndata, Ndata_loc
  INTEGER :: i, istart, iend

  CALL MPI_Init( ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )

  Ndata = 10
  Ndata_loc = Ndata/nprocs

  istart = my_rank*Ndata_loc + 1
  iend   = (my_rank+1)*Ndata_loc
  WRITE(*,*) 'my_rank, istart, iend = ', my_rank, istart, iend

  ALLOCATE( vec1_loc(Ndata_loc) )
  ALLOCATE( vec2_loc(Ndata_loc) )
  DO i = 1, Ndata_loc
    vec1_loc(i) = 0.2d0*istart/Ndata
    vec2_loc(i) = -0.3d0*istart/Ndata
    WRITE(*,'(1x,A,2I4,2F18.10)') 'my_rrank, i, vec1, vec2: ', my_rank, i, vec1_loc(i), vec2_loc(i)
    istart = istart + 1
  ENDDO

  CALL MPI_Finalize( ierr )
END PROGRAM

