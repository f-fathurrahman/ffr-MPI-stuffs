! ffr, June 2016

PROGRAM hello_v3
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER :: ierr
  INTEGER :: nprocs
  INTEGER :: my_rank
  INTEGER :: tag, dest, siz, source
  INTEGER :: iproc
  INTEGER :: stat(MPI_STATUS_SIZE)
  CHARACTER(100) :: message
  CHARACTER :: my_rank_str ! works only for nprocs < 10

  CALL MPI_Init( ierr )

  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )

  IF( my_rank /= 0 ) THEN
    ! This only works for nproc < 10
    WRITE(my_rank_str,'(I1)') my_rank
    message = 'This is hello message from process ' // my_rank_str
    !WRITE(*,*) len_trim(message)
    !WRITE(*,*) message
    dest = 0
    tag  = 0
    CALL MPI_Send( message, len_trim(message), MPI_CHARACTER, dest, tag, &
            MPI_COMM_WORLD, ierr )
  ELSE
    !WRITE(*,*) 'This is proc = ', my_rank
    DO iproc = 1, nprocs-1
      tag = 0
      source = iproc
      CALL MPI_Recv( message, 100, MPI_CHARACTER, source, tag, &
            MPI_COMM_WORLD, stat, ierr )
      CALL MPI_Get_count( stat, MPI_CHARACTER, siz, ierr )
      !WRITE(*,*) 'siz = ', siz
      WRITE(*,'(A)') message(1:siz)
    ENDDO
  ENDIF

  CALL MPI_Finalize( ierr )

END PROGRAM

