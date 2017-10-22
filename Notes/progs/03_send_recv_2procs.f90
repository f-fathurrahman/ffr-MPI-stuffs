PROGRAM test
  IMPLICIT NONE 
  INCLUDE 'mpif.h'
  INTEGER, PARAMETER :: Nbuf = 10
  REAL(8) :: buf(Nbuf)
  INTEGER :: ierr, dest, source, tag
  INTEGER :: Nprocs, rank
  INTEGER :: stat(MPI_STATUS_SIZE)
  INTEGER :: i

  CALL MPI_Init( ierr )

  CALL MPI_Comm_size( MPI_COMM_WORLD, Nprocs, ierr )
  IF( Nprocs /= 2 ) THEN 
    WRITE(*,*) 'ERROR: This program is only intended for two processors'
    WRITE(*,*) 'while this Nprocs = ', Nprocs
    WRITE(*,*) 'Please use Nprocs = 2'
    WRITE(*,*)
    CALL MPI_Finalize( ierr )
    STOP 
  ENDIF 

  CALL MPI_Comm_rank( MPI_COMM_WORLD, rank, ierr )

  ! rank 0 send data
  IF( rank == 0 ) THEN 
    ! buf is initialized here
    buf(:) = (/ 1.1d0, 2.1d0, 3.2d0, 4.1d0, 5.3d0, &
                6.0d0, 7.3d0, 8.8d0, 9.5d0, 10.3d0 /)
    dest = 1  ! to rank 1
    tag  = 0
    CALL MPI_Send( buf, Nbuf, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, ierr )
    !
    ! modify buf after sending data
    buf(1) = 1111d0
    buf(10) = 19191d0
    !
    WRITE(*,*)
    WRITE(*,*) 'Rank 0: After sending buf:'
    DO i = 1, Nbuf
      WRITE(*,'(1x,I4,I8,F18.10)') rank, i, buf(i)
    ENDDO 
  !
  ! rank 1 receive data
  !
  ELSE 
    source = 0  ! from rank 0
    tag    = 0
    !
    WRITE(*,*)
    WRITE(*,*) 'Rank 1: Before receiving buf:'
    DO i = 1, Nbuf
      WRITE(*,'(1x,I4,I8,F18.10)') rank, i, buf(i)
    ENDDO 
    !
    CALL MPI_Recv( buf, Nbuf, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, stat, ierr )
    !
    WRITE(*,*)
    WRITE(*,*) 'Rank 1: After receiving buf:'
    DO i = 1, Nbuf
      WRITE(*,'(1x,I4,I8,F18.10)') rank, i, buf(i)
    ENDDO 
  ENDIF 


  CALL MPI_Finalize( ierr )
  
  WRITE(*,*) 'Program ended normally'
END PROGRAM 
