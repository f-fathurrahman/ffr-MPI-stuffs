! ffr, June 2016

! use get_data subroutine to get input from user and send the input
! to other processes

PROGRAM integ_trapez
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  REAL(8) :: a, b
  INTEGER :: N
  REAL(8) :: integral, x, h
  INTEGER :: ierr, nprocs, dest, tag, stat, my_rank
  REAL(8) :: local_a, local_b, total
  INTEGER :: local_N
  INTEGER :: iproc

  dest = 0
  tag  = 0

  CALL MPI_Init( ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )

  IF( my_rank == 0 ) THEN
    WRITE(*,*) 'A program to calculate definite integral from a to b'
    WRITE(*,*) 'using trapezoidal rule with N division'
  ENDIF

  CALL get_data( a, b, N, my_rank, nprocs )

  ! for debugging purpose
  WRITE(*,'(1x,A,I4,2F18.10,I4)') 'my_rank, a, b, N', my_rank, a, b, N

  h = (b - a)/N

  local_N = N/nprocs
  local_a = a + my_rank * local_N * h
  local_b = local_a + local_N * h
  integral = Trapez( local_a, local_b, local_N, h )

  ! for debugging purpose
  WRITE(*,'(1x,A,I4,F18.10)') 'my_rank, integral = ', my_rank, integral
 
  IF( my_rank == 0 ) THEN
    total = integral
    DO iproc = 1, nprocs-1
      CALL MPI_Recv( integral, 1, MPI_REAL8, iproc, tag, MPI_COMM_WORLD, stat, ierr )
      total = total + integral
    ENDDO
  ELSE
    CALL MPI_Send( integral, 1, MPI_REAL8, dest, tag, MPI_COMM_WORLD, ierr )
  ENDIF

  IF( my_rank == 0 ) THEN
    WRITE(*,*) 'With N = ', N, ' trapezoids, our estimate'
    WRITE(*,'(1x,A,F18.10,A,F18.10,A,F18.10)') 'of the integral from ', &
        a, ' to ', b, ' = ', total
  ENDIF

  CALL MPI_Finalize( ierr )

CONTAINS


SUBROUTINE get_data( a, b, n, my_rank, nprocs )
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  REAL(8) :: a, b
  INTEGER :: n, my_rank, nprocs
  !
  INTEGER :: iproc
  INTEGER :: dest, tag, ierr, source
  INTEGER :: stat(MPI_STATUS_SIZE)

  source = 0 ! master process
  IF( my_rank == 0 ) THEN
    WRITE(*,*)
    WRITE(*,*) 'Enter a, b, and N'
    READ(*,*) a, b, N
    !
    DO iproc = 1, nprocs-1
      tag = 0
      CALL MPI_Send( a, 1, MPI_REAL8, iproc, tag, MPI_COMM_WORLD, ierr )
      tag = 1
      CALL MPI_Send( b, 1, MPI_REAL8, iproc, tag, MPI_COMM_WORLD, ierr )
      tag = 2
      CALL MPI_Send( N, 1, MPI_INTEGER, iproc, tag, MPI_COMM_WORLD, ierr )
    ENDDO
  ELSE
    tag = 0
    CALL MPI_Recv( a, 1, MPI_REAL8, source, tag, MPI_COMM_WORLD, stat, ierr )
    tag = 1
    CALL MPI_Recv( b, 1, MPI_REAL8, source, tag, MPI_COMM_WORLD, stat, ierr )
    tag = 2
    CALL MPI_Recv( N, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, stat, ierr )
  ENDIF
END SUBROUTINE


REAL(8) FUNCTION Trapez( local_a, local_b, local_N, h )
  IMPLICIT NONE
  !
  REAL(8) :: local_a, local_b, h
  INTEGER :: local_N
  !
  REAL(8) :: integral, x
  INTEGER :: i

  integral = 0.5d0 * ( func(local_a) + func(local_b) )
  x = local_a
  DO i = 1, local_N-1
    x = x + h
    integral = integral + func(x)
  ENDDO
  Trapez = integral*h
END FUNCTION

REAL(8) FUNCTION func(x)
  IMPLICIT NONE
  REAL(8) :: x
  func = x*x
END FUNCTION

END PROGRAM

