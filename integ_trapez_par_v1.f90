! ffr, June 2016

PROGRAM integ_trapez
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  REAL(8) :: a, b
  INTEGER :: N, i
  REAL(8) :: integral, x, h
  INTEGER :: ierr, nprocs, dest, tag, stat(MPI_STATUS_SIZE), my_rank
  REAL(8) :: local_a, local_b, total
  INTEGER :: local_N
  INTEGER :: iproc

  a = 0.d0
  b = 1.d0
  N = 1024
  dest = 0
  tag  = 0


  CALL MPI_Init( ierr )
  CALL MPI_Comm_rank( MPI_COMM_WORLD, my_rank, ierr )
  CALL MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr )

  IF( my_rank == 0 ) THEN
    WRITE(*,*) 'A program to calculate definite integral from a to b'
    WRITE(*,*) 'using trapezoidal rule with N division'
  ENDIF

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

