! ffr, June 2016

! use 3 calls to MPI_Bcast to distribute input data
! also use MPI_Reduce to compute final sum

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

  CALL get_data_bcast( a, b, N, my_rank )

  ! for debugging purpose
  !WRITE(*,'(1x,A,I4,2F18.10,I4)') 'my_rank, a, b, N', my_rank, a, b, N

  h = (b - a)/N

  local_N = N/nprocs
  local_a = a + my_rank * local_N * h
  local_b = local_a + local_N * h
  integral = Trapez( local_a, local_b, local_N, h )

  ! for debugging purpose
  !WRITE(*,'(1x,A,I4,F18.10)') 'my_rank, integral = ', my_rank, integral
 
  ! Add up integrals calculated by each process
  CALL MPI_Reduce( integral, total, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

  IF( my_rank == 0 ) THEN
    WRITE(*,*) 'With N = ', N, ' trapezoids, our estimate'
    WRITE(*,'(1x,A,F18.10,A,F18.10,A,F18.10)') 'of the integral from ', &
        a, ' to ', b, ' = ', total
  ENDIF

  CALL MPI_Finalize( ierr )

CONTAINS


SUBROUTINE get_data_bcast( a, b, n, my_rank )
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  REAL(8) :: a, b
  INTEGER :: n, my_rank
  !
  INTEGER :: iproc
  INTEGER :: dest, tag, ierr, source
  INTEGER :: stat(MPI_STATUS_SIZE)

  source = 0 ! master process
  IF( my_rank == 0 ) THEN
    WRITE(*,*)
    WRITE(*,*) 'Enter a, b, and N'
    READ(*,*) a, b, N
  ENDIF

  CALL MPI_Bcast( a, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_Bcast( b, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_Bcast( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

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

