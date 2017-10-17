! ffr, June 2016

! use hand-coded tree-structured broadcast
!
! TODO: not working yet

PROGRAM integ_trapez
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  REAL(8) :: a, b
  INTEGER :: Ne
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

  CALL get_data_tree( a, b, N, my_rank, nprocs )

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


! use unsigned so that right shift will fill the leftmost bit with 0
INTEGER FUNCTION ceiling_log2(x)
  IMPLICIT NONE
  INTEGER :: x
  INTEGER :: res, temp

  temp = x - 1
  res = 0
  DO WHILE( temp /= 0 )
    temp = temp / 2
    res = res + 1
  ENDDO
  ceiling_log2 = res
END FUNCTION


INTEGER FUNCTION I_receive( stage, my_rank, source )
  IMPLICIT NONE
  INTEGER :: stage, my_rank, source
  !
  INTEGER :: power_2_stage
  INTEGER :: i

  ! 2^stage = 1 << stage
  power_2_stage = 1
  DO i = 1, stage
    power_2_stage = power_2_stage * 2
  ENDDO

  IF( power_2_stage <= my_rank .AND. &
      my_rank < 2*power_2_stage ) THEN
    source = my_rank - power_2_stage
    I_receive = 1
  ELSE
    I_receive = 0
  ENDIF

END FUNCTION


INTEGER FUNCTION I_send( stage, my_rank, nprocs, dest )
  IMPLICIT NONE
  INTEGER :: stage, my_rank, nprocs, dest
  !
  INTEGER :: power_2_stage
  INTEGER :: i

  ! 2^stage = 1 << stage
  power_2_stage = 1
  DO i = 1, stage
    power_2_stage = power_2_stage * 2
  ENDDO
  !
  IF( my_rank < power_2_stage ) THEN
    dest = my_rank + power_2_stage
    IF( dest >= 0 ) THEN
      I_send = 0
    ELSE
      I_send = 1
    ENDIF
  ELSE
    I_send = 0
  ENDIF
END FUNCTION


SUBROUTINE send( a, b, n, dest )
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  !
  REAL(8) :: a, b
  INTEGER :: n
  INTEGER :: dest
  !
  INTEGER :: ierr

  CALL MPI_Send( a, 1, MPI_REAL8, dest, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_Send( b, 1, MPI_REAL8, dest, 1, MPI_COMM_WORLD, ierr )
  CALL MPI_Send( n, 1, MPI_INTEGER, dest, 2, MPI_COMM_WORLD, ierr )
END SUBROUTINE


SUBROUTINE receive( a, b, n, source )
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  !
  REAL(8) :: a, b
  INTEGER :: n
  INTEGER :: source
  !
  INTEGER :: stat(MPI_STATUS_SIZE)
  INTEGER :: ierr

  CALL MPI_Recv( a, 1, MPI_REAL8, source, 0, MPI_COMM_WORLD, stat, ierr )
  CALL MPI_Recv( b, 1, MPI_REAL8, source, 1, MPI_COMM_WORLD, stat, ierr )
  CALL MPI_Recv( n, 1, MPI_INTEGER, source, 2, MPI_COMM_WORLD, stat, ierr )
END SUBROUTINE


! Process 0 prompts user for input and reads in the values
! Process 0 sends input values to other process using hand-coded
! tree-structured broadcast
SUBROUTINE get_data_tree( a, b, n, my_rank, nprocs )
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  REAL(8) :: a, b
  INTEGER :: n, my_rank, nprocs
  !
  INTEGER :: iproc
  INTEGER :: dest, tag, ierr, source, stage
  INTEGER :: stat(MPI_STATUS_SIZE)

  source = 0 ! master process
  IF( my_rank == 0 ) THEN
    WRITE(*,*)
    WRITE(*,*) 'Enter a, b, and N'
    READ(*,*) a, b, N
  ENDIF
  !
  DO stage = 1, ceiling_log2(nprocs)-1
    IF( I_receive( stage, my_rank, source ) == 1 ) THEN
      CALL receive( a, b, n, source )
    ELSE
      CALL send( a, b, n, dest )
    ENDIF
  ENDDO
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

