! ffr, June 2016
!
! calculate dot product between two vectors

PROGRAM dot
  IMPLICIT NONE
  INTEGER :: Ndata
  REAL(8), ALLOCATABLE :: vec1(:), vec2(:)
  INTEGER :: i
  REAL(8) :: res

  Ndata = 1000

  ALLOCATE( vec1(Ndata) )
  ALLOCATE( vec2(Ndata) )
  
  DO i = 1, Ndata
    vec1(i) = 0.2d0*i/Ndata
    vec2(i) = -0.3d0*i/Ndata
  ENDDO

  res = serial_dot( Ndata, vec1, vec2 )

  WRITE(*,*) 'result = ', res
  
  ! is this safe ?
  DEALLOCATE( vec1 )
  DEALLOCATE( vec2 )

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


