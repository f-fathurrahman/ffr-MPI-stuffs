! ffr, June 2016

PROGRAM integ_trapez
  IMPLICIT NONE
  REAL(8) :: a, b
  INTEGER :: N, i
  REAL(8) :: integral, x, h

  WRITE(*,*) 'A program to calculate definite integral from a to b'
  WRITE(*,*) 'using trapezoidal rule with N division'
  WRITE(*,*)
  WRITE(*,*) 'Please enter a, b, N'
  READ(*,*) a, b, N

  h = (b - a)/N
  integral = ( func(a) + func(b) ) /2.d0
  x = a
  DO i = 1, N-1
    x = x + h
    integral = integral + func(x)
  ENDDO
  integral = integral*h

  WRITE(*,*) 'With N = ', N, ' trapezoids, our estimate'
  WRITE(*,'(1x,A,F18.10,A,F18.10,A,F18.10)') 'of the integral from ', a, ' to ', b, ' = ', integral


CONTAINS

REAL(8) FUNCTION func(x)
  IMPLICIT NONE
  REAL(8) :: x
  func = x*x
END FUNCTION

END PROGRAM

