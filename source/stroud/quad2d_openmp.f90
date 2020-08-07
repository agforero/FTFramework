program main

!*****************************************************************************80
!
!! MAIN is the main program for QUAD2D_OPENMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2011
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  double precision a
  double precision b
  double precision error
  double precision exact
  external f
  double precision f
  integer i
  integer j
  integer n
  integer nx
  integer ny
  double precision pi
  double precision total
  double precision wtime
  double precision x
  double precision y

  a = 0.0
  b = 1.0
  nx = 32768
  ny = 32768
  n = nx * ny
  pi = 3.141592653589793D+00
  exact = pi * pi / 6.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD2D_OPENMP:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Estimate the integral of f(x,y) over [0,1]x[0,1].'
  write ( *, '(a)' ) '  f(x,y) = 1 / ( 1 - x * y ).'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  A        = ', a
  write ( *, * ) '  B        = ', b
  write ( *, * ) '  NX       = ', nx
  write ( *, * ) '  NY       = ', ny
  write ( *, * ) '  N        = ', n
  write ( *, * ) '  Exact    = ', exact

  wtime = omp_get_wtime ( )

  total = 0.0D+00
!$omp parallel shared ( a, b, nx, ny ) private ( i, j, x, y ) 

!$omp do reduction ( + : total )
  do i = 1, nx
    x = ( ( 2 * nx - 2 * i + 1 ) * a + ( 2 * i - 1 ) * b ) / ( 2 * nx )
    do j = 1, ny
      y = ( ( 2 * ny - 2 * j + 1 ) * a + ( 2 * j - 1 ) * b ) / ( 2 * ny )
      total = total + f ( x, y )
    end do
  end do
!$omp end do

!$omp end parallel

  wtime = omp_get_wtime ( ) - wtime

  total = ( b - a ) * ( b - a ) * total / dble ( nx ) / dble ( ny )
  error = abs ( total - exact )
 
  write ( *, '(a)' ) ' '
  write ( *, * ) '  Estimate = ', total
  write ( *, * ) '  Error    = ', error
  write ( *, * ) '  Time     = ', wtime
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD2D_OPENMP'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
function f ( x, y )

!*****************************************************************************80
!
!! F evaluates the function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  double precision f
  double precision x
  double precision y

  f = 1.0D+00 / ( 1.0D+00 - x * y )

  return
end
