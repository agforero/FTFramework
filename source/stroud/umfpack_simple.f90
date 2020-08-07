program main

!*****************************************************************************80
!
!  Purpose:
!
!    MAIN is the main program for UMFPACK_SIMPLE.
!
!  Discussion:
!
!    This program uses UMFPACK to solve the 5x5 linear system A*X=B:
!
!        2  3  0  0  0        1.0         8.0
!        3  0  4  0  6        2.0        45.0
!    A = 0 -1 -3  2  0    X = 3.0    B = -3.0
!        0  0  1  0  0        4.0         3.0
!        0  4  2  0  1        5.0        10.0
!
!    The matrix contains 12 nonzero values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Timothy Davis,
!    UMFPACK User Guide,
!    Version 5.6.2, 25 April 2013
!    http://suitesparse.com
!
  implicit none
 
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: nnz = 12

  integer ( kind = 4 ), dimension ( nnz ) :: ai = (/ &
    0, 1, 
    0, 2, 4, 
    1, 2, 3, 4, 
    2, 
    1, 4 /)
  integer ( kind = 4 ), dimension ( n + 1 ) :: ap = (/ &
    0, 2, 5, 9, 10, 12 /)
  real ( kind = 8 ), dimension ( nnz ) :: ax = (/ &
    2.0,  3.0, 
    3.0, -1.0, 4.0, 
    4.0, -3.0, 1.0, 2.0, 
    2.0, 
    6.0, 1.0 /)
  real ( kind = 8 ), dimension ( n ) :: b = (/ &
    8.0, 45.0, -3.0, 3.0, 19.0 /)
  real ( kind = 8 ) control(20)
  integer ( kind = 4 ) filenum
  integer ( kind = 4 ) i
  real ( kind = 8 ) info(90)
  integer ( kind = 8 ) numeric
  integer ( kind = 4 ) status
  integer ( kind = 8 ) symbolic
  integer ( kind = 4 ) sys
  real ( kind = 8 ) x(n)

  call timestamp ( );
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UMFPACK_SIMPLE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Use UMFPACK for the sparse linear system A*x=b.'
!
!  Set the default control parameters.
!
  call umf4def ( control )
!
!  From the matrix data, create the symbolic factorization information.
!
  call umf4sym ( n, n, ap, ai, ax, symbolic, control, info )

  if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4SYM returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  From the symbolic factorization information, carry out the numeric factorization.
!
  call umf4num ( ap, ai, ax, symbolic, numeric, control, info )

  if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4NUM returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  Free the memory associated with the symbolic factorization.
!
  call umf4fsym ( symbolic )
!
!  Solve the linear system.
!
  sys = 0
  call umf4sol ( sys, x, b, numeric, control, info )

  if ( info(1) < 0.0D+00 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'UMFPACK_SIMPLE - Fatal error!'
    write ( *, '(a,g14.6)' ) '  UMF4SOL returns INFO(1) = ', info(1)
    stop 1
  end if
!
!  Free the memory associated with the numeric factorization.
!
  call umf4fnum ( numeric )
!
!  Print the solution.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Computed solution:'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,g14.6)' ) x(i)
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'UMFPACK_SIMPLE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
