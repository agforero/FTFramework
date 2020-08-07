function point_in_polygon ( n, x, y, x0, y0 )

!*****************************************************************************80
!
!! POINT_IN_POLYGON determines if a point is inside a polygon
!
!  Discussion:
!
!    If the points ( x(i), y(i) ) ( i = 1, 2, ..., n ) are,
!    in this cyclic order, the vertices of a simple closed polygon and
!    (x0,y0) is a point not on any side of the polygon, then the
!    procedure determines, by setting "point_in_polygon" to TRUE or FALSE,
!    whether (x0,y0) lies in the interior of the polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2016
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Moshe Shimrat,
!    ACM Algorithm 112,
!    Position of Point Relative to Polygon,
!    Communications of the ACM,
!    Volume 5, Number 8, page 434, August 1962.
!
!    Richard Hacker,
!    Certification of Algorithm 112,
!    Communications of the ACM,
!    Volume 5, Number 12, page  606, December 1962.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes or vertices in 
!    the polygon.  N must be at least 3.
!
!    Input, real ( kind = 8 ) V(2,N), the vertices of the polygon.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of the point to be tested.
!
!    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is
!    inside the polygon.
!
  implicit none

  integer ( kind = 4 ) n

  logical b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip1
  logical point_in_polygon
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
  
  b = .false.

  do i = 1, n
    ip1 = mod ( i, n ) + 1
!
!   if ( ( y(ip1) < y0 .and. y0 <= y(i)   ) .or. &
!        ( y(i)   < y0 .and. y0 <= y(ip1) ) ) then
!
    if ( y(ip1) < y0 .eqv. y0 <= y(i) ) then
      t = x0 - x(i) - ( y0 - y(i) ) * ( x(ip1) - x(i) ) / ( y(ip1) - y(i) )
      if ( t < 0.0D+00 ) then
        b = .not. b
      end if
    end if
  end do

  point_in_polygon = b

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
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
!    06 August 2005
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

