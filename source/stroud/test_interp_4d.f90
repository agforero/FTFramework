function p00_fun_4d ( prob, x )

!*****************************************************************************80
!
!! P00_FUN_4D evaluates the function for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the desired test problem.
!
!    Input, real ( kind = 8 ) X(4), the point at which the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) P00_FUN_4D, the value of the function at X.
!
  implicit none

  integer ( kind = 4 ) prob
  real ( kind = 8 ) p00_fun_4d
  real ( kind = 8 ) p01_fun_4d
  real ( kind = 8 ) p02_fun_4d
  real ( kind = 8 ) value
  real ( kind = 8 ) x(4)

  if ( prob == 1 ) then
    value = p01_fun_4d ( x )
  else if ( prob == 2 ) then
    value = p02_fun_4d ( x )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FUN_4D - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', prob
    value = 0.0D+00
    stop
  end if

  p00_fun_4d = value

  return
end
subroutine p00_lim_4d ( prob, a, b )

!*****************************************************************************80
!
!! P00_LIM_2D returns the interpolation interval for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the desired test problem.
!
!    Output, real ( kind = 8 ) A(4), B(4), the interpolation interval limits.
!
  implicit none

  real ( kind = 8 ) a(4)
  real ( kind = 8 ) b(4)
  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_lim_4d ( a, b )
  else if ( prob == 2 ) then
    call p02_lim_4d ( a, b )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_LIM_4D - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
subroutine p00_prob_num ( prob_num )

!*****************************************************************************80
!
!! P00_PROB_NUM returns the number of problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROB_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) prob_num

  prob_num = 2

  return
end
subroutine p00_story ( prob )

!*****************************************************************************80
!
!! P00_STORY prints the "story" for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
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

  integer ( kind = 4 ) prob

  if ( prob == 1 ) then
    call p01_story ( )
  else if ( prob == 2 ) then
    call p02_story ( )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_STORY - Fatal error!'
    write ( *, '(a)' ) '  Unexpected input value of PROB.'
    stop
  end if

  return
end
subroutine p00_title ( prob, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title of any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROB, the number of the desired test problem.
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  integer ( kind = 4 ) prob
  character ( len = * ) title

  if ( prob == 1 ) then
    call p01_title ( title )
  else if ( prob == 2 ) then
    call p02_title ( title )
  else
    title = ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal problem number = ', prob
    stop
  end if

  return
end
function p01_fun_4d ( x )

!*****************************************************************************80
!
!! P01_FUN_4D evaluates the function for problem 1.
!
!  Discussion:
!
!    This is a 4D version of the Runge example.
!
!  Interval:
!
!    -5 <= X(1:4) <= +5
!
!  Function:
!
!    1 / ( 1 + X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) P01_FUN_4D, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) p01_fun_4d
  real ( kind = 8 ) value
  real ( kind = 8 ) x(4)

  value = 1.0D+00 / ( x(1)**2 + x(2)**2 + x(3)**2 + x(4)**2 + 1.0D+00 )

  p01_fun_4d = value

  return
end
subroutine p01_lim_4d ( a, b )

!*****************************************************************************80
!
!! P01_LIM_4D returns the interpolation interval for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A(4), B(4), the interpolation interval limits.
!
  implicit none

  real ( kind = 8 ) a(4)
  real ( kind = 8 ) b(4)

  a(1:4) = -5.0D+00
  b(1:4) =  5.0D+00

  return
end
subroutine p01_story ( )

!*****************************************************************************80
!
!! P01_STORY prints the "story" for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This is a 4D version of Runge''s function.'
  write ( *, '(a)' ) '  In 1D, equally spaced interpolation nodes'
  write ( *, '(a)' ) '  result in a sequence of interpolants that'
  write ( *, '(a)' ) '  become highly oscillatory.'

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the title of problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = '4D Runge example, f(x,y) = 1 / ( w^2 + x^2 + y^2 + x^2 + 1 )'

  return
end
function p02_fun_4d ( x )

!*****************************************************************************80
!
!! P02_FUN_4D evaluates the function for problem 2.
!
!  Discussion:
!
!    This is an example from William Press.
!
!  Interval:
!
!    0 <= X(1:4) <= 1
!
!  Function:
!
!    F(X) = 514.1890 * exp ( - 2.0 * norm ( X - (0.3,0.3,0.3,0.3) ) )
!      * x(1) * ( 1 - x(1) ) * x(2) * ( 1 - x(2) ) 
!      * x(3) * ( 1 - x(3) ) * x(4) * ( 1 - x(4) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which the function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) P02_FUN_4D, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) p02_fun_4d
  real ( kind = 8 ) value
  real ( kind = 8 ) x(4)
  real ( kind = 8 ) xn

  xn = sqrt ( sum ( ( x(1:4) - 0.3D+00 )**2 ) )

  value = 514.1890D+00 * exp ( - 2.0D+00 * xn ) &
    * product ( x(1:4) ) * product ( 1.0D+00 - x(1:4) )

  p02_fun_4d = value

  return
end
subroutine p02_lim_4d ( a, b )

!*****************************************************************************80
!
!! P02_LIM_4D returns the interpolation interval for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A(4), B(4), the interpolation interval limits.
!
  implicit none

  real ( kind = 8 ) a(4)
  real ( kind = 8 ) b(4)

  a(1:4) = 0.0D+00
  b(1:4) = 1.0D+00

  return
end
subroutine p02_story ( )

!*****************************************************************************80
!
!! P02_STORY prints the "story" for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This is an example from William Press.'
  write ( *, '(a)' ) '  It is an offcenter Gaussian tapered in the'
  write ( *, '(a)' ) '  unit hypercube, zero at the edges.'

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns the title of problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'Press offcenter Gaussian'

  return
end
function r8_round ( x )

!*****************************************************************************80
!
!! R8_ROUND sets an R8 to the nearest integral value.
!
!  Example:
!
!        X        R8_ROUND
!
!      1.3         1.0
!      1.4         1.0
!      1.5         1.0 or 2.0
!      1.6         2.0
!      0.0         0.0
!     -0.7        -1.0
!     -1.1        -1.0
!     -1.6        -2.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value.
!
!    Output, real ( kind = 8 ) R8_ROUND, the rounded value.
!
  implicit none

  real ( kind = 8 ) r8_round
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    value = - real ( int ( - x + 0.5D+00 ), kind = 8 )
  else
    value =   real ( int ( + x + 0.5D+00 ), kind = 8 )
  end if

  r8_round = value

  return
end
subroutine r8vec_even_select ( n, xlo, xhi, ival, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
!
!    Unless N = 1, X(1) = XLO and X(N) = XHI.
!
!    If N = 1, then X(1) = 0.5*(XLO+XHI).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) XLO, XHI, the low and high values.
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any integer value.
!
!    Output, real ( kind = 8 ) XVAL, the IVAL-th of N evenly spaced values
!    between XLO and XHI.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ival
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xval

  if ( n == 1 ) then

    xval = 0.5D+00 * ( xlo + xhi )

  else

    xval = ( real ( n - ival,     kind = 8 ) * xlo   &
           + real (     ival - 1, kind = 8 ) * xhi ) &
           / real ( n        - 1, kind = 8 )

  end if

  return
end
subroutine r8vec_even2_select ( n, xlo, xhi, ival, xval )

!*****************************************************************************80
!
!! R8VEC_EVEN2_SELECT returns the I-th of N evenly spaced midpoint values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    This function returns the I-th of N evenly spaced midpoints of N
!    equal subintervals of [XLO,XHI].
!
!    XVAL = ( ( 2 * N - 2 * IVAL + 1 ) * XLO 
!           + (         2 * IVAL - 1 ) * XHI ) 
!           / ( 2 * N                )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) XLO, XHI, the low and high values.
!
!    Input, integer ( kind = 4 ) IVAL, the index of the desired point.
!    IVAL is normally between 1 and N, but may be any integer value.
!
!    Output, real ( kind = 8 ) XVAL, the IVAL-th of N evenly spaced midpoints
!    between XLO and XHI.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ival
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xval

  xval = ( real ( 2 * n - 2 * ival + 1, kind = 8 ) * xlo   &
         + real (         2 * ival - 1, kind = 8 ) * xhi ) &
         / real ( 2 * n, kind = 8 )

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
