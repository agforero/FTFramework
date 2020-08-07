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
subroutine vandermonde_coef_1d ( nd, xd, yd, cd )

!*****************************************************************************80
!
!! VANDERMONDE_COEF_1D computes coefficients of a 1D Vandermonde interpolant.
!
!  Discussion:
!
!    We assume the interpolant has the form
!
!      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
!
!    We have n data values (x(i),y(i)) which must be interpolated:
!
!      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
!
!    This can be cast as an NxN linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
!      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
!      [ ...................... ] [ ... ] = [ ... ]
!      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
!
!    and if the x values are distinct, the system is theoretically
!    invertible, so we can retrieve the coefficient vector c and
!    evaluate the interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), YD(ND), the data values.
!
!    Output, real ( kind = 8 ) CD(0:ND-1), the coefficients of the interpolating
!    polynomial.  CD(0) is the constant term, and CD(ND-1) multiplies X^(ND-1).
!
  implicit none

  integer ( kind = 4 ) nd

  real ( kind = 8 ) ad(nd,nd)
  real ( kind = 8 ) cd(0:nd-1)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)

  call vandermonde_matrix_1d ( nd, xd, ad )

  call qr_solve ( nd, nd, ad, yd, cd )

  return
end
subroutine vandermonde_matrix_1d ( nd, xd, ad )

!*****************************************************************************80
!
!! VANDERMONDE_MATRIX_1D computes a Vandermonde 1D interpolation matrix.
!
!  Discussion:
!
!    We assume the interpolant has the form
!
!      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
!
!    We have n data values (x(i),y(i)) which must be interpolated:
!
!      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
!
!    This can be cast as an NxN linear system for the polynomial
!    coefficients:
!
!      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
!      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
!      [ ...................... ] [ ... ] = [ ... ]
!      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
!
!    and if the x values are distinct, the matrix A is theoretically
!    invertible (though in fact, generally badly conditioned).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the data values.
!
!    Output, real ( kind = 8 ) AD(ND,ND), the Vandermonde matrix for X.
!
  implicit none

  integer ( kind = 4 ) nd

  real ( kind = 8 ) ad(nd,nd)
  integer ( kind = 4 ) j
  real ( kind = 8 ) xd(nd)

  ad(1:nd,1) = 1.0D+00
  do j = 2, nd
    ad(1:nd,j) = ad(1:nd,j-1) * xd(1:nd)
  end do

  return
end
subroutine vandermonde_value_1d ( nd, cd, ni, xi, yi )

!*****************************************************************************80
!
!! VANDERMONDE_VALUE_1D evaluates a Vandermonde 1D interpolant.
!
!  Discussion:
!
!    The polynomial 
!
!      p(x) = cd0 + cd1 * x + cd2 * x^2 + ... + cd(nd-1) * x^(nd-1)
!
!    is to be evaluated at the vector of NI values XI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data values.
!
!    Input, real ( kind = 8 ) CD(0:ND-1), the polynomial coefficients.  
!    CD(I) is the coefficient of X^I.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolation values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) cd(0:nd-1)
  integer ( kind = 4 ) i
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yi(ni)

  yi(1:ni) = cd(nd-1)
  do i = nd - 2, 0, -1
    yi(1:ni) = yi(1:ni) * xi(1:ni) + cd(i)
  end do

  return
end
