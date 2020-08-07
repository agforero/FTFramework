subroutine cheby_van1 ( m, a, b, n, x, v )

!*****************************************************************************80
!
!! CHEBY_VAN1 returns the Chebyshev Vandermonde-like matrix for [A,B].
!
!  Discussion:
!
!    Normally, the Chebyshev polynomials are defined on -1 <= XI <= +1.
!    Here, we assume the Chebyshev polynomials have been defined on the
!    interval A <= X <= B, using the mapping
!      XI = ( - ( B - X ) + ( X - A ) ) / ( B - A )
!    so that
!      ChebyAB(A,B;X) = Cheby(XI).
!
!    if ( I == 1 ) then
!      V(1,1:N) = 1;
!    elseif ( I == 2 ) then
!      V(2,1:N) = XI(1:N);
!    else
!      V(I,1:N) = 2.0 * XI(1:N) * V(I-1,1:N) - V(I-2,1:N);
!
!  Example:
!
!    M = 5, N = 5, X = ( 1, 2, 3, 4, 5 )
!
!    1  1   1    1    1
!    1  2   3    4    5
!    1  7  17   31   49
!    1 26  99  244  485
!    1 97 577 1921 4801
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Nicholas Higham,
!    Stability analysis of algorithms for solving confluent
!    Vandermonde-like systems,
!    SIAM Journal on Matrix Analysis and Applications,
!    Volume 11, 1990, pages 23-41.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, real ( kind = 8 ) A, B, the interval.
!
!    Input, integer ( kind = 4 ) N, the number of values in X, and the number
!    of columns in the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector that defines A.
!
!    Output, real ( kind = 8 ) V(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xi(n)
  real ( kind = 8 ) v(m,n)
!
!  Compute the normalized abscissas in [-1,+1].
!
  xi(1:n) = ( - 1.0D+00 * ( b - x(1:n)     )   &
              + 1.0D+00 * (     x(1:n) - a ) ) &
            /             ( b          - a )

  if ( 1 <= m ) then
    v(1,1:n) = 1.0D+00
  end if

  if ( 2 <= m ) then
    v(2,1:n) = xi(1:n)
  end if

  do i = 3, m
    v(i,1:n) = 2.0D+00 * xi(1:n) * v(i-1,1:n) - v(i-2,1:n)
  end do

  return
end
subroutine legendre_van ( m, a, b, n, x, v )

!*****************************************************************************80
!
!! LEGENDRE_VAN returns the LEGENDRE_VAN matrix.
!
!  Discussion:
!
!    The LEGENDRE_VAN matrix is the Legendre Vandermonde-like matrix.
!
!    Normally, the Legendre polynomials are defined on -1 <= XI <= +1.
!    Here, we assume the Legendre polynomials have been defined on the
!    interval A <= X <= B, using the mapping
!      XI = ( - ( B - X ) + ( X - A ) ) / ( B - A )
!    so that
!      Lab(A,B;X) = L(XI).
!
!    if ( I = 1 ) then
!      V(1,1:N) = 1
!    else if ( I = 2 ) then
!      V(2,1:N) = XI(1:N)
!    else
!      V(I,1:N) = ( (2*I-1) * XI(1:N) * V(I-1,1:N) - (I-1)*V(I-2,1:N) ) / I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, real ( kind = 8 ) A, B, the interval.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector that defines the matrix.
!
!    Output, real ( kind = 8 ) V(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) v(m,n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xi

  do j = 1, n

    xi = ( - ( b - x(j) ) + ( x(j) - a ) ) / ( b - a )

    do i = 1, m

      if ( i == 1 ) then
        v(i,j) = 1.0D+00
      else if ( i == 2 ) then
        v(i,j) = xi
      else
        v(i,j) = ( real ( 2 * i - 1, kind = 8 ) * xi * v(i-1,j) + &
                   real (   - i + 1, kind = 8 ) *      v(i-2,j) ) &
                 / real (     i, kind = 8 )
      end if

    end do
  end do

  return
end
subroutine line_fekete_chebyshev ( m, a, b, n, x, nf, xf, wf )

!*****************************************************************************80
!
!! LINE_FEKETE_CHEBYSHEV: approximate Fekete points in an interval [A,B].
!
!  Discussion:
!
!    We use the Chebyshev basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Len Bos, Norm Levenberg,
!    On the calculation of approximate Fekete points: the univariate case,
!    Electronic Transactions on Numerical Analysis, 
!    Volume 30, pages 377-397, 2008.
!    
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of basis polynomials.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!    M <= N.
!
!    Input, real ( kind = 8 ) X(N), the coordinates of the sample points.
!
!    Output, integer ( kind = 4 ) NF, the number of Fekete points.
!    If the computation is successful, NF = M.
!
!    Output, real ( kind = 8 ) XF(NF), the coordinates of the Fekete points.
!
!    Output, real ( kind = 8 ) WF(NF), the weights of the Fekete points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) j
  real ( kind = 8 ) mom(m)
  integer ( kind = 4 ) nf
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) v(m,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) wf(m)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xf(m)

  if ( n < m ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LINE_FEKETE_CHEBYSHEV - Fatal error!'
    write ( *, '(a)' ) '  N < M.'
    stop 1
  end if
!
!  Compute the Chebyshev-Vandermonde matrix.
!
  call cheby_van1 ( m, a, b, n, x, v )
!
!  MOM(I) = Integral ( A <= x <= B ) Tab(A,B,I;x) dx
!
  mom(1) = r8_pi * ( b - a ) / 2.0D+00
  mom(2:m) = 0.0D+00
!
!  Solve the system for the weights W.
!
  call qr_solve ( m, n, v, mom, w )
!
!  Extract the data associated with the nonzero weights.
!
  nf = 0
  do j = 1, n
    if ( w(j) /= 0.0D+00 ) then
      if ( nf < m ) then
        nf = nf + 1
        xf(nf) = x(j)
        wf(nf) = w(j)
      end if
    end if
  end do

  return
end
subroutine line_fekete_legendre ( m, a, b, n, x, nf, xf, wf )

!*****************************************************************************80
!
!! LINE_FEKETE_LEGENDRE computes approximate Fekete points in an interval [A,B].
!
!  Discussion:
!
!    We use the uniform weight and the Legendre basis:
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Len Bos, Norm Levenberg,
!    On the calculation of approximate Fekete points: the univariate case,
!    Electronic Transactions on Numerical Analysis, 
!    Volume 30, pages 377-397, 2008.
!    
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of basis polynomials.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!    M <= N.
!
!    Input, real ( kind = 8 ) X(N), the coordinates of the sample points.
!
!    Output, integer ( kind = 4 ) NF, the number of Fekete points.
!    If the computation is successful, NF = M.
!
!    Output, real ( kind = 8 ) XF(NF), the coordinates of the Fekete points.
!
!    Output, real ( kind = 8 ) WF(NF), the weights of the Fekete points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) j
  real ( kind = 8 ) mom(m)
  integer ( kind = 4 ) nf
  real ( kind = 8 ) v(m,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) wf(m)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xf(m)

  if ( n < m ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LINE_FEKETE_LEGENDRE - Fatal error!'
    write ( *, '(a)' ) '  N < M.'
    stop 1
  end if
!
!  Compute the Legendre-Vandermonde matrix.
!
  call legendre_van ( m, a, b, n, x, v )
!
!  MOM(i) = integral ( A <= X <= B ) Lab(A,B,I;X) dx
!
  mom(1) = b - a
  mom(2:m) = 0.0D+00
!
!  Solve the system for the weights W.
!
  call qr_solve ( m, n, v, mom, w )
!
!  Extract the data associated with the nonzero weights.
!
  nf = 0
  do j = 1, n
    if ( w(j) /= 0.0D+00 ) then
      if ( nf < m ) then
        nf = nf + 1
        xf(nf) = x(j)
        wf(nf) = w(j)
      end if
    end if
  end do

  return
end
subroutine line_fekete_monomial ( m, a, b, n, x, nf, xf, wf )

!*****************************************************************************80
!
!! LINE_FEKETE_MONOMIAL computes approximate Fekete points in an interval [A,B].
!
!  Discussion:
!
!    We use the uniform weight and the monomial basis:
!
!      P(j) = x^(j-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alvise Sommariva, Marco Vianello,
!    Computing approximate Fekete points by QR factorizations of Vandermonde 
!    matrices,
!    Computers and Mathematics with Applications,
!    Volume 57, 2009, pages 1324-1336.
!    
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of basis polynomials.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!    M <= N.
!
!    Input, real ( kind = 8 ) X(N), the coordinates of the sample points.
!
!    Output, integer ( kind = 4 ) NF, the number of Fekete points.
!    If the computation is successful, NF = M.
!
!    Output, real ( kind = 8 ) XF(NF), the coordinates of the Fekete points.
!
!    Output, real ( kind = 8 ) WF(NF), the weights of the Fekete points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) ai
  real ( kind = 8 ) b
  real ( kind = 8 ) bi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) mom(m)
  integer ( kind = 4 ) nf
  real ( kind = 8 ) v(m,n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) wf(m)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xf(m)

  if ( n < m ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'LINE_FEKETE_MONOMIAL - Fatal error!'
    write ( *, '(a)' ) '  N < M.'
    stop 1
  end if
!
!  Form the moments.
!
  call line_monomial_moments ( a, b, m, mom )
!
!  Form the rectangular Vandermonde matrix V for the polynomial basis.
!
  v(1,1:n) = 1.0D+00
  do i = 2, m
    v(i,1:n) = v(i-1,1:n) * x(1:n)
  end do
!
!  Solve the system for the weights W.
!
  call qr_solve ( m, n, v, mom, w )
!
!  Extract the data associated with the nonzero weights.
!
  nf = 0
  do j = 1, n
    if ( w(j) /= 0.0D+00 ) then
      if ( nf < m ) then
        nf = nf + 1
        xf(nf) = x(j)
        wf(nf) = w(j)
      end if
    end if
  end do

  return
end
subroutine line_monomial_moments ( a, b, m, mom )

!*****************************************************************************80
!
!! LINE_MONOMIAL_MOMENTS computes monomial moments in [A,B].
!
!  Discussion:
!
!    We use the uniform weight and the shifted and scaled monomial basis:
!
!      P(a,b,i;x) = xi(a,b;x)^(i-1)
!       xi(a,b;x) = ( - ( b - x ) + ( x - a ) ) / ( b - a )
!
!    The i-th moment is
!
!      mom(i) = integral ( a <= x <= b ) P(a,b,i;x) dx
!             = integral ( a <= x <= b ) xi(a,b;x)^(i-1) dx
!             = 0.5 * ( b - a ) * integral ( -1 <= xi <= +1 ) xi^(i-1) dxi
!             = 0.5 * ( b - a ) * xi^i / i | ( -1 <= xi <= +1 )
!             = 0.5 * ( b - a ) * ( 1 - (-1)^i ) / i
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) M, the number of basis polynomials.
!
!    Output, real ( kind = 8 ) MOM(M), the moments.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) mom(m)

  do i = 1, m
    mom(i) = ( b - a ) * real ( mod ( i, 2 ), kind = 8 ) / real ( i, kind = 8 )
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = 8 ) * a   &
             + real (     i - 1, kind = 8 ) * b ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
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
