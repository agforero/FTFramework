function c8_le_l2 ( x, y )

!*****************************************************************************80
!
!! C8_LE_L2 := X <= Y for the L2 norm on C8 values.
!
!  Discussion:
!
!    The L2 norm can be defined here as:
!
!      C8_NORM2(X) = sqrt ( ( real (X) )^2 + ( imag (X) )^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
! 
!  Parameters:
!
!    Input, complex ( kind = 8 ) X, Y, the values to be compared.
!
!    Output, logical C8_LE_L2, is TRUE if X <= Y.
!
  implicit none

  complex ( kind = 8 ) x
  complex ( kind = 8 ) y
  logical c8_le_l2

  if ( ( real ( x, kind = 8 ) )**2 + ( imag ( x ) )**2 <= &
       ( real ( y, kind = 8 ) )**2 + ( imag ( y ) )**2 ) then
    c8_le_l2 = .true.
  else
    c8_le_l2 = .false.
  end if

  return
end
subroutine c8_swap ( x, y )

!*****************************************************************************80
!
!! C8_SWAP swaps two C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, complex ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  complex ( kind = 8 ) x
  complex ( kind = 8 ) y
  complex ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine c83_cr_fa ( n, a, a_cr )

!*****************************************************************************80
!
!! C83_CR_FA decomposes a C83 matrix using cyclic reduction.
!
!  Discussion:
!
!    The C83 storage format is used for a real tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Once C83_CR_FA has decomposed a matrix A, then C83_CR_SL may be used 
!    to solve linear systems A * x = b.
!
!    C83_CR_FA does not employ pivoting.  Hence, the results can be more
!    sensitive to ill-conditioning than standard Gauss elimination.  In
!    particular, C83_CR_FA will fail if any diagonal element of the matrix
!    is zero.  Other matrices may also cause C83_CR_FA to fail.
!
!    C83_CR_FA can be guaranteed to work properly if the matrix is strictly
!    diagonally dominant, that is, if the absolute value of the diagonal
!    element is strictly greater than the sum of the absolute values of
!    the offdiagonal elements, for each equation.
!
!    The algorithm may be illustrated by the following figures:
!
!    The initial matrix is given by:
!
!          D1 U1
!          L1 D2 U2
!             L2 D3 U3
!                L3 D4 U4
!                   L4 D5 U5
!                      L5 D6
!
!    Rows and columns are permuted in an odd/even way to yield:
!
!          D1       U1
!             D3    L2 U3
!                D5    L4 U5
!          L1 U2    D2
!             L3 U4    D4
!                L5       D6
!
!    A block LU decomposition is performed to yield:
!
!          D1      |U1
!             D3   |L2 U3
!                D5|   L4 U5
!          --------+--------
!                  |D2'F3
!                  |F1 D4'F4
!                  |   F2 D6'
!
!    For large systems, this reduction is repeated on the lower right hand
!    tridiagonal subsystem until a completely upper triangular system
!    is obtained.  The system has now been factored into the product of a
!    lower triangular system and an upper triangular one, and the information
!    defining this factorization may be used by C83_CR_SL to solve linear
!    systems.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(3,N), the matrix.
!
!    Output, complex ( kind = 8 ) A_CR(3,0:2*N), factorization information 
!    needed by C83_CR_SL.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) a_cr(3,0:2*n)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulp
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ilp
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C83_CR_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop 1
  end if

  if ( n == 1 ) then
    a_cr(1,0:2) = 0.0D+00
    a_cr(2,0) = 0.0D+00
    a_cr(2,1) = 1.0D+00 / a(2,1)
    a_cr(2,2) = 0.0D+00
    a_cr(3,0:2) = 0.0D+00
    return
  end if
!
!  Zero out the workspace entries.
!
  a_cr(1,0) = 0.0D+00
  a_cr(1,1:n-1) = a(1,2:n)
  a_cr(1,n:2*n) = 0.0D+00

  a_cr(2,0) = 0.0D+00
  a_cr(2,1:n) = a(2,1:n)
  a_cr(2,n+1:2*n) = 0.0D+00

  a_cr(3,0) = 0.0D+00
  a_cr(3,1:n-1) = a(3,1:n-1)
  a_cr(3,n:2*n) = 0.0D+00

  il = n
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    if ( mod ( il, 2 ) == 1 ) then
      inc = il + 1
    else
      inc = il
    end if

    incr = inc / 2
    il = il / 2
    ihaf = ipntp + incr + 1
    ifulp = ipnt + inc + 2

!dir$ ivdep
    do ilp = incr, 1, -1
      ifulp = ifulp - 2
      iful = ifulp - 1
      ihaf = ihaf - 1
      a_cr(2,iful) = 1.0D+00 / a_cr(2,iful)
      a_cr(3,iful)  = a_cr(3,iful)  * a_cr(2,iful)
      a_cr(1,ifulp) = a_cr(1,ifulp) * a_cr(2,ifulp+1)
      a_cr(2,ihaf)  = a_cr(2,ifulp) - a_cr(1,iful)  * a_cr(3,iful) &
                                  - a_cr(1,ifulp) * a_cr(3,ifulp)
      a_cr(3,ihaf) = -a_cr(3,ifulp) * a_cr(3,ifulp+1)
      a_cr(1,ihaf) = -a_cr(1,ifulp) * a_cr(1,ifulp+1)
    end do

  end do

  a_cr(2,ipntp+1) = 1.0D+00 / a_cr(2,ipntp+1)

  return
end
subroutine c83_cr_sls ( n, a_cr, nb, b, x )

!*****************************************************************************80
!
!! C83_CR_SLS solves several linear systems factored by C83_CR_FA.
!
!  Discussion:
!
!    The matrix A must be tridiagonal.  C83_CR_FA is called to compute the
!    LU factors of A.  It does so using a form of cyclic reduction.  If
!    the factors computed by C83_CR_FA are passed to C83_CR_SL, then one or many
!    linear systems involving the matrix A may be solved.
!
!    Note that C83_CR_FA does not perform pivoting, and so the solution 
!    produced by C83_CR_SL may be less accurate than a solution produced 
!    by a standard Gauss algorithm.  However, such problems can be 
!    guaranteed not to occur if the matrix A is strictly diagonally 
!    dominant, that is, if the absolute value of the diagonal coefficient 
!    is greater than the sum of the absolute values of the two off diagonal 
!    coefficients, for each row of the matrix.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Hockney,
!    A fast direct solution of Poisson's equation using Fourier Analysis,
!    Journal of the ACM,
!    Volume 12, Number 1, pages 95-113, January 1965.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A_CR(3,0:2*N), factorization information 
!    computed by C83_CR_FA.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input, real ( kind = 8 ) B(N,NB), the right hand sides.
!
!    Output, real ( kind = 8 ) X(N,NB), the solutions of the linear systems.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  complex ( kind = 8 ) a_cr(3,0:2*n)
  complex ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) iful
  integer ( kind = 4 ) ifulm
  integer ( kind = 4 ) ihaf
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ipnt
  integer ( kind = 4 ) ipntp
  integer ( kind = 4 ) ndiv
  complex ( kind = 8 ) rhs(0:2*n,nb)
  complex ( kind = 8 ) x(n,nb)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C83_CR_SLS - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonpositive N = ', n
    stop 1
  end if

  if ( n == 1 ) then
    x(1,1:nb) = a_cr(2,1) * b(1,1:nb)
    return
  end if
!
!  Set up RHS.
!
  rhs(0,1:nb) = 0.0D+00
  rhs(1:n,1:nb) = b(1:n,1:nb)
  rhs(n+1:2*n,1:nb) = 0.0D+00

  il = n
  ndiv = 1
  ipntp = 0

  do while ( 1 < il )

    ipnt = ipntp
    ipntp = ipntp + il
    il = il / 2
    ndiv = ndiv * 2
    ihaf = ipntp

!dir$ ivdep
    do iful = ipnt + 2, ipntp, 2
      ihaf = ihaf + 1
      rhs(ihaf,1:nb) = rhs(iful,1:nb) &
        - a_cr(3,iful-1) * rhs(iful-1,1:nb) &
        - a_cr(1,iful)   * rhs(iful+1,1:nb)
    end do

  end do

  rhs(ihaf,1:nb) = rhs(ihaf,1:nb) * a_cr(2,ihaf)
  ipnt = ipntp

  do while ( 0 < ipnt )

    ipntp = ipnt
    ndiv = ndiv / 2
    il = n / ndiv
    ipnt = ipnt - il
    ihaf = ipntp

!dir$ ivdep
    do ifulm = ipnt + 1, ipntp, 2
      iful = ifulm + 1
      ihaf = ihaf + 1
      rhs(iful,1:nb) = rhs(ihaf,1:nb)
      rhs(ifulm,1:nb) = a_cr(2,ifulm) &
        * (                     rhs(ifulm,1:nb) &
            - a_cr(3,ifulm-1) * rhs(ifulm-1,1:nb) &
            - a_cr(1,ifulm)   * rhs(iful,1:nb) )
    end do

  end do

  x(1:n,1:nb) = rhs(1:n,1:nb)

  return
end
subroutine c83_indicator ( n, a )

!*****************************************************************************80
!
!! C83_INDICATOR sets up a C83 indicator matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Output, complex ( kind = 8 ) A(3,N), the indicator matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  a(1,1) = 0.0D+00
  do j = 2, n
    i = j - 1
    a(1,j) = cmplx ( i, j, kind = 8 )
  end do

  do j = 1, n
    i = j
    a(2,j) = cmplx ( i, j, kind = 8 )
  end do

  do j = 1, n - 1
    i = j + 1
    a(3,j) = cmplx ( i, j, kind = 8 )
  end do
  a(3,n) = 0.0D+00

  return
end
subroutine c83_jac_sl ( n, a, b, x, tol, it_max, job, it, diff )

!*****************************************************************************80
!
!! C83_JAC_SL solves a C83 system using Jacobi iteration.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A(3,N), the matrix.
!
!    Input, complex ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Input/output, complex ( kind = 8 ) X(N), an approximate solution 
!    to the system.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.  If the maximum change in
!    the solution is less than TOL, the iteration is terminated early.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer ( kind = 4 ) IT, the number of iterations taken.
!
!    Output, real ( kind = 8 ) DIFF, the maximum change in the solution
!    on the last iteration.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) it_num
  integer ( kind = 4 ) job
  real ( kind = 8 ) tol
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) x_new(n)
  real ( kind = 8 ) x_norm
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C83_JAC_SL - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
      stop 1
    end if
  end do

  if ( job == 0 ) then

    do it_num = 1, it_max

      it = it_num

      x_new(1) =   b(1)                   - a(3,1) * x(2)
      do i = 2, n - 1
        x_new(i) = b(i) - a(1,i) * x(i-1) - a(3,i) * x(i+1)
      end do
      x_new(n) =   b(n) - a(1,n) * x(n-1)
!
!  Divide by diagonal terms.
!
      x_new(1:n) = x_new(1:n) / a(2,1:n)
!
!  Measure norms of solution and change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update.
!
      x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  else

    do it_num = 1, it_max

      it = it_num

      x_new(1) =   b(1)                     - a(1,2) * x(2)
      do i = 2, n - 1
        x_new(i) = b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1)
      end do
      x_new(n) =   b(n) - a(3,n-1) * x(n-1)
!
!  Divide by diagonal terms.
!
      x_new(1:n) = x_new(1:n) / a(2,1:n)
!
!  Measure norms of solution and change:
!
      x_norm = maxval ( abs ( x(1:n) ) )
      diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
!
!  Update:
!
      x(1:n) = x_new(1:n)
!
!  Test for early termination.
!
      if ( diff <= tol * ( x_norm + 1.0D+00 ) ) then
        exit
      end if

    end do

  end if

  return
end
subroutine c83_mtv ( n, a, x, b )

!*****************************************************************************80
!
!! C83_MTV multiplies a C8VEC by a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A'.
!
!    Output, complex ( kind = 8 ) B(N), the product A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) x(n)

  b(1:n)   =            a(2,1:n)   * x(1:n)
  b(1:n-1) = b(1:n-1) + a(3,1:n-1) * x(2:n) 
  b(2:n)   = b(2:n)   + a(1,2:n)   * x(1:n-1)

  return
end
subroutine c83_mv ( n, a, x, b )

!*****************************************************************************80
!
!! C83_MV multiplies a C83 matrix times a C8VEC.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, complex ( kind = 8 ) A(3,N), the matrix.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) x(n)

  b(1:n)   =            a(2,1:n)   * x(1:n)
  b(1:n-1) = b(1:n-1) + a(1,2:n)   * x(2:n)
  b(2:n)   = b(2:n)   + a(3,1:n-1) * x(1:n-1)

  return
end
subroutine c83_np_det ( n, a_lu, det )

!*****************************************************************************80
!
!! C83_NP_DET returns the determinant of a C83 system factored by C83_NP_FA.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A_LU(3,N), the LU factors computed 
!    by C83_NP_FA.
!
!    Output, complex ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a_lu(3,n)
  complex ( kind = 8 ) det

  det = product ( a_lu(2,1:n) )

  return
end
subroutine c83_np_fa ( n, a, info )

!*****************************************************************************80
!
!! C83_NP_FA factors a C83 matrix without pivoting.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output, complex ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.  On output, factorization information.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  info = 0

  do i = 1, n - 1

    if ( a(2,i) == 0.0D+00 ) then
      info = i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C83_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Store the multiplier in L.
!
    a(3,i) = a(3,i) / a(2,i)
!
!  Modify the diagonal entry in the next column.
!
    a(2,i+1) = a(2,i+1) - a(3,i) * a(1,i+1)

  end do

  if ( a(2,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C83_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop 1
  end if

  return
end
subroutine c83_np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! C83_NP_FS factors and solves a C83 system.
!
!  Discussion:
!
!    The C83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, complex ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, complex ( kind = 8 ) B(N), the right hand side.
!
!    Output, complex ( kind = 8 ) X(N), the solution.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C83_NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      stop 1
    end if
  end do

  x(1:n) = b(1:n)

  do i = 2, n
    xmult = a(3,i-1) / a(2,i-1)
    a(2,i) = a(2,i) - xmult * a(1,i)
    x(i)   = x(i)   - xmult * x(i-1)
  end do

  x(n) = x(n) / a(2,n)
  do i = n - 1, 1, -1
    x(i) = ( x(i) - a(1,i+1) * x(i+1) ) / a(2,i)
  end do

  return
end
subroutine c83_np_ml ( n, a_lu, x, b, job )

!*****************************************************************************80
!
!! C83_NP_ML computes A * x or x * A, where A has been factored by C83_NP_FA.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A(3,N), the LU factors from C83_FA.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product.
!
!    Input, integer ( kind = 4 ) JOB, specifies the product to find.
!    0, compute A * x.
!    nonzero, compute A' * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a_lu(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  complex ( kind = 8 ) x(n)

  b(1:n) = x(1:n)

  if ( job == 0 ) then
!
!  Compute X := U * X
!
    do i = 1, n

      b(i) = a_lu(2,i) * b(i)

      if ( i < n ) then
        b(i) = b(i) + a_lu(1,i+1) * b(i+1)
      end if

    end do
!
!  Compute X: = L * X.
!
    do i = n, 2, -1
      b(i) = b(i) + a_lu(3,i-1) * b(i-1)
    end do

  else
!
!  Compute X: = L' * X.
!
    do i = 1, n - 1
      b(i) = b(i) + a_lu(3,i) * b(i+1)
    end do
!
!  Compute X: = U' * X.
!
    do i = n, 2, -1
      b(i) = a_lu(2,i) * b(i)
      b(i) = b(i) + a_lu(1,i) * b(i-1)
    end do
    b(1) = a_lu(2,1) * b(1)

  end if

  return
end
subroutine c83_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! C83_NP_SL solves a C83 system factored by C83_NP_FA.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a C83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A(3,N), the factor information from C83_NP_FA.
!
!    Input/output, complex ( kind = 8 ) B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a_lu(3,n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job

  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    do i = 2, n
      b(i) = b(i) - a_lu(3,i-1) * b(i-1)
    end do
!
!  Solve U * X = Y.
!
    do i = n, 1, -1
      b(i) = b(i) / a_lu(2,i)
      if ( 1 < i ) then
        b(i-1) = b(i-1) - a_lu(1,i) * b(i)
      end if
    end do

  else
!
!  Solve U' * Y = B
!
    do i = 1, n
      b(i) = b(i) / a_lu(2,i)
      if ( i < n ) then
        b(i+1) = b(i+1) - a_lu(1,i+1) * b(i)
      end if
    end do
!
!  Solve L' * X = Y.
!
    do i = n - 1, 1, -1
      b(i) = b(i) - a_lu(3,i) * b(i+1)
    end do

  end if

  return
end
subroutine c83_print ( n, a, title )

!*****************************************************************************80
!
!! C83_PRINT prints a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  character ( len = * ) title

  call c83_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine c83_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C83_PRINT_SOME prints some of a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column, to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 3
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  character ( len = 12 ) citemp(incx)
  character ( len = 12 ) crtemp(incx)
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
  real ( kind = 8 ) xi
  real ( kind = 8 ) xr
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( crtemp(j2), '(i8,6x)' ) j
      write ( citemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col:  '',6a12)' ) ( crtemp(j2), citemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - 1 )
    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + 1 )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( 1 < i - j .or. 1 < j - i ) then

          crtemp(j2) = ' '
          citemp(j2) = ' '

        else

          if ( j == i - 1 ) then
            xr = real ( a(1,i), kind = 8 )
            xi = imag ( a(1,i) )
          else if ( j == i ) then
            xr = real ( a(2,i), kind = 8 )
            xi = imag ( a(2,i) )
          else if ( j == i + 1 ) then
            xr = real ( a(3,i), kind = 8 )
            xi = imag ( a(3,i) )
          end if

          if ( xr == 0.0D+00 .and. xi == 0.0D+00 ) then
            crtemp(j2) = '    0.0'
            citemp(j2) = ' '
          else if ( xr == 0.0D+00 .and. xi /= 0.0D+00 ) then
            crtemp(j2) = ' '
            write ( citemp(j2), '(g12.5)' ) xi
          else if ( xr /= 0.0D+00 .and. xi == 0.0D+00 ) then
            write ( crtemp(j2), '(g12.5)' ) xr
            citemp(j2) = ' '
          else
            write ( crtemp(j2), '(g12.5)' ) xr
            write ( citemp(j2), '(g12.5)' ) xi
          end if

        end if

      end do

      write ( *, '(i5,1x,6a12)' ) i, ( crtemp(j2), citemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c83_random ( n, seed, a )

!*****************************************************************************80
!
!! C83_RANDOM randomizes a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) seed

  a(1,1) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
  call c8vec_uniform_01 ( n - 1, seed, a(1,2:n) )

  call c8vec_uniform_01 ( n,     seed, a(2,1:n) )

  call c8vec_uniform_01 ( n - 1, seed, a(3,1:n-1) )
  a(3,n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c83_to_c8ge ( n, a1, a2 )

!*****************************************************************************80
!
!! C83_TO_C8GE copies a C83 matrix into a C8GE matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    The C8GE storage format is used for a complex general M by N matrix.  
!    A storage space is made for each logical entry.  The two 
!    dimensional logical array is mapped to a vector, in which storage 
!    is by columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, complex ( kind = 8 ) A1(3,N), the C83 matrix.
!
!    Output, complex ( kind = 8 ) A2(N,N), the C8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a1(3,n)
  complex ( kind = 8 ) a2(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, n

      if ( j == i - 1 ) then
        a2(i,j) = a1(1,i)
      else if ( i == j ) then
        a2(i,j) = a1(2,i)
      else if ( j == i + 1 ) then
        a2(i,j) = a1(3,i)
      else
        a2(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine c83_zeros ( n, a )

!*****************************************************************************80
!
!! C83_ZEROS zeros a C83 matrix.
!
!  Discussion:
!
!    The C83 storage format is used for a complex tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Output, complex ( kind = 8 ) A(3,N), the C83 matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(3,n)

  a(1:3,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c8ci_eval ( n, a, lambda )

!*****************************************************************************80
!
!! C8CI_EVAL returns the eigenvalues of a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis,
!    Circulant Matrices,
!    Second Edition,
!    Chelsea, 1994,
!    ISBN: 0828403384,
!    LC: QA188.D37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Output, complex ( kind = 8 ) LAMBDA(N), the eigenvalues.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  complex ( kind = 8 ) lambda(n)
  complex ( kind = 8 ) w(n)

  call c8vec_unity ( n, w )

  lambda(1:n) = a(n)
  do i = n - 1, 1, -1
    lambda(1:n) = lambda(1:n) * w(1:n) + a(i)
  end do

  call c8vec_sort_a_l2 ( n, lambda )

  return
end
subroutine c8ci_mtv ( n, a, x, b )

!*****************************************************************************80
!
!! C8CI_MTV multiplies a C8VEC by a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, complex X(N), the vector to be multiplied by A.
!
!    Output, complex B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  b(1:n) = zero

  do i = 1, n
    do j = 1, i
      b(i) = b(i) + a(i+1-j) * x(j)
    end do
    do j = i + 1, n
      b(i) = b(i) + a(n+i+1-j) * x(j)
    end do
  end do

  return
end
subroutine c8ci_mv ( n, a, x, b )

!*****************************************************************************80
!
!! C8CI_MV multiplies a C8CI matrix times a C8VEC.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  do i = 1, n
    b(i) = zero
    do j = 1, i - 1
      b(i) = b(i) + a(n+j+1-i) * x(j)
    end do
    do j = i, n
      b(i) = b(i) + a(j+1-i) * x(j)
    end do
  end do

  return
end
subroutine c8ci_print ( n, a, title )

!*****************************************************************************80
!
!! C8CI_PRINT prints a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  character ( len = * ) title

  call c8ci_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine c8ci_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8CI_PRINT_SOME prints some of a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) aij
  character ( len = 20 ) ctemp(incx)
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
  complex ( kind = 8 ) zero
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+j+1-i)
        end if

        if ( aij == zero ) then
          ctemp(j2) = '     0.0            '
        else if ( imag ( aij ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( aij, kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) aij
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8ci_random ( n, seed, a )

!*****************************************************************************80
!
!! C8CI_RANDOM randomizes a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, complex ( kind = 8 ) A(N), the C8CI matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) seed

  call c8vec_uniform_01 ( n, seed, a )

  return
end
subroutine c8ci_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! C8CI_SL solves a C8CI system.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Input, complex ( kind = 8 ) B(N), the right hand side.
!
!    Output, complex ( kind = 8 ) X(N), the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nsub
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  complex ( kind = 8 ) r5
  complex ( kind = 8 ) r6
  complex ( kind = 8 ) work(2*n-2)
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( job == 0 ) then
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = zero
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(n+2-nsub)
      r6 = a(nsub)

      if ( 2 < nsub ) then

        work(nsub-1) = r2

        do i = 1, nsub - 2
          r5 = r5 + a(n+1-i) * work(nsub-i)
          r6 = r6 + a(i+1) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( 2 < nsub ) then

        r6 = work(n)
        work(n-1+nsub-1) = zero
        do i = 2, nsub - 1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = zero
      do i = 1, nsub - 1
        r5 = r5 + a(n+1-i) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub - 1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  else
!
!  Solve the system with the principal minor of order 1.
!
    r1 = a(1)
    x(1) = b(1) / r1

    r2 = zero
!
!  Recurrent process for solving the system.
!
    do nsub = 2, n
!
!  Compute multiples of the first and last columns of
!  the inverse of the principal minor of order N.
!
      r5 = a(nsub)
      r6 = a(n+2-nsub)

      if ( 2 < nsub ) then

        work(nsub-1) = r2

        do i = 1, nsub - 2
          r5 = r5 + a(i+1) * work(nsub-i)
          r6 = r6 + a(n+1-i) * work(n-1+i)
        end do

      end if

      r2 = - r5 / r1
      r3 = - r6 / r1
      r1 = r1 + r5 * r3

      if ( 2 < nsub ) then

        r6 = work(n)
        work(n-1+nsub-1) = zero
        do i = 2, nsub - 1
          r5 = work(n-1+i)
          work(n-1+i) = work(i) * r3 + r6
          work(i) = work(i) + r6 * r2
          r6 = r5
        end do

      end if

      work(n) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
      r5 = zero
      do i = 1, nsub - 1
        r5 = r5 + a(i+1) * x(nsub-i)
      end do

      r6 = ( b(nsub) - r5 ) / r1
      do i = 1, nsub - 1
        x(i) = x(i) + work(n-1+i) * r6
      end do

      x(nsub) = r6

    end do

  end if

  return
end
subroutine c8ci_to_c8ge ( n, a, b )

!*****************************************************************************80
!
!! C8CI_TO_C8GE copies a C8CI matrix to a C8GE matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!    The C8GE storage format is used for a complex general M by N matrix.  
!    A storage space is made for each logical entry.  The two 
!    dimensional logical array is mapped to a vector, in which storage 
!    is by columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(N), the C8CI matrix.
!
!    Output, complex ( kind = 8 ) B(N,N), the C8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  complex ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, i - 1
      b(i,j) = a(n+j+1-i)
    end do
    do j = i, n
      b(i,j) = a(j+1-i)
    end do
  end do

  return
end
subroutine c8ci_zeros ( n, a )

!*****************************************************************************80
!
!! C8CI_ZEROS zeros a C8CI matrix.
!
!  Discussion:
!
!    The C8CI storage format is used for a complex N by N circulant matrix.
!    An N by N circulant matrix A has the property that the entries on
!    row I appear again on row I+1, shifted one position to the right,
!    with the final entry of row I appearing as the first of row I+1.
!    The C8CI format simply records the first row of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, complex ( kind = 8 ) A(N), the C8CI matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)

  a(1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c8ge_print ( m, n, a, title )

!*****************************************************************************80
!
!! C8GE_PRINT prints a C8GE.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call c8ge_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine c8ge_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8GE_PRINT_SOME prints some of a C8GE.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
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
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( imag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8ge_random ( m, n, seed, c )

!*****************************************************************************80
!
!! C8GE_RANDOM returns a unit pseudorandom C8GE.
!
!  Discussion:
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, double complex C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if

      theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

    end do

  end do

  return
end
subroutine c8ge_zeros ( m, n, c )

!*****************************************************************************80
!
!! C8GE_ZEROS zeros a C8GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the matrix.
!
!    Output, double complex C(M,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(m,n)

  c(1:m,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c8to_mtv ( n, a, x, b )

!*****************************************************************************80
!
!! C8TO_MTV multiplies a C8VEC by a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the entries of the first row of the
!    Toeplitz matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product A' * X.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  do i = 1, n

    b(i) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

    do j = 1, i
      b(i) = b(i) + a(i+1-j) * x(j)
    end do

    do j = i + 1, n
      b(i) = b(i) + a(n+j-i) * x(j)
    end do

  end do

  return
end
subroutine c8to_mv ( n, a, x, b )

!*****************************************************************************80
!
!! C8TO_MV multiplies a C8TO matrix times a C8VEC.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the entries of the first row of 
!    the Toeplitz matrix, followed by the entries of the first column, beginning
!    with the second row.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, complex ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  do i = 1, n

    b(i) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

    do j = 1, i - 1
      b(i) = b(i) + a(n+i-j) * x(j)
    end do

    do j = i, n
      b(i) = b(i) + a(j+1-i) * x(j)
    end do

  end do

  return
end
subroutine c8to_print ( n, a, title )

!*****************************************************************************80
!
!! C8TO_PRINT prints a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  character ( len = * ) title

  call c8to_print_some ( n, a, 1, 1, n, n, title )

  return
end
subroutine c8to_print_some ( n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8TO_PRINT_SOME prints some of a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the N by N Toeplitz matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) aij
  character ( len = 20 ) ctemp(incx)
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
  complex ( kind = 8 ) zero
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, n )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i <= j ) then
          aij = a(j+1-i)
        else
          aij = a(n+i-j)
        end if

        if ( aij == zero ) then
          ctemp(j2) = '    0.0'
        else if ( imag ( aij ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( aij, kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) aij
        end if

      end do

      write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8to_random ( n, seed, a )

!*****************************************************************************80
!
!! C8TO_RANDOM randomizes a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, complex ( kind = 8 ) A(2*N-1), the randomized matrix, with 
!    entries between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) seed

  n2 = 2 * n - 1

  call c8vec_uniform_01 ( n2, seed, a )

  return
end
subroutine c8to_sl ( n, a, b, x, job )

!*****************************************************************************80
!
!! C8TO_SL solves a C8TO system.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the first row of the Toeplitz
!    matrix, followed by the first column of the Toeplitz matrix, beginning 
!    with the second element.
!
!    Input, complex ( kind = 8 ) B(N) the right hand side vector.
!
!    Output, complex ( kind = 8 ) X(N), the solution vector.  X and B may 
!    share the same storage.
!
!    Input, integer ( kind = 4 ) JOB,
!    0 to solve A*X=B,
!    nonzero to solve A'*X=B.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n)
  complex ( kind = 8 ) c1(n-1)
  complex ( kind = 8 ) c2(n-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job
  integer ( kind = 4 ) nsub
  complex ( kind = 8 ) r1
  complex ( kind = 8 ) r2
  complex ( kind = 8 ) r3
  complex ( kind = 8 ) r5
  complex ( kind = 8 ) r6
  complex ( kind = 8 ) x(n)
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  if ( n < 1 ) then
    return
  end if
!
!  Solve the system with the principal minor of order 1.
!
  r1 = a(1)
  x(1) = b(1) / r1

  if ( n == 1 ) then
    return
  end if
!
!  Recurrent process for solving the system with the Toeplitz matrix.
!
  do nsub = 2, n
!
!  Compute multiples of the first and last columns of the inverse of
!  the principal minor of order NSUB.
!
    if ( job == 0 ) then
      r5 = a(n+nsub-1)
      r6 = a(nsub)
    else
      r5 = a(nsub)
      r6 = a(n+nsub-1)
    end if

    if ( 2 < nsub ) then

      c1(nsub-1) = r2

      do i = 1, nsub - 2
        if ( job == 0 ) then
          r5 = r5 + a(n+i) * c1(nsub-i)
          r6 = r6 + a(i+1) * c2(i)
        else
          r5 = r5 + a(i+1) * c1(nsub-i)
          r6 = r6 + a(n+i) * c2(i)
        end if
      end do

    end if

    r2 = -r5 / r1
    r3 = -r6 / r1
    r1 = r1 + r5 * r3

    if ( 2 < nsub ) then

      r6 = c2(1)
      c2(nsub-1) = zero
      do i = 2, nsub - 1
        r5 = c2(i)
        c2(i) = c1(i) * r3 + r6
        c1(i) = c1(i) + r6 * r2
        r6 = r5
      end do

    end if

    c2(1) = r3
!
!  Compute the solution of the system with the principal minor of order NSUB.
!
    r5 = zero

    do i = 1, nsub - 1
      if ( job == 0 ) then
        r5 = r5 + a(n+i) * x(nsub-i)
      else
        r5 = r5 + a(i+1) * x(nsub-i)
      end if
    end do

    r6 = ( b(nsub) - r5 ) / r1

    do i = 1, nsub - 1
      x(i) = x(i) + c2(i) * r6
    end do

    x(nsub) = r6

  end do

  return
end
subroutine c8to_to_c8ge ( n, a, b )

!*****************************************************************************80
!
!! C8TO_TO_C8GE copies a C8TO matrix to a C8GE matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!    The C8GE storage format is used for a complex general M by N matrix.  
!    A storage space is made for each logical entry.  The two 
!    dimensional logical array is mapped to a vector, in which storage 
!    is by columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, complex ( kind = 8 ) A(2*N-1), the C8TO matrix.
!
!    Output, complex ( kind = 8 ) B(N,N), the C8GE matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)
  complex ( kind = 8 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    do j = 1, i - 1
      b(i,j) = a(n+i-j)
    end do
    do j = i, n
      b(i,j) = a(j-i+1)
    end do
  end do

  return
end
subroutine c8to_zeros ( n, a )

!*****************************************************************************80
!
!! C8TO_ZEROS zeros a C8TO matrix.
!
!  Discussion:
!
!    The C8TO storage format is used for a complex Toeplitz matrix, which 
!    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
!    there are at most 2*N-1 distinct entries.  The format stores the 
!    N elements of the first row, followed by the N-1 elements of the 
!    first column (skipping the entry in the first row).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2015
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Output, complex ( kind = 8 ) A(2*N-1), the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(2*n-1)

  a(1:2*n-1) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  return
end
subroutine c8vec_indicator ( n, a )

!*****************************************************************************80
!
!! C8VEC_INDICATOR sets a C8VEC to an "indicator" vector.
!
!  Discussion:
!
!    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = cmplx ( i, -i, kind = 8 )
  end do

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! C8VEC_PRINT prints a C8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, complex ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine c8vec_print_some ( n, x, max_print, title )

!*****************************************************************************80
!
!! C8VEC_PRINT_SOME prints some of a C8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 8 ) X(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines 
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title
  complex ( kind = 8 ) x(n)

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i8,2x,2g14.6)' ) i, x(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,2g14.6)' ) i, x(i)
    end do
    i = max_print
    write ( *, '(i8,2x,2g14.6,2x,a)' ) i, x(i), '...more entries...'

  end if

  return
end
subroutine c8vec_sort_a_l2 ( n, x )

!*****************************************************************************80
!
!! C8VEC_SORT_A_L2 ascending sorts a C8VEC by L2 norm.
!
!  Discussion:
!
!    The L2 norm of A+Bi is sqrt ( A**2 + B**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input/output, complex ( kind = 8 ) X(N).
!    On input, an unsorted array.
!    On output, X has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  logical c8_le_l2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  complex ( kind = 8 ) x(n)

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      call c8_swap ( x(i), x(j) )

    else if ( indx < 0 ) then

      if ( c8_le_l2 ( x(i), x(j) ) ) then
        isgn = - 1
      else
        isgn = + 1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, double complex C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    k = seed / 127773 

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  end do

  return
end
subroutine c8vec_unity ( n, a )

!*****************************************************************************80
!
!! C8VEC_UNITY returns the N roots of unity as a C8VEC.
!
!  Discussion:
!
!    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
!
!    X(1:N)**N = ( (1,0), (1,0), ..., (1,0) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, complex ( kind = 8 ) A(N), the N roots of unity.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) theta

  do i = 1, n
    theta = r8_pi * real ( 2 * ( i - 1 ), kind = 8 ) / real ( n, kind = 8 )
    a(i) = cmplx ( cos ( theta ), sin ( theta ), kind = 8 )
  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I precedes J, ISGN = +1 if J precedes I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements 
!    I and J.  (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I precedes J;
!    ISGN => 0 means J precedes I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
