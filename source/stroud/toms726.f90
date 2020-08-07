function alga_r8 ( x )

!*****************************************************************************80
!
!! ALGA_R8 evaluates the logarithm of the gamma function.
!
!  Discussion:
!
!    A combination of recurrence and asymptotic approximation is used.
!
!    The entries in the data statement are the numerators and
!    denominators, respectively, of the quantities B[16]/(16*15),
!    B[14]/(14*13),..., B[2]/(2*1), where B[2n] are the Bernoulli
!    numbers.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the
!    Gamma Function,
!    Mathematics of Computation,
!    Volume 21, Number 98, April 1967, pages 198-203.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) ALGA_R8, the logarithm of gamma(X).
!
  implicit none

  real ( kind = 8 ) alga_r8
  real ( kind = 8 ), dimension ( 8 ) :: dbden = (/ &
    1.224D+05, &
    1.56D+02, &
    3.6036D+05, &
    1.188D+03, &
    1.68D+03, &
    1.26D+03, &
    3.6D+02, &
    1.2D+01 /)
  real ( kind = 8 ), dimension ( 8 ) :: dbnum = (/ &
    -3.617D+03, &
     1.0D+00, &
    -6.91D+02, &
     1.0D+00, &
    -1.0D+00, &
     1.0D+00, &
    -1.0D+00, &
     1.0D+00 /)
  real ( kind = 8 ) dc
  real ( kind = 8 ) dp
  real ( kind = 8 ) dprec
  real ( kind = 8 ) ds
  real ( kind = 8 ) dt
  real ( kind = 8 ) dy
  integer ( kind = 4 ) i
  real ( kind = 8 ) x
  real ( kind = 4 ) y
  real ( kind = 8 ) y0
!
!  The quantity dprec in the next statement is the number of decimal
!  digits carried in double-precision floating-point arithmetic.
!
  dprec = - log10 ( epsilon ( dbnum ) )
  dc = 0.5D+00 * log ( 8.0D+00 * atan ( 1.0D+00 ) )
  dp = 1.0D+00
  dy = x
  y = real ( dy, kind = 8 )
!
!  The quantity  y0  below is the threshold value beyond which asymptotic
!  evaluation gives sufficient accuracy; see Eq. 6.1.42 in M. Abramowitz
!  and I.A. Stegun,Handbook of Mathematical Functions''. The constants
!  are .12118868...  =  ln(10)/19 and .05390522... = ln(|B[20]|/190)/19.
!
  y0 = exp ( 0.121189D+00 * dprec + 0.053905D+00 )

  do while ( y <= y0 )
    dp = dy * dp
    dy = dy + 1.0D+00
    y = real ( dy, kind = 8 )
  end do

  dt = 1.0D+00 / ( dy * dy )
!
!  The right-hand side of the next assignment statement is B[18]/(18*17).
!
  ds = 4.3867D+04 / 2.44188D+05
  do i = 1, 8
    ds = dt * ds + dbnum(i) / dbden(i)
  end do

  alga_r8 = ( dy - 0.5D+00 ) * log ( dy ) - dy + dc + ds / dy - log ( dp )

  return
end
subroutine cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, ierr )

!*****************************************************************************80
!
!! CHEB_R8 implements the modified Chebyshev algorithm.
!
!  Discussion:
!
!    The routine generates recursion coefficients ALPHA and BETA.
!
!    Given a set of polynomials  p(0),p(1),...,p(2*n-1)  satisfying
!
!      p(k+1)(x) = (x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
!      k = 0,1,...,2*n-2,
!
!      p(-1)(x) = 0,  p(0)(x)=1,
!
!    and associated modified moments
!
!      fnu(k) = integral of p(k)(x)*dlambda(x),
!      k = 0,1,...,2*n-1,
!
!    this subroutine uses the modified Chebyshev algorithm to generate the
!    recursion coefficients  alpha(k),beta(k), k = 0,1,...,n-1, for the
!    polynomials  pi(k)  orthogonal with respect to the integration
!    measure  dlambda(x), i.e.,
!
!      pi(k+1)(x) = (x-alpha(k))*pi(k)(x)-beta(k)*pi(k-1)(x),
!      k = 0,1,...,n-1,
!
!      pi(-1)(x) = 0,  pi(0)(x)=1.
!
!    On machines with limited exponent range, the occurrence of underflow
!    [overflow] in the computation of the  alpha's  and  beta's  can often
!    be avoided by multiplying all modified moments by a sufficiently large
!    [small] scaling factor and dividing the new  beta(0)  by the same
!    scaling factor.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of recursion coefficients 
!    desired.
!
!    Input, real ( kind = 8 ) A(2*N-1), B(2*N-1), the values of A(k-1), b(k-1),
!    k = 1,2,...,2*n-1.
!
!    Input, real ( kind = 8 ) FNU(2*N), the values of the modified moments 
!    fnu(k-1), k = 1,2,...,2*n
!
!    Output, real ( kind = 8 ) ALPHA(N), BETA(N), the
!    recursion coefficients  alpha(k-1),beta(k-1),
!    k = 1,2,...,n, where  beta(0)  is the total mass.
!
!    Output, real ( kind = 8 ) S(N), the normalization factors
!    s(k) = integral [pi(k)(x)]**2 dlambda(x), k=0,1,
!    2,...,n-1.
!
!    Output, integer ( kind = 4 ) IERR, an error flag.
!    0, on normal return,
!    1, if  abs(fnu(0))  is less than the machine zero,
!    2, if  n  is out of range,
!    -k  if S(K), k = 0,1,2,...,n-1, is about to underflow,
!    +k  if S(K) is about to overflow.
!
  implicit none

  integer n

  real ( kind = 8 ) da(2*n-1)
  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db(2*n-1)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dnu(2*n)
  real ( kind = 8 ) ds(n)
  real ( kind = 8 ) ds0(2*n)
  real ( kind = 8 ) ds1(2*n)
  real ( kind = 8 ) ds2(2*n)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lk

  ierr = 0

  if ( abs ( dnu(1) ) < 10.0D+00 * tiny ( dnu(1) ) ) then
    ierr = 1
    return
  end if

  if ( n < 1 ) then
    ierr = 2
    return
  end if

  dalpha(1) = da(1) + dnu(2) / dnu(1)
  dbeta(1) = dnu(1)

  if ( n == 1 ) then
    return
  end if

  ds(1) = dnu(1)
  ds0(1:2*n) = 0.0D+00
  ds1(1:2*n) = dnu(1:2*n)

  do k = 2, n

    lk = 2 * n - k + 1

    do l = k, 2 * n - k + 1

      ds2(l) = ds1(l+1) - ( dalpha(k-1) - da(l) ) * ds1(l) &
        - dbeta(k-1) * ds0(l) + db(l) * ds1(l-1)

      if ( l == k ) then
        ds(k) = ds2(k)
      end if

    end do

    if ( abs ( ds(k) ) < 10.0D+00 * tiny ( ds(k) ) ) then
      ierr = -( k - 1 )
      return
    end if

    if ( 0.1D+00 * huge ( ds(k) ) < abs ( ds(k) ) ) then
      ierr = k - 1
      return
    end if

    dalpha(k) = da(k) + ( ds2(k+1) / ds2(k) ) - ( ds1(k) / ds1(k-1) )
    dbeta(k) = ds2(k) / ds1(k-1)

    ds0(k:lk) = ds1(k:lk)
    ds1(k:lk) = ds2(k:lk)

  end do

  return
end
subroutine chri_r8 ( n, iopt, da, db, dx, dy, dhr, dhi, dalpha, dbeta, ierr )

!*****************************************************************************80
!
!! CHRI_R8 implements the Christoffel or generalized Christoffel theorem.
!
!  Discussion:
!
!    In all cases except  iopt = 7, it uses nonlinear recurrence
!    algorithms described in W. Gautschi,An algorithmic implementation
!    of the generalized Christoffel theorem'', Numerical Integration
!    (G. Haemmerlin, ed.), Birkhaeuser, Basel, 1982, pp. 89-106. The case
!    iopt = 7  incorporates a QR step with shift  x  in the manner of
!    J. Kautsky and G.H. Golub, On the calculation of Jacobi matrices'',
!    Linear Algebra Appl. 52/53, 1983, 439-455, using the algorithm of
!    Eq. (67.11) on p. 567 in J.H. Wilkinson,The Algebraic Eigenvalue
!    Problem'', Clarendon Press, Oxford, 1965. Given the recursion
!    coefficients  a(k),b(k), k = 0,1,...,n, for the (monic) orthogonal
!    polynomials with respect to some measure  dlambda(t), it generates
!    the recursion coefficients  alpha(k),beta(k), k = 0,1,...,n-1, for the
!    measure
!
!      (t-x)dlambda(t)               if  iopt = 1
!      [(t-x)**2+y**2]dlambda(t)     if  iopt = 2
!      (t**2+y**2)dlambda(t) with    if  iopt = 3
!      dlambda(t) and supp(dlambda) symmetric  with respect to
!      the origin
!      dlambda(t)/(t-x)              if  iopt = 4
!      dlambda(t)/[(t-x)**2+y**2]    if  iopt = 5
!      dlambda(t)/(t**2+y**2) with   if  iopt = 6
!      dlambda(t) and supp(dlambda) symmetric with respect to
!      the origin
!      [(t-x)**2]dlambda(t)          if  iopt = 7
!
!    It is assumed that  n  is larger than or equal to 2. Otherwise, the
!    routine exits immediately with the error flag  ierr  set equal to 1.
!    If  iopt  is not between 1 and 7, the routine exits with  ierr = 2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of recurrence coefficients 
!    desired.
!
!    Input, integer ( kind = 4 ) IOPT, the desired weight distribution.
!
!    Input, real ( kind = 8 ) A(N+1), B(N+1), the recursion coefficients
!    a(k-1),b(k-1),k = 1,2,...,n+1, of the polynomials orthogonal with
!    respect to the given measure  dlambda(t)
!
!    Input, real ( kind = 8 ) X, Y, the linear and quadratic factors, 
!    or divisors, of dlambda(t).
!
!    Input, real ( kind = 8 ) HR, HI, the real and imaginary part, 
!    respectively, of the integral of dlambda(t)/(z-t), where z = x+iy;
!    the parameter  hr  is used only if  iopt = 4 or
!    5, the parameter  hi  only if  iopt = 5 or 6
!
!    Output, real ( kind = 8 ) ALPHA(N), BETA(N), the desired recursion 
!    coefficients alpha(k-1), beta(k-1), k = 1,2,...,n
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) IERR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) da(n+1)
  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db(n+1)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dc
  real ( kind = 8 ) dc0
  real ( kind = 8 ) dcm1
  real ( kind = 8 ) dd
  real ( kind = 8 ) de
  real ( kind = 8 ) dei
  real ( kind = 8 ) deio
  real ( kind = 8 ) deioo
  real ( kind = 8 ) deo
  real ( kind = 8 ) deoo
  real ( kind = 8 ) deps
  real ( kind = 8 ) der
  real ( kind = 8 ) dero
  real ( kind = 8 ) deroo
  real ( kind = 8 ) dgam
  real ( kind = 8 ) dhi
  real ( kind = 8 ) dhr
  real ( kind = 8 ) dp2
  real ( kind = 8 ) dq
  real ( kind = 8 ) ds
  real ( kind = 8 ) dso
  real ( kind = 8 ) dt
  real ( kind = 8 ) du
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) k

  deps = 5.0D+00 * epsilon ( deps )
  ierr = 0

  if ( n < 2 ) then
    ierr = 1
    return
  end if

  if ( iopt == 1 ) then

    de = 0.0D+00
    do k = 1, n
      dq = da(k) - de - dx
      dbeta(k) = dq * de
      de = db(k+1) / dq
      dalpha(k) = dx + dq + de
    end do

    dbeta(1) = db(1) * ( da(1) - dx )

  else if ( iopt == 2 ) then

    ds = dx - da(1)
    dt = dy
    deio = 0.0D+00

    do k = 1, n
      dd = ds * ds + dt * dt
      der = - db(k+1) * ds / dd
      dei = db(k+1) * dt / dd
      ds = dx + der - da(k+1)
      dt = dy + dei
      dalpha(k) = dx + dt * der / dei - ds * dei / dt
      dbeta(k) = dt * deio * ( 1.0D+00 + ( der / dei )**2 )
      deio = dei
    end do

    dbeta(1) = db(1) * ( db(2) + ( da(1) - dx )**2 + dy * dy )

  else if ( iopt == 3 ) then

    dt = dy
    deio = 0.0D+00

    do k = 1, n
      dei = db(k+1) / dt
      dt = dy + dei
      dalpha(k) = 0.0D+00
      dbeta(k) = dt * deio
      deio = dei
    end do

    dbeta(1) = db(1) * ( db(2) + dy * dy )

  else if ( iopt == 4 ) then

    dalpha(1) = dx - db(1) / dhr
    dbeta(1) = - dhr
    dq = - db(1) / dhr

    do k = 2, n
      de = da(k-1) - dx - dq
      dbeta(k) = dq * de
      dq = db(k) / de
      dalpha(k) = dq + de + dx
    end do

  else if ( iopt == 5 ) then

    dd = dhr * dhr + dhi * dhi
    deroo = da(1) - dx + db(1) * dhr / dd
    deioo = - db(1) * dhi / dd - dy
    dalpha(1) = dx + dhr * dy / dhi
    dbeta(1) = - dhi / dy
    dalpha(2) = dx - db(1) * dhi * deroo / ( dd * deioo ) + dhr * deioo / dhi
    dbeta(2) = dy * deioo * ( 1.0D+00 + ( dhr / dhi )**2 )

    if ( n == 2 ) then
      return
    end if

    dso = db(2) / ( deroo**2 + deioo**2 )
    dero = da(2) - dx - dso * deroo
    deio = dso * deioo - dy
    dalpha(3) = dx + deroo * deio / deioo + dso * deioo * dero / deio
    dbeta(3) = - db(1) * dhi * deio * ( 1.0D+00 + ( deroo / deioo )**2 ) / dd

    do k = 3, n - 1
      ds = db(k) / ( dero**2 + deio**2 )
      der = da(k) - dx - ds * dero
      dei = ds * deio - dy
      dalpha(k+1) = dx + dero * dei / deio + ds * deio * der / dei
      dbeta(k+1) = dso * deioo * dei * ( 1.0D+00 + ( dero / deio )**2 )
      deroo = dero
      deioo = deio
      dero = der
      deio = dei
      dso = ds
    end do

  else if ( iopt == 6 ) then

    deoo = - db(1) / dhi - dy
    deo = db(2) / deoo - dy
    dalpha(1) = 0.0D+00
    dbeta(1) = - dhi / dy
    dalpha(2) = 0.0D+00
    dbeta(2) = dy * deoo

    if ( n == 2 ) then
      return
    end if

    dalpha(3) = 0.0D+00
    dbeta(3) = - db(1) * deo / dhi

    do k = 3, n - 1
      de = db(k) / deo - dy
      dbeta(k+1) = db(k-1) * de / deoo
      dalpha(k+1) = 0.0D+00
      deoo = deo
      deo = de
    end do

  else if ( iopt == 7 ) then

    du = 0.0D+00
    dc = 1.0D+00
    dc0 = 0.0D+00

    do k = 1, n

      dgam = da(k) - dx - du
      dcm1 = dc0
      dc0 = dc

      if ( deps < abs ( dc0 ) ) then
        dp2 = ( dgam**2 ) / dc0
      else
        dp2 = dcm1*db(k)
      end if

      if ( 1 < k ) then
        dbeta(k) = ds * ( dp2 + db(k+1) )
      end if

      ds = db(k+1) / ( dp2 + db(k+1) )
      dc = dp2 / ( dp2 + db(k+1) )
      du = ds * ( dgam + da(k+1) - dx )
      dalpha(k) = dgam + du + dx

    end do

    dbeta(1) = db(1) * ( db(2) + ( dx - da(1) )**2 )

  else

    ierr = 2

  end if

  return
end
subroutine fejer_r8 ( n, x, w )

!*****************************************************************************80
!
!! FEJER_R8 generates a Fejer quadrature rule.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of quadrature nodes.
!
!    Output, real ( kind = 8 ) X(N), W(N), the quadrature nodes and weights.
!    The nodes are listed in increasing order.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) np1h
  real ( kind = 8 ) pi
  real ( kind = 8 ) t
  real ( kind = 8 ) total
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  pi = 4.0D+00 * atan ( 1.0D+00 )
  nh = n / 2
  np1h = ( n + 1 ) / 2
  do k = 1, nh
    x(n+1-k) = cos ( 0.5D+00 * real ( 2 * k - 1, kind = 8 ) * pi &
      / real ( n, kind = 8 ) )
    x(k) = -x(n+1-k)
  end do

  if ( 2 * nh /= n ) then
    x(np1h) = 0.0D+00
  end if

  do k = 1, np1h

    c1 = 1.0D+00
    c0 = 2.0D+00 * x(k) * x(k) - 1.0D+00
    t = 2.0D+00 * c0
    total = c0 / 3.0D+00

    do m = 2, nh
      c2 = c1
      c1 = c0
      c0 = t * c1 - c2
      total = total + c0 / real ( 4 * m * m - 1, kind = 8 )
    end do

    w(k) = 2.0D+00 * ( 1.0D+00 - 2.0D+00 * total ) / real ( n, kind = 8 )
    w(n+1-k) = w(k)

  end do

  return
end
function gamma_r8 ( x, ierr )

!*****************************************************************************80
!
!! GAMMA_R8 evaluates the gamma function for real positive argument.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IERR, an error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Output, real ( kind = 8 ) GAMMA_R4, the value of the Gamma
!   function at X.
!
  implicit none

  real ( kind = 8 ) alga_r8
  real ( kind = 8 ) almach
  real ( kind = 8 ) gamma_r8
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  almach = log ( huge ( almach ) )
  ierr = 0
  t = alga_r8 ( x )

  if ( almach <= t ) then
    ierr = 2
    gamma_r8 = huge ( gamma_r8 )
  else
    gamma_r8 = exp ( t )
  end if

  return
end
subroutine gauss_r8 ( n, dalpha, dbeta, deps, dzero, dweigh, ierr, de )

!*****************************************************************************80
!
!! GAUSS_R8 generates an N-point Gaussian quadrature formula.
!
!  Discussion:
!
!    Given N and a measure DLAMBDA, this routine generates the N-point
!    Gaussian quadrature formula:
!
!      Integral over supp(dlambda) of f(x) dlambda(x)
!
!        =  sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k) and the weights as
!    weight(k) = w(k), k=1,2,...,n. The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,2,...,n-1, for the measure
!    dlambda. The routine computes the nodes as eigenvalues, and the
!    weights in term of the first component of the respective normalized
!    eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
!    It uses a translation and adaptation of the algol procedure  imtql2,
!    Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
!    by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
!    Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
!    routine  imtql2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the Gaussian 
!    quadrature formula.
!
!    Input, real ( kind = 8 ) ALPHA(N), BETA(N), the values of alpha(k-1), 
!    beta(k-1), k = 1,2,...,n
!
!    Input, real ( kind = 8 ) EPS, the relative accuracy desired in the nodes 
!    and weights.
!
!    Output, real ( kind = 8 ) ZERO(N), the Gaussian nodes in increasing order.
!
!    Output, real ( kind = 8 ) WEIGHT(N), the Gaussian weights.
!
!    Output, integer ( kind = 4 ) IERR, an error flag equal to 0 on normal 
!    return, equal to I if the QR algorithm does not converge within 30
!    iterations on evaluating the i-th eigenvalue, equal to -1 if N is not
!    in range, and equal to -2 if one of the BETA's is negative.
!
!    Workspace, real ( kind = 8 ) DE(N).
!
  implicit  none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dc
  real ( kind = 8 ) de(n)
  real ( kind = 8 ) deps
  real ( kind = 8 ) df
  real ( kind = 8 ) dg
  real ( kind = 8 ) dp
  real ( kind = 8 ) dr
  real ( kind = 8 ) ds
  real ( kind = 8 ) dweigh(n)
  real ( kind = 8 ) dzero(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) mml

  if ( n < 1 ) then
    ierr = -1
    return
  end if

  ierr = 0

  dzero(1) = dalpha(1)

  if ( dbeta(1) < 0.0D+00 ) then
    ierr = -2
    return
  end if

  dweigh(1) = dbeta(1)

  if ( n == 1 ) then
    return
  end if

  dweigh(1) = 1.0D+00
  de(n) = 0.0D+00

  do k = 2, n

    dzero(k) = dalpha(k)

    if ( dbeta(k) < 0.0D+00 ) then
      ierr = -2
      return
    end if

    de(k-1) = sqrt ( dbeta(k) )
    dweigh(k) = 0.0D+00

  end do

  do l = 1, n

    j = 0

    do

      do m2 = l, n

        m = m2

        if ( m == n ) then
          exit
        end if

        if ( abs ( de(m) ) <= deps &
          * ( abs ( dzero(m) ) + abs ( dzero(m+1) ) ) ) then
          exit
        end if

      end do

      dp = dzero(l)

      if ( m == l ) then
        exit
      end if

      if ( 30 <= j ) then
        ierr = 1
        return
      end if

      j = j + 1
      dg = ( dzero(l+1) - dp ) / ( 2.0D+00 * de(l) )
      dr = sqrt ( dg * dg + 1.0D+00 )
      dg = dzero(m) - dp + de(l) / ( dg + sign ( dr, dg ) )
      ds = 1.0D+00
      dc = 1.0D+00
      dp = 0.0D+00
      mml = m - l

      do i = m - 1, 1, -1

        df = ds * de(i)
        db = dc * de(i)

        if ( abs ( df ) < abs ( dg ) ) then

          ds = df / dg
          dr = sqrt ( ds * ds + 1.0D+00 )
          de(i+1) = dg * dr
          dc = 1.0D+00 / dr
          ds = ds * dc

        else

          dc = dg / df
          dr = sqrt ( dc * dc + 1.0D+00 )
          de(i+1) = df * dr
          ds = 1.0D+00 / dr
          dc = dc * ds

        end if

        dg = dzero(i+1) - dp
        dr = ( dzero(i) - dg ) * ds + 2.0D+00 * dc * db
        dp = ds * dr
        dzero(i+1) = dg + dp
        dg = dc * dr - db
        df = dweigh(i+1)
        dweigh(i+1) = ds * dweigh(i) + dc * df
        dweigh(i) = dc * dweigh(i) - ds * df

      end do

      dzero(l) = dzero(l) - dp
      de(l) = dg
      de(m) = 0.0D+00

    end do

  end do

  do i = 1, n - 1

    k = i
    dp = dzero(i)

    do j = i + 1, n
      if ( dzero(j) < dp ) then
        k = j
        dp = dzero(j)
      end if
    end do

    if ( k /= i ) then
      dzero(k) = dzero(i)
      dzero(i) = dp
      dp = dweigh(i)
      dweigh(i) = dweigh(k)
      dweigh(k) = dp
    end if

  end do

  dweigh(1:n) = dbeta(1) * dweigh(1:n)**2

  return
end
subroutine gchri_r8 ( n, iopt, nu0, numax, deps, da, db, dx, dy, dalpha, &
  dbeta, nu, ierr, ierrc, drhor, drhoi, droldr, droldi )

!*****************************************************************************80
!
!! GCHRI_R8 implements the generalized Christoffel theorem.
!
!  Discussion:
!
!    The routine uses the method of modified moments.
!
!    Given the recursion coefficients  a(k), b(k), k = 0,1,...n, for the monic
!    orthogonal polynomials with respect to some measure  dlambda(t), it
!    generates the recursion coefficients  alpha(k), beta(k), k = 0,1,2,...,n-1
!    for the measure
!
!      dlambda(t)/(t-x)        if iopt = 1
!      dlambda(t)/{(t-x)**2+y**2} if iopt = 2
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!    Walter Gautschi, 
!    Minimal solutions of three-term recurrence relations and orthogonal 
!    polynomials, 
!    Mathematics of Computation,
!    Volume 36, 1981, pages 547-554.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of recurrence coefficients 
!    desired.
!
!    Input, integer ( kind = 4 ) IOPT, selects the weight distribution.
!
!    Input, integer ( kind = 4 ) NU0, estimates the starting backward 
!    recurrence index.  If no good estimate is available, set NU0 = 3 * N.
!
!    Input, integer ( kind = 4 ) NUMAX, controls the termination of the backward
!    recursion in case of nonconvergence; a conservative choice is
!    NUMAX = 500.
!
!    Input, real ( kind = 8 ) DEPS, a relative error tolerance.
!
!    Input, real ( kind = 8 ) DA(NUMAX), DB(NUMAX), the recursion coefficients 
!    A(k) = alpha(k-1), B(k) = beta(k), k = 1,2,...,numax, for the measure  
!    dlambda
!
!    Input, real ( kind = 8 ) DX, DY, define the linear and quadratic
!    divisors of dlambda
!
!   Output: 
!
!    alpha,beta - arrays of dimension  n  containing the desired
!                   recursion coefficients  alpha(k-1), beta(k-1), k = 1,
!                   2,...,n
!
!           nu  - - the backward recurrence index yielding convergence;
!                   in case of nonconvergence,  nu  will have the value
!                   numax
!
!           ierr  - an error flag, where
!                   ierr = 0     on normal return
!                   ierr = 1     if  iopt  is neither 1 nor 2
!                   ierr = nu0   if  numax < nu0.
!                   ierr = numax if the backward recurrence algorithm does
!                              not converge
!                   ierr = -1    if  n  is not in range
!
!    Output, integer ( kind = 4 ) IERRC, an error flag from CHEB_R8.
!    0, no error was detected.
!    nonzero, an error was detected.
!
!    DRHOR, DRHOI, DROLDR, DROLDI, ?
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) numax

  real ( kind = 8 ) da(numax)
  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db(numax)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) deps
  real ( kind = 8 ) dnu(2*n)
  real ( kind = 8 ) drhoi(2*n)
  real ( kind = 8 ) drhor(2*n)
  real ( kind = 8 ) droldi
  real ( kind = 8 ) droldr
  real ( kind = 8 ) ds(n)
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ierrc
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) nu
  integer ( kind = 4 ) nu0

  if ( n < 1 ) then
    ierr = -1
    return
  end if

  ierr = 0

  if ( iopt == 1 ) then

    call knum_r8 ( 2 * n - 1, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, &
      nu, ierr )

    dnu(1:2*n) = - drhor(1:2*n)

    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, ierrc )

  else if ( iopt == 2 ) then

    dy = abs ( dy )

    call knum_r8 ( 2 * n - 1, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, &
      nu, ierr )

    dnu(1:2*n) = - drhoi(1:2*n) / dy

    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, ierrc )

  else

    ierr = 1

  end if

  return
end
subroutine kern_r8 ( n, nu0, numax, dx, dy, deps, da, db, dkerr, dkeri, nu, &
  ierr )

!*****************************************************************************80
!
!! KERN_R8 generates the kernels in the Gauss quadrature remainder term.
!
!  Discussion:
!
!    This routine was written under the assumption that double precision complex
!    arithmetic was NOT available.
!
!    The kernels to be generated are
!
!      K(k)(z) = rho(k)(z)/pi(k)(z), k=0,1,2,...,n,
!
!    where rho(k) are the output quantities of the routine KNUM_R8, and
!    pi(k) the (monic) orthogonal polynomials.  The results are returned
!    in the array ker as ker(k) = K(k-1)(z), k=1,2,...,n+1.  All the other
!    input and output parameters have the same meaning as in the routine
!    KNUM_R8.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, indicates the number of kernels to
!    be computed.
!
!    Output, integer ( kind = 4 ) IERR, an error flag from KNUM_R8.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) numax

  real ( kind = 8 ) da(numax)
  real ( kind = 8 ) db(numax)
  real ( kind = 8 ) dden
  real ( kind = 8 ) deps
  real ( kind = 8 ) dkeri(n+1)
  real ( kind = 8 ) dkerr(n+1)
  real ( kind = 8 ) dp0i
  real ( kind = 8 ) dp0r
  real ( kind = 8 ) dpi
  real ( kind = 8 ) dpm1i
  real ( kind = 8 ) dpm1r
  real ( kind = 8 ) dpr
  real ( kind = 8 ) dt
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nu
  integer ( kind = 4 ) nu0

  call knum_r8 ( n, nu0, numax, dx, dy, deps, da, db, dkerr, dkeri, nu, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  dp0r = 0.0D+00
  dp0i = 0.0D+00
  dpr = 1.0D+00
  dpi = 0.0D+00

  do k = 1, n

    dpm1r = dp0r
    dpm1i = dp0i
    dp0r = dpr
    dp0i = dpi

    dpr = ( dx - da(k) ) * dp0r - dy * dp0i - db(k) * dpm1r
    dpi = ( dx - da(k) ) * dp0i + dy * dp0r - db(k) * dpm1i

    dden = dpr**2 + dpi**2
    dt = ( dkerr(k+1) * dpr + dkeri(k+1) * dpi ) / dden

    dkeri(k+1) = ( dkeri(k+1) * dpr - dkerr(k+1) * dpi ) / dden
    dkerr(k+1) = dt

  end do

  return
end
subroutine knum_r8 ( n, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, nu, &
  ierr )

!*****************************************************************************80
!
!! KNUM_R8 integrates certain rational polynomials.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, indicates that RHO(0:N) are desired.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) numax

  real ( kind = 8 ) da(numax)
  real ( kind = 8 ) db(numax)
  real ( kind = 8 ) dden
  real ( kind = 8 ) deps
  logical done
  real ( kind = 8 ) drhoi(n+1)
  real ( kind = 8 ) drhor(n+1)
  real ( kind = 8 ) droldi(n+1)
  real ( kind = 8 ) droldr(n+1)
  real ( kind = 8 ) dri
  real ( kind = 8 ) drr
  real ( kind = 8 ) dt
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nu
  integer ( kind = 4 ) nu0

  ierr = 0

  if ( numax < nu0 ) then
    ierr = nu0
    return
  end if

  nu0 = max ( nu0, n + 1 )

  nu = nu0 - 5

  drhor(1:n+1) = 0.0D+00
  drhoi(1:n+1) = 0.0D+00

  do

    nu = nu + 5

    if ( numax < nu ) then
      ierr = numax
      exit
    end if

    droldr(1:n+1) = drhor(1:n+1)
    droldi(1:n+1) = drhoi(1:n+1)

    drr = 0.0D+00
    dri = 0.0D+00

    do j = 1, nu

      j1 = nu - j + 1
      dden = ( dx - da(j1) - drr )**2 + ( dy - dri )**2
      drr = db(j1) * ( dx - da(j1) - drr ) / dden
      dri = - db(j1) * ( dy - dri ) / dden

      if ( j1 <= n + 1 ) then
        drhor(j1) = drr
        drhoi(j1) = dri
      end if

    end do

    done = .true.

    do k = 1, n + 1

      if ( ( deps**2 ) * ( drhor(k)**2 + drhoi(k)**2 ) < &
        ( drhor(k) - droldr(k) )**2 + ( drhoi(k) - droldi(k) )**2 ) then
        done = .false.
        exit
      end if

    end do

    if ( done ) then
      exit
    end if

  end do

  do k = 2, n + 1
    dt = drhor(k) * drhor(k-1) - drhoi(k) * drhoi(k-1)
    drhoi(k) = drhor(k) * drhoi(k-1) + drhoi(k) * drhor(k-1)
    drhor(k) = dt
  end do

  return
end
subroutine lancz_r8 ( n, ncap, dx, dw, dalpha, dbeta, ierr, dp0, dp1 )

!*****************************************************************************80
!
!! LANCZ_R8 applies Stieltjes's procedure, using the Lanczos method.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncap

  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dgam
  real ( kind = 8 ) dp0(ncap)
  real ( kind = 8 ) dp1(ncap)
  real ( kind = 8 ) dpi
  real ( kind = 8 ) drho
  real ( kind = 8 ) dsig
  real ( kind = 8 ) dt
  real ( kind = 8 ) dtk
  real ( kind = 8 ) dtmp
  real ( kind = 8 ) dtsig
  real ( kind = 8 ) dw(ncap)
  real ( kind = 8 ) dx(ncap)
  real ( kind = 8 ) dxlam
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k

  if ( n <= 0 .or. ncap < n ) then
    ierr = 1
    return
  end if

  ierr = 0

  dp0(1:ncap) = dx(1:ncap)
  dp1(1) = dw(1)
  dp1(2:ncap) = 0.0D+00

  do i = 1, ncap - 1

    dpi = dw(i+1)
    dgam = 1.0D+00
    dsig = 0.0D+00
    dt = 0.0D+00
    dxlam = dx(i+1)

    do k = 1, i + 1

      drho = dp1(k) + dpi
      dtmp = dgam * drho
      dtsig = dsig

      if ( drho <= 0.0D+00 ) then
        dgam = 1.0D+00
        dsig = 0.0D+00
      else
        dgam = dp1(k) / drho
        dsig = dpi / drho
      end if

      dtk = dsig * ( dp0(k) - dxlam ) - dgam * dt
      dp0(k) = dp0(k) - ( dtk - dt )
      dt = dtk

      if ( dsig <= 0.0D+00 ) then
        dpi = dtsig * dp1(k)
      else
        dpi = ( dt**2 ) / dsig
      end if

      dtsig = dsig
      dp1(k) = dtmp

    end do
  end do

  dalpha(1:n) = dp0(1:n)
  dbeta(1:n) = dp1(1:n)

  return
end
subroutine lob_r8 ( n, alpha, beta, left, right, zero, weight, ierr )

!*****************************************************************************80
!
!! LOB_R8 generates a Gauss-Lobatto quadrature rule.
!
!  Discussion:
!
!    Given N and a measure DLAMBDA, this routine generates the
!    (N+2)-point Gauss-Lobatto quadrature formula
!
!      Integral over support(DLAMBDA) of F(X) DLAMBDA(X)
!
!      = w(0) f(x(0)) + sum from k=1 to k=n of w(k) f(x(k))
!      + w(n+1) f(x(n+1)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k), the weights as  weight(k)
!    = w(k), k=0,1,...,n,n+1.  The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,...,n,n+1, for the measure
!    dlambda.  The nodes and weights are computed in terms of the
!    eigenvalues and first component of the normalized eigenvectors of
!    a slightly modified Jacobi matrix of order  n+2.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of interior points in the 
!    Gauss-Lobatto formula.
!
!    Input, real ( kind = 8 ) ALPHA(N+2), BETA(N+2), the recursion coefficients
!    alpha(k-1), beta(k-1), k = 1,2,...,n+2, of the underlying measure;
!    the routine does not use alpha(n+2), beta(n+2).
!
!    Input, real ( kind = 8 ) LEFT, RIGHT, the prescribed left and right
!    endpoints x(0) and x(n+1) of the Gauss-Lobatto formula
!
!    Output, real ( kind = 8 ) ZERO(N+2), the nodes in increasing order.  
!    zero(k) = x(k), k=0,1,...,n,n+1
!
!    Output, real ( kind = 8 ) WEIGHT(N+2), the weights
!    weight(k) = w(k), k=0,1,...,n,n+1
!
!    Output, integer ( kind = 4 ) IERR, an error flag from the routine GAUSS_R8.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n+2)
  real ( kind = 8 ) alpha(n+2)
  real ( kind = 8 ) b(n+2)
  real ( kind = 8 ) beta(n+2)
  real ( kind = 8 ) det
  real ( kind = 8 ) e(n+2)
  real ( kind = 8 ) epsma
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  real ( kind = 8 ) left
  real ( kind = 8 ) p0l
  real ( kind = 8 ) p0r
  real ( kind = 8 ) p1l
  real ( kind = 8 ) p1r
  real ( kind = 8 ) pm1l
  real ( kind = 8 ) pm1r
  real ( kind = 8 ) right
  real ( kind = 8 ) weight(n+2)
  real ( kind = 8 ) zero(n+2)

  epsma = epsilon ( epsma )

  a(1:n+1) = alpha(1:n+1)
  b(1:n+1) = beta(1:n+1)

  p0l = 0.0D+00
  p0r = 0.0D+00
  p1l = 1.0D+00
  p1r = 1.0D+00

  do k = 1, n + 1
    pm1l = p0l
    p0l = p1l
    pm1r = p0r
    p0r = p1r
    p1l = ( left - a(k) ) * p0l - b(k) * pm1l
    p1r = ( right - a(k) ) * p0r - b(k) * pm1r
  end do

  det = p1l * p0r - p1r * p0l
  a(n+2) = ( left * p1l * p0r - right * p1r * p0l ) / det
  b(n+2) = ( right - left ) * p1l * p1r / det

  call gauss_r8 ( n + 2, a, b, epsma, zero, weight, ierr, e )

  return
end
subroutine mccheb_r8 ( n, ncapm, mc, mp, dxp, dyp, quad_r8, deps, iq, &
  idelta, finld, finrd, dendl, dendr, dxfer, dwfer, da, db, dnu, dalpha, &
  dbeta, ncap, kount, ierrd )

!*****************************************************************************80
!
!! MCCHEB_R8 is a multiple-component discretized modified Chebyshev algorithm.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERRD, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer ( kind = 4 ) mc
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncapm

  real ( kind = 8 ) da(2*n-1)
  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) db(2*n-1)
  real ( kind = 8 ) dbe(n)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dendl(mc)
  real ( kind = 8 ) dendr(mc)
  real ( kind = 8 ) deps
  real ( kind = 8 ) dnu(2*n)
  logical done
  real ( kind = 8 ) dp
  real ( kind = 8 ) dp1
  real ( kind = 8 ) dpm1
  real ( kind = 8 ) ds(n)
  real ( kind = 8 ) dsum
  real ( kind = 8 ) dw(ncapm)
  real ( kind = 8 ) dwfer(ncapm)
  real ( kind = 8 ) dwm(mc*ncapm+mp)
  real ( kind = 8 ) dx(ncapm)
  real ( kind = 8 ) dxfer(ncapm)
  real ( kind = 8 ) dxm(mc*ncapm+mp)
  real ( kind = 8 ) dxp(mp)
  real ( kind = 8 ) dyp(mp)
  logical finld
  logical finrd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idelta
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ierrd
  integer ( kind = 4 ) im1tn
  integer ( kind = 4 ) incap
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mcd
  integer ( kind = 4 ) mtncap
  integer ( kind = 4 ) mtnpmp
  integer ( kind = 4 ) ncap
  integer ( kind = 4 ) nd
  external quad_r8

  nd = 2 * n

  if ( idelta <= 0 ) then
    idelta = 1
  end if

  if ( n < 1 ) then
    ierrd = -1
    return
  end if

  incap = 1
  kount = -1
  ierrd = 0
  dbeta(1:n) = 0.0D+00
  ncap = ( 2 * n - 1 ) / idelta

  do

    dbe(1:n) = dbeta(1:n)

    kount = kount + 1

    if ( 1 < kount ) then
      incap = 2**( kount / 5 ) * n
    end if

    ncap = ncap + incap

    if ( ncapm < ncap ) then
      ierrd = ncapm
      return
    end if

    mtncap = mcd * ncap

    do i = 1, mc

      im1tn = ( i - 1 ) * ncap
      if ( iq == 1 ) then
        call quad_r8 ( ncap, dx, dw, i, ierr )
      else
        call qgp_r8 ( ncap, dx, dw, i, ierr, mc, finld, finrd, dendl, dendr, &
          dxfer, dwfer )
      end if

      if ( ierr /= 0 ) then
        ierrd = i
        return
      end if

      dxm(im1tn+1:im1tn+ncap) = dx(1:ncap)
      dwm(im1tn+1:im1tn+ncap) = dw(1:ncap)

    end do

    dxm(mtncap+1:mtncap+mp) = dxp(1:mp)
    dwm(mtncap+1:mtncap+mp) = dyp(1:mp)

    mtnpmp = mtncap + mp

    do k = 1, 2 * n

      dsum = 0.0D+00

      do i = 1, mtnpmp

        dp1 = 0.0D+00
        dp = 1.0D+00

        if ( 1 < k ) then

          do l = 1, k - 1
            dpm1 = dp1
            dp1 = dp
            dp = ( dxm(i) - da(l) ) * dp1 - db(l) * dpm1
          end do

        end if

        dsum = dsum + dwm(i) * dp

      end do

      dnu(k) = dsum

    end do

    call cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, ierr )

    done = .true.

    do k = 1, n
      if ( deps * abs ( dbeta(k) ) < abs ( dbeta(k) - dbe(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do

  return
end
subroutine mcdis_r8 ( n, ncapm, mc, mp, dxp, dyp, quad_r8, deps, iq, idelta, &
  irout, finld, finrd, dendl, dendr, dxfer, dwfer, dalpha, dbeta, ncap, &
  kount, ierrd, ied )

!*****************************************************************************80
!
!! MCDIS_R8 is a multiple-component discretization procedure.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) mc
  integer ( kind = 4 ) mp
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncapm

  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) dbe(n)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dendl(mc)
  real ( kind = 8 ) dendr(mc)
  real ( kind = 8 ) deps
  logical done
  real ( kind = 8 ) dp0(mc*ncapm+mp)
  real ( kind = 8 ) dp1(mc*ncapm+mp)
  real ( kind = 8 ) dp2(mc*ncapm+mp)
  real ( kind = 8 ) dw(ncapm)
  real ( kind = 8 ) dwfer(ncapm)
  real ( kind = 8 ) dwm(mc*ncapm+mp)
  real ( kind = 8 ) dx(ncapm)
  real ( kind = 8 ) dxfer(ncapm)
  real ( kind = 8 ) dxm(mc*ncapm+mp)
  real ( kind = 8 ) dxp(mp)
  real ( kind = 8 ) dyp(mp)
  logical finld
  logical finrd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idelta
  integer ( kind = 4 ) ied
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ierrd
  integer ( kind = 4 ) im1tn
  integer ( kind = 4 ) incap
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) irout
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) mtncap
  integer ( kind = 4 ) ncap
  external quad_r8

  if ( idelta <= 0 ) then
    idelta = 1
  end if

  if ( n < 1 ) then
    ierrd = - 1
    return
  end if

  incap = 1
  kount = -1
  ierr = 0
  dbeta(1:n) = 0.0D+00

  ncap = ( 2 * n - 1 ) / idelta

  do

    dbe(1:n) = dbeta(1:n)

    kount = kount + 1

    if ( 1 < kount ) then
      incap = 2**( kount / 5 ) * n
    end if

    ncap = ncap + incap

    if ( ncapm < ncap ) then
      ierrd = ncapm
      return
    end if

    mtncap = mc * ncap

    do i = 1, mc

      im1tn = ( i - 1 ) * ncap

      if ( iq == 1 ) then
        call quad_r8 ( ncap, dx, dw, i, ierr )
      else
        call qgp_r8 ( ncap, dx, dw, i, ierr, mc, finld, finrd, dendl, &
          dendr, dxfer, dwfer )
      end if

      if ( ierr /= 0 ) then
        ierrd = i
        return
      end if

      dxm(im1tn+1:im1tn+ncap) = dx(1:ncap)
      dwm(im1tn+1:im1tn+ncap) = dw(1:ncap)

    end do

    dxm(mtncap+1:mtncap+mp) = dxp(1:mp)
    dwm(mtncap+1:mtncap+mp) = dyp(1:mp)

    if ( irout == 1 ) then
      call sti_r8 ( n, mtncap+mp, dxm, dwm, dalpha, dbeta, ied, dp0, dp1, dp2 )
    else
      call lancz_r8 ( n, mtncap+mp, dxm, dwm, dalpha, dbeta, ied, dp0, dp1 )
    end if

    done = .true.

    do k = 1, n
      if ( deps * abs ( dbeta(k) ) < abs ( dbeta(k) - dbe(k) ) ) then
        done = .false.
        exit
      end if
    end do

    if ( done ) then
      exit
    end if

  end do

  return
end
function nu0her ( n, z, eps )

!*****************************************************************************80
!
!! NU0HER estimates a starting index for recursion with the Hermite measure.
!
!  Discussion:
!
!    This routine provides a starting backward recurrence index for the Hermite
!    measure that can be used in place of NU0 in the routines KNUM_R4 and
!    KNUM_R8.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) eps
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nu0her
  complex ( kind = 8 ) z

  nu0her = int ( 2.0D+00 * ( sqrt ( 0.5D+00 * real ( n + 1, kind = 8 ) ) &
    + 0.25D+00 * log ( 1.0D+00 / eps ) / abs ( aimag ( z ) ) )**2 )

  return
end
function nu0jac ( n, z, eps )

!*****************************************************************************80
!
!! NU0JAC estimates a starting index for recursion with the Jacobi measure.
!
!  Discussion:
!
!    This routine provides a starting backward recurrence index for the Jacobi
!    measure that can be used in place of NU0 in the routines KNUM_R4
!    and KNUM_R8.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) eps
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nu0jac
  real ( kind = 8 ) pi
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) y2
  complex ( kind = 8 ) z

  pi = 4.0D+00 * atan ( 1.0D+00 )
  x = real ( z, kind = 8 )
  y = abs ( aimag ( z ) )

  if ( x < -1.0D+00 ) then
    angle = 0.5D+00 * ( 2.0D+00 * pi + atan ( y / ( x - 1.0D+00 ) ) &
      + atan ( y / ( x + 1.0D+00 ) ) )
  else if ( x == -1.0D+00 ) then
    angle = 0.5D+00 * ( 1.5D+00 * pi - atan ( 0.5D+00 * y ) )
  else if ( x < 1.0D+00 ) then
    angle = 0.5D+00 * ( pi + atan ( y / ( x - 1.0D+00 ) ) &
      + atan ( y / ( x + 1.0D+00 ) ) )
  else if ( x == 1.0D+00 ) then
    angle = 0.5D+00 * ( 0.5D+00 * pi + atan ( 0.5D+00 * y ) )
  else if ( 1.0D+00 < x ) then
    angle = 0.5D+00 * ( atan ( y / ( x - 1.0D+00 ) ) &
      + atan ( y / ( x + 1.0D+00 ) ) )
  end if

  x2 = x * x
  y2 = y * y
  r = ( ( x2 - y2 - 1.0D+00 )**2 + 4.0D+00 * x2 * y2 )**0.25D+00
  r = sqrt ( ( x + r * cos ( angle ) )**2 + ( y + r * sin ( angle ) )**2)

  nu0jac = int ( real ( n + 1, kind = 8 ) + 0.5D+00 * log ( 1.0D+00 / eps ) &
    / log ( r ) )

  return
end
function nu0lag ( n, z, al, eps )

!*****************************************************************************80
!
!! NU0LAG estimates a starting index for recursion with the Laguerre measure.
!
!  Discussion:
!
!    This routine provides a starting backward recurrence index for the
!    Laguerre measure that can be used in place of NU0 in the routines KNUM_R4
!    and KNUM_R8.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) al
  real ( kind = 8 ) eps
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nu0lag
  real ( kind = 8 ) phi
  real ( kind = 8 ) pi
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  complex ( kind = 8 ) z

  pi = 4.0D+00 * atan ( 1.0D+00 )
  x = real ( z, kind = 8 )
  y = aimag ( z )

  if ( y < 0.0D+00 ) then
    phi = 1.5D+00 * pi
  else
    phi = 0.5D+00 * pi
  end if

  if ( x /= 0.0D+00 ) then

    phi = atan ( y / x )

    if ( y <= 0.0D+00 .or. x <= 0.0D+00 ) then

      phi = phi + pi

      if ( 0.0D+00 <= x ) then
        phi = phi + pi
      end if

    end if

  end if

  nu0lag = int ( &
    ( sqrt ( real ( n + 1, kind = 8 ) + 0.5D+00 * ( al + 1.0D+00 ) ) &
    + log ( 1.0D+00 / eps ) &
    / ( 4.0D+00 * ( x * x + y * y )**0.25D+00 &
    * cos ( 0.5D+00 * ( phi - pi ) ) ) &
    )**2 - 0.5D+00 * ( al + 1.0D+00 ) )

  return
end
subroutine qgp_r8 ( n, dx, dw, i, ierr, mcd, finld, finrd, dendl, dendr, &
  dxfer, dwfer )

!*****************************************************************************80
!
!! QGP_R8 is a general-purpose discretization routine.
!
!  Discussion:
!
!    The user has to supply the routine
!
!       function wf_r8 ( dx, i )
!
!    which evaluates the weight function in real ( kind = 8 ) at the
!    point  dx  on the i-th component interval.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) mcd
  integer ( kind = 4 ) n

  real ( kind = 8 ) dendl(mcd)
  real ( kind = 8 ) dendr(mcd)
  real ( kind = 8 ) dphi
  real ( kind = 8 ) dphi1
  real ( kind = 8 ) dw(n)
  real ( kind = 8 ) dwfer(*)
  real ( kind = 8 ) dx(n)
  real ( kind = 8 ) dxfer(*)
  logical finld
  logical finrd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  real ( kind = 8 ) wf_r8
!
!  The arrays  dxfer,dwfer  are dimensioned in the routine  dmcdis.
!
  ierr = 0

  if ( i == 1 ) then
    call fejer_r8 ( n, dxfer, dwfer )
  end if

  if ( 1 < i .and. i < mcd ) then
    go to 60
  end if

  if ( mcd == 1 ) then

    if ( finld .and. finrd ) go to 60
    if ( finld ) go to 20
    if ( finrd ) go to 40

    do k = 1, n
      call symtr_r8 ( dxfer(k), dphi, dphi1 )
      dx(k) = dphi
      dw(k) = dwfer(k) * wf_r8 ( dphi, i ) * dphi1
    end do

    return

  else

    if ( (i == 1 .and. finld ) .or. ( i == mcd .and. finrd ) ) go to 60
    if ( i == 1 ) go to 40

  end if

   20 continue

  do k = 1, n
    call tr_r8 ( dxfer(k), dphi, dphi1 )
    dx(k) = dendl(mcd) + dphi
    dw(k) = dwfer(k) * wf_r8 ( dx(k), mcd ) * dphi1
  end do

  return

   40 continue

  do k = 1, n
    call tr_r8 ( -dxfer(k), dphi, dphi1 )
    dx(k) = dendr(1) - dphi
    dw(k) = dwfer(k) * wf_r8 ( dx(k), 1 ) * dphi1
  end do

  return

   60 continue

  do k = 1, n
    dx(k) = 0.5D+00 * ( ( dendr(i) - dendl(i) ) &
      * dxfer(k) + dendr(i) + dendl(i) )
    dw(k) = 0.5D+00 * ( dendr(i) - dendl(i) ) * dwfer ( k ) * wf_r8 ( dx(k), i )
  end do

  return
end
subroutine radau_r8 ( n, alpha, beta, endl, zero, weigh, ierr )

!*****************************************************************************80
!
!! RADAU_R8 generates a Gauss-Radau quadrature formula.
!
!  Discussion:
!
!    Given N and a measure  dlambda, this routine generates the
!    (N+1)-point Gauss-Radau quadrature formula
!
!      integral over supp(dlambda) of f(t)dlambda(t)
!
!      =  w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!
!    The nodes are returned as  zero(k) = x(k), the weights as  weight(k)
!    = w(k), k=0,1,2,...,n. The user has to supply the recursion
!    coefficients  alpha(k), beta(k), k = 0,1,2,...,n, for the measure
!    dlambda. The nodes and weights are computed as eigenvalues and
!    in terms of the first component of the respective normalized
!    eigenvectors of a slightly modified Jacobi matrix of order  n+1.
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of interior points in the 
!    Gauss-Radau formula.
!
!    Input, real ( kind = 8 ) ALPHA(N+1), BETA(N+1), the recursion coefficients
!    alpha(k-1), beta(k-1), k = 1,2,...,n+1; the coefficient  alpha(n+1)  is not
!    used by the routine.
!
!    Input, real ( kind = 8 ) ENDL, the prescribed endpoint x(0) of the 
!    Gauss-Radau formula.
!
!    Output, real ( kind = 8 ) ZERO(N+1), the nodes (in increasing order)
!    zero(k) = x(k), k=0,1,2,...,n
!
!    Output, real ( kind = 8 ) WEIGHT(N+1), the weights weight(k) = w(k),
!    k=0,1,2,...,n
!
!    Output, integer ( kind = 4 ) IERR, an error flag from GAUSS_R8.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n+1)
  real ( kind = 8 ) alpha(n+1)
  real ( kind = 8 ) b(n+1)
  real ( kind = 8 ) beta(n+1)
  real ( kind = 8 ) e(n+1)
  real ( kind = 8 ) endl
  real ( kind = 8 ) epsma
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  real ( kind = 8 ) p0
  real ( kind = 8 ) p1
  real ( kind = 8 ) pm1
  real ( kind = 8 ) weigh(n+1)
  real ( kind = 8 ) zero(n+1)

  epsma = epsilon ( epsma )

  a(1:n) = alpha(1:n)
  b(1:n+1) = beta(1:n+1)
!
!  Determine A(N+1).
!
  p0 = 0.0D+00
  p1 = 1.0D+00
  do k = 1, n
    pm1 = p0
    p0 = p1
    p1 = ( endl - a(k) ) * p0 - b(k) * pm1
  end do

  a(n+1) = endl - b(n+1) * p0 / p1
!
!  Call GAUSS_R8.
!
  call gauss_r8 ( n+1, a, b, epsma, zero, weigh, ierr, e )

  return
end
subroutine recur_r8 ( n, ipoly, dal, dbe, da, db, ierr )

!*****************************************************************************80
!
!! RECUR_R8 generates recursion coefficients for orthogonal polynomials.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of recursion coefficients
!    desired.
!
!    Input, integer ( kind = 4 ) IPOLY, identifies the polynomial.
!    1, Legendre polynomial on [-1,+1];
!    2, Legendre polynomial on [0,+1];
!    3, Chebyshev polynomial of the first kind;
!    4, Chebyshev polynomial of the second kind;
!    5, Chebyshev polynomial of the third kind.
!    6, Jacobia polynomial with parameters AL, BE.
!    7, Generalized Laguerre polynomial with parameter AL.
!    8, Hermite polynomial.
!    9, Chebyshev polynomial of the fourth kind.
!
!    Input, real ( kind = 8 ) DAL, a parameter needed if IPOLY is 6 or 7.
!
!    Input, real ( kind = 8 ) DBE, a parameter needed if IPOLY is 6.
!
!    Output, real ( kind = 8 ) DA(N), DB(N), the recursion coefficients.
!
!    Output, integer ( kind = 4 ) IERR, the error flag.
!    0, normal return.
!    1, AL or BE is out of range.
!    2, there is a potential overflow evaluating BETA0 for IPOLY = 6 or 7.
!    3, N is out of range.
!    4, IPOLY is not an admissible value.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alga_r8
  real ( kind = 8 ) da(n)
  real ( kind = 8 ) dal
  real ( kind = 8 ) dal2
  real ( kind = 8 ) dalpbe
  real ( kind = 8 ) db(n)
  real ( kind = 8 ) dbe
  real ( kind = 8 ) dbe2
  real ( kind = 8 ) dkm1
  real ( kind = 8 ) dlmach
  real ( kind = 8 ) dt
  real ( kind = 8 ), external :: gamma_r8
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) k

  if ( n < 1 ) then
    ierr = 3
    return
  end if

  dlmach = log ( huge ( dlmach ) )
  ierr = 0
  da(1:n) = 0.0D+00
!
!  Legendre on [-1,+1].
!
  if ( ipoly == 1 ) then

    db(1) = 2.0D+00
    do k = 2, n
      dkm1 = real ( k - 1, kind = 8 )
      db(k) = 1.0D+00 / ( 4.0D+00 - 1.0D+00 / ( dkm1 * dkm1 ) )
    end do
!
!  Legendre on [0,+1].
!
  else if ( ipoly == 2 ) then

    da(1) = 0.5D+00
    db(1) = 1.0D+00
    do k = 2, n
      da(k) = 0.5D+00
      dkm1 = real ( k - 1, kind = 8 )
      db(k) = 0.25D+00 / ( 4.0D+00 - 1.0D+00 / ( dkm1 * dkm1 ) )
    end do
!
!  Chebyshev polynomial of the first kind.
!
  else if ( ipoly == 3 ) then

    db(1) = 4.0D+00 * atan ( 1.0D+00 )
    if ( n == 1 ) then
      return
    end if
    db(2) = 0.5D+00
    db(3:n) = 0.25D+00
!
!  Chebyshev polynomial of the second kind.
!
  else if ( ipoly == 4 ) then

    db(1) = 2.0D+00 * atan ( 1.0D+00 )
    db(2:n) = 0.25D+00
!
!  Chebyshev polynomial of the third kind.
!
  else if ( ipoly == 5 ) then

    db(1) = 4.0D+00 * atan ( 1.0D+00 )
    da(1) = 0.5D+00
    db(2:n) = 0.25D+00
!
!  Jacobi polynomial with parameters AL and BE.
!
  else if ( ipoly == 6 ) then

    if ( dal <= -1.0D+00 .or. dbe <= -1.0D+00 ) then
      ierr = 1
      return
    end if

    dalpbe = dal + dbe
    da(1) = ( dbe - dal ) / ( dalpbe + 2.0D+00 )
    dt = ( dalpbe + 1.0D+00 ) * log ( 2.0D+00 ) + alga_r8 ( dal + 1.0D+00 ) &
      + alga_r8 ( dbe + 1.0D+00 ) - alga_r8 ( dalpbe + 2.0D+00 )

    if ( dlmach < dt ) then
      ierr = 2
      db(1) = huge ( db(1) )
    else
      db(1) = exp ( dt )
    end if

    if ( n == 1 ) then
      return
    end if

    dal2 = dal * dal
    dbe2 = dbe * dbe
    da(2) = ( dbe2 - dal2 ) / ( ( dalpbe + 2.0D+00 ) * ( dalpbe + 4.0D+00 ) )
    db(2) = 4.0D+00 * ( dal + 1.0D+00 ) * ( dbe + 1.0D+00 ) &
      / ( ( dalpbe + 3.0D+00 ) * ( dalpbe + 2.0D+00 )**2 )

    do k = 3, n

      dkm1 = real ( k - 1, kind = 8 )

      da(k) = 0.25D+00 * ( dbe2 - dal2 ) &
        / ( dkm1 * dkm1 * ( 1.0D+00 + 0.5D+00 * dalpbe / dkm1 ) &
        * ( 1.0D+00 + 0.5D+00 * ( dalpbe + 2.0D+00 ) / dkm1 ) )

      db(k) = 0.25D+00 * ( 1.0D+00 + dal / dkm1 ) &
        * ( 1.0D+00 + dbe / dkm1 ) * ( 1.0D+00 + dalpbe / dkm1 ) &
        / ( ( 1.0D+00 + 0.5D+00 * ( dalpbe + 1.0D+00 ) / dkm1 ) &
        * ( 1.0D+00 + 0.5D+00 * ( dalpbe -1.0D+00 ) / dkm1 ) &
        * ( 1.0D+00 + 0.5D+00 * dalpbe / dkm1 )**2 )

    end do
!
!  Generalized Laguerre polynomial with parameter AL.
!
  else if ( ipoly == 7 ) then

    if ( dal <= -1.0D+00 ) then
      ierr = 1
      return
    end if

    da(1) = dal + 1.0D+00
    db(1) = gamma_r8 ( dal + 1.0D+00, ierr )

    if ( ierr == 2 ) then
      db(1) = huge ( db(1) )
      return
    end if

    do k = 2, n
      dkm1 = real ( k - 1, kind = 8 )
      da(k) = 2.0D+00 * dkm1 + dal + 1.0D+00
      db(k) = dkm1 * ( dkm1 + dal )
    end do
!
!  Hermite polynomial.
!
  else if ( ipoly == 8 ) then

    db(1) = sqrt ( 4.0D+00 * atan ( 1.0D+00 ) )

    do k = 2, n
      db(k) = 0.5D+00 * real ( k - 1, kind = 8 )
    end do
!
!  Chebyshev polynomial of the fourth kind.
!
  else if ( ipoly == 9 ) then

    db(1) = 4.0D+00 * atan ( 1.0D+00 )
    da(1) = - 0.5D+00
    db(2:n) = 0.25D+00
!
!  Inadmissible value of IPOLY.
!
  else

    ierr = 4

  end if

  return
end
subroutine sti_r8 ( n, ncap, dx, dw, dalpha, dbeta, ierr, dp0, dp1, dp2 )

!*****************************************************************************80
!
!! STI_R8 applies Stieltjes's procedure.
!
!  Modified:
!
!    05 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    On Generating Orthogonal Polynomials,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 3, Number 3, 1982, pages 289-317.
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncap

  real ( kind = 8 ) dalpha(n)
  real ( kind = 8 ) dbeta(n)
  real ( kind = 8 ) dp0(ncap)
  real ( kind = 8 ) dp1(ncap)
  real ( kind = 8 ) dp2(ncap)
  real ( kind = 8 ) dsum0
  real ( kind = 8 ) dsum1
  real ( kind = 8 ) dsum2
  real ( kind = 8 ) dt
  real ( kind = 8 ) dw(ncap)
  real ( kind = 8 ) dx(ncap)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m

  ierr = 0

  if ( n <= 0 .or. ncap < n ) then
    ierr = 1
    return
  end if

  dsum0 = sum ( dw(1:ncap) )
  dsum1 = dot_product ( dw(1:ncap), dx(1:ncap) )

  dalpha(1) = dsum1 / dsum0
  dbeta(1) = dsum0

  if ( n == 1 ) then
    return
  end if

  dp1(1:ncap) = 0.0D+00
  dp2(1:ncap) = 1.0D+00

  do k = 1, n - 1

    dsum1 = 0.0D+00
    dsum2 = 0.0D+00

    do m = 1, ncap

      if ( dw(m) /= 0.0D+00 ) then

        dp0(m) = dp1(m)
        dp1(m) = dp2(m)
        dp2(m) = ( dx(m) - dalpha(k) ) * dp1(m) - dbeta(k) * dp0(m)

        if ( 0.1D+00 * huge ( dp2(m) ) < abs ( dp2(m) ) .or. &
          0.1D+00 * huge ( dsum2 ) < abs ( dsum2 ) ) then
          ierr = k
          return
        end if

        dt = dw(m) * dp2(m) * dp2(m)
        dsum1 = dsum1 + dt
        dsum2 = dsum2 + dt * dx(m)

      end if

    end do

    if ( abs ( dsum1 ) < 10.0D+00 * tiny ( dsum1 ) ) then
      ierr = - k
      return
    end if

    dalpha(k+1) = dsum2 / dsum1
    dbeta(k+1) = dsum1 / dsum0
    dsum0 = dsum1

  end do

  return
end
subroutine symtr_r8 ( t, phi, phi1 )

!*****************************************************************************80
!
!! SYMTR_R8 maps T in [-1,1] to X in (-oo,oo).
!
!  Discussion:
!
!    X = T / ( 1 - T * T )
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the point in [-1,1] to be mapped.
!    -1 < T < 1.
!
!    Output, real ( kind = 8 ) PHI, the value of X(T).
!
!    Output, real ( kind = 8 ) PHI1, the derivative of the mapping.
!
  implicit none

  real ( kind = 8 ) phi
  real ( kind = 8 ) phi1
  real ( kind = 8 ) t

  phi = t / ( 1.0D+00 - t * t )
  phi1 = ( t * t + 1.0D+00 ) / ( t * t - 1.0D+00 )**2

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
subroutine tr_r8 ( t, phi, phi1 )

!*****************************************************************************80
!
!! TR_R8 maps T in [-1,1] to X in [0,oo).
!
!  Modified:
!
!    10 August 2007
!
!  Author:
!
!    Walter Gautschi
!
!  Reference:
!
!    Walter Gautschi,
!    Algorithm 726:
!    ORTHPOL - A Package of Routines for Generating Orthogonal
!    Polynomials and Gauss-Type Quadrature Rules,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 1, March 1994, pages 21-62.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T.
!
!    Output, real ( kind = 8 ) PHI, the value of X(T).
!
!    Output, real ( kind = 8 ) PHI1, the derivative of the mapping.
!
  implicit none

  real ( kind = 8 ) phi
  real ( kind = 8 ) phi1
  real ( kind = 8 ) t

  phi = ( 1.0D+00 + t ) / ( 1.0D+00 - t )
  phi1 = 2.0D+00 / ( t - 1.0D+00 )**2

  return
end
