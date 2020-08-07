subroutine diffusivity_1d_pwc ( nc, xc, vc, np, xp, vp )

!*****************************************************************************80
!
!! DIFFUSIVITY_1D_PWC: piecewise constant diffusivity function in 1D.
!
!  Discussion:
!
!    A piecewise constant function is defined over NC intervals, 
!    with interval IC associated with constant value VC(I),
!    separated by NC-1 ascending sorted breakpoints.
!
!    The function is to be evaluated at NP points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2019
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NC, the number of function values.
!
!    Input, real ( kind = 8 ) XC(NC-1), the breakpoints, in ascending order.
!
!    Input, real ( kind = 8 ) VC(NC), the function values over each interval.
!
!    Input, integer ( kind = 4 ) NP, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XP(NP), the evaluation points.
!
!    Output, real ( kind = 8 ) VP(NP), the function value at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) nc
  integer ( kind = 4 ) np

  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) kc
  real ( kind = 8 ) vc(nc)
  real ( kind = 8 ) vp(np)
  real ( kind = 8 ) xc(nc-1)
  real ( kind = 8 ) xp(np)

  do ip = 1, np

    kc = 1
    do ic = 1, nc - 1
      if ( xp(ip) < xc(ic) ) then
        exit
      end if
      kc = kc + 1
    end do

    vp(ip) = vc(kc)

  end do

  return
end
subroutine diffusivity_1d_xk ( dc0, m, omega, n, x, dc )

!*****************************************************************************80
!
!! DIFFUSIVITY_1D_XK evaluates a 1D stochastic diffusivity function.
!
!  Discussion:
!
!    The 1D diffusion equation has the form
!
!      - d/dx ( DC(X) Del U(X) ) = F(X)
!
!    where DC(X) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - d/dx ( DC(X;OMEGA) d/dx U(X) ) = F(X).
!
!    In this function, the domain is assumed to be the unit interval [0.1].
!
!
!    For DC0 = 1 and F(X) = 0, with boundary conditions U(0:OMEGA) = 0,
!    U(1;OMEGA) = 1, the exact solution is
!
!    If OMEGA ~= 0:
!
!      U(X;OMEGA) = log ( 1 + OMEGA * X ) / log ( 1 + OMEGA )
!
!    If OMEGA = 0:
!
!      U(X;OMEGA) = X
!
!    In the numerical experiments described in the paper, OMEGA was taken
!    to be a random variable with a Beta, or Uniform, or Gaussian or 
!    Poisson or Binomial distribution.
!
!    For the Gaussian and Poisson distributions, the positivity requirement 
!    could not be guaranteed, and the experiments were simply made with a 
!    "small" variance of 0.1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu, George Karniadakis,
!    Modeling uncertainty in steady state diffusion problems via
!    generalized polynomial chaos,
!    Computer Methods in Applied Mechanics and Engineering,
!    Volume 191, 2002, pages 4927-4948.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.
!
!    Input, integer ( kind = 4 ) M, the number of stochastic parameters.
!
!    Input, real ( kind = 8 ) OMEGA(M), the stochastic parameters.  
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), the point where the diffusion coefficient 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) DC(N), the value of the diffusion coefficient 
!    at X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) dc(n)
  real ( kind = 8 ) dc0
  integer ( kind = 4 ) k
  real ( kind = 8 ) omega(m)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) w
  real ( kind = 8 ) x(n)

  k = 0
  w = 1.0D+00

  dc(1:n) = 0.0D+00

  do while ( k < m )

    if ( k < m ) then
      k = k + 1
      dc(1:n) = dc(1:n) + omega(k) * sin ( w * r8_pi * x(1:n) )
    end if

    if ( k < m ) then
      k = k + 1
      dc(1:n) = dc(1:n) + omega(k) * cos ( w * r8_pi * x(1:n) )    
    end if

    w = w + 1.0D+00

  end do

  dc(1:n) = exp ( - 0.125D+00 ) * dc(1:n)

  dc(1:n) = dc0 + exp ( dc(1:n) )

  return
end
subroutine diffusivity_2d_bnt ( dc0, omega, n, x, y, dc )

!*****************************************************************************80
!
!! DIFFUSIVITY_2D_BNT evaluates a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The 2D diffusion equation has the form
!
!      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
!
!    where DC(X,Y) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
!
!    In this function, the domain is the rectangle [-1.5,0]x[-0.4,0.8].
!
!    The four stochastic parameters OMEGA(1:4) are assumed to be independent
!    identically distributed random variables with mean value zero and 
!    variance 1.  The distribution is typically taken to be Gaussian or
!    uniform.
!
!    A collocation approach to this problem would then use the roots of
!    Hermite or Legendre polynomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ivo Babuska, Fabio Nobile, Raul Tempone,
!    A stochastic collocation method for elliptic partial differential equations
!    with random input data,
!    SIAM Journal on Numerical Analysis,
!    Volume 45, Number 3, 2007, pages 1005-1034.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.  Take DC0 = 10.
!
!    Input, real ( kind = 8 ) OMEGA(4), the stochastic parameters.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the points where the diffusion 
!    coefficient is to be evaluated.
!
!    Output, real ( kind = 8 ) DC(N), the value of the diffusion coefficient 
!    at (X,Y).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg(n)
  real ( kind = 8 ) dc(n)
  real ( kind = 8 ) dc0
  real ( kind = 8 ) omega(4)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  arg(1:n) = omega(1) * cos ( r8_pi * x(1:n) ) &
           + omega(2) * sin ( r8_pi * x(1:n) ) &
           + omega(3) * cos ( r8_pi * y(1:n) ) &
           + omega(4) * sin ( r8_pi * y(1:n) )

  arg(1:n)= exp ( - 0.125D+00 ) * arg(1:n)

  dc(1:n) = dc0 + exp ( arg(1:n) )

  return
end
subroutine diffusivity_2d_elman ( a, cl, dc0, m_1d, omega, n1, n2, x, y, dc )

!*****************************************************************************80
!
!! DIFFUSIVITY_2D_ELMAN evaluates a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The 2D diffusion equation has the form
!
!      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
!
!    where DC(X,Y) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
!
!    In this function, the domain is assumed to be the square [-A,+A]x[-A,+A].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Elman, Darran Furnaval,
!    Solving the stochastic steady-state diffusion problem using multigrid,
!    IMA Journal on Numerical Analysis,
!    Volume 27, Number 4, 2007, pages 675-688.
!
!    Roger Ghanem, Pol Spanos,
!    Stochastic Finite Elements: A Spectral Approach,
!    Revised Edition,
!    Dover, 2003,
!    ISBN: 0486428184,
!    LC: TA347.F5.G56.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the "radius" of the square region.  The region
!    is assumed to be [-A,+A]x[-A,+A].
!    0 < A.
!
!    Input, real ( kind = 8 ) CL, the correlation length.
!    0 < CL.
!
!    Input, real ( kind = 8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.  Take DC0 = 10.
!
!    Input, integer ( kind = 4 ) M_1D, the first and second dimensions of the
!    stochastic parameter array.
!
!    Input, real ( kind = 8 ) OMEGA(M_1D,M_1D), the stochastic parameters.
!
!    Input, integer ( kind = 4 ) N1, N2, the dimensions of the X and Y arrays.
!
!    Input, real ( kind = 8 ) X(N1,N2), Y(N1,N2), the points where the
!    diffusion coefficient is to be evaluated.
!
!    Output, real ( kind = 8 ) DC(N1,N2), the value of the diffusion 
!    coefficient at X.
!
  implicit none

  integer ( kind = 4 ) m_1d
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a
  real ( kind = 8 ) c_1dx(m_1d,n1,n2)
  real ( kind = 8 ) c_1dy(m_1d,n1,n2)
  real ( kind = 8 ) cl
  real ( kind = 8 ) dc(n1,n2)
  real ( kind = 8 ) dc0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) lambda_1d(m_1d)
  integer ( kind = 4 ) m
  real ( kind = 8 ) omega(m_1d,m_1d)
  real ( kind = 8 ) theta_1d(m_1d)
  real ( kind = 8 ) x(n1,n2)
  real ( kind = 8 ) y(n1,n2)

  m = m_1d * m_1d
!
!  Compute THETA.
!
  call theta_solve ( a, cl, m_1d, theta_1d )
!
!  Compute LAMBDA_1D.
!
  lambda_1d(1:m_1d) = 2.0D+00 * cl &
    / ( 1.0D+00 + cl * cl * theta_1d(1:m_1d) ** 2 )
!
!  Compute C_1DX(1:M1D) and C_1DY(1:M1D) at (X,Y).
!
  c_1dx(1:m_1d,1:n1,1:n2) = 0.0D+00
  c_1dy(1:m_1d,1:n1,1:n2) = 0.0D+00

  k = 0

  do

    if ( m_1d <= k ) then
      exit
    end if

    k = k + 1

    c_1dx(k,1:n1,1:n2) = cos ( theta_1d(k) * a * x(1:n1,1:n2) ) &
      / sqrt ( a + sin ( 2.0D+00 * theta_1d(k) * a ) &
      / ( 2.0D+00 * theta_1d(k) ) )

    c_1dy(k,1:n1,1:n2) = cos ( theta_1d(k) * a * y(1:n1,1:n2) ) &
      / sqrt ( a + sin ( 2.0D+00 * theta_1d(k) * a ) &
      / ( 2.0D+00 * theta_1d(k) ) )

    if ( m_1d <= k ) then
      exit
    end if

    k = k + 1

    c_1dx(k,1:n1,1:n2) = sin ( theta_1d(k) * a * x(1:n1,1:n2) ) &
      / sqrt ( a - sin ( 2.0D+00 * theta_1d(k) * a ) &
      / ( 2.0D+00 * theta_1d(k) ) )

    c_1dy(k,1:n1,1:n2) = sin ( theta_1d(k) * a * y(1:n1,1:n2) ) &
      / sqrt ( a - sin ( 2.0D+00 * theta_1d(k) * a ) &
      / ( 2.0D+00 * theta_1d(k) ) )

  end do
!
!  Evaluate the diffusion coefficient DC at (X,Y).
!
  dc(1:n1,1:n2) = dc0
  do j = 1, m_1d
    do i = 1, m_1d
      dc(1:n1,1:n2) = dc(1:n1,1:n2) + sqrt ( lambda_1d(i) * lambda_1d(j) ) &
        * c_1dx(i,1:n1,1:n2) * c_1dy(j,1:n1,1:n2) * omega(i,j)
    end do
  end do

  return
end
subroutine diffusivity_2d_ntw ( cl, dc0, m, omega, n, x, y, dc )

!*****************************************************************************80
!
!! DIFFUSIVITY_2D_NTW evaluates a 2D stochastic diffusivity function.
!
!  Discussion:
!
!    The 2D diffusion equation has the form
!
!      - Del ( DC(X,Y) Del U(X,Y) ) = F(X,Y)
!
!    where DC(X,Y) is a function called the diffusivity.
!
!    In the stochastic version of the problem, the diffusivity function
!    includes the influence of stochastic parameters:
!
!      - Del ( DC(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
!
!    In this function, the domain is the rectangle [0,D]x[0,D] where D = 1.
!
!    Note that in this problem the diffusivity has a one-dimensional
!    spatial dependence on X, but not on Y!
!
!    The random variables OMEGA are independent, have zero mean and unit
!    variance, and are uniformly distributed in [-sqrt(3),+sqrt(3)].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Xiang Ma, Nicholas Zabaras,
!    An adaptive hierarchical sparse grid collocation algorithm for the solution
!    of stochastic differential equations,
!    Journal of Computational Physics,
!    Volume 228, pages 3084-3113, 2009.
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) CL, the desired physical correlation length for 
!    the coefficient.
!
!    Input, real ( kind = 8 ) DC0, the constant term in the expansion of the 
!    diffusion coefficient.  Take DC0 = 0.5.
!
!    Input, integer ( kind = 4 ) M, the number of terms in the expansion.
!
!    Input, real ( kind = 8 ) OMEGA(M), the stochastic parameters.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the points where the diffusion 
!    coefficient is to be evaluated.
!
!    Output, real ( kind = 8 ) DC(N), the value of the diffusion coefficient
!    at (X,Y).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) cl
  real ( kind = 8 ) d
  real ( kind = 8 ) dc(n)
  real ( kind = 8 ) dc_arg(n)
  real ( kind = 8 ) dc0
  integer ( kind = 4 ) i
  real ( kind = 8 ) ihalf_r8
  real ( kind = 8 ) l
  real ( kind = 8 ) lp
  real ( kind = 8 ) omega(m)
  real ( kind = 8 ) phi(n)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) zeta
  real ( kind = 8 ) zeta_arg

  d = 1.0D+00
  lp = max ( d, 2.0D+00 * cl )
  l = cl / lp

  dc_arg(1:n) = 1.0D+00 + omega(1) * sqrt ( sqrt ( r8_pi ) * l / 2.0D+00 )

  do i = 2, m

    ihalf_r8 = real ( i / 2, kind = 8 )
    zeta_arg = - ( ihalf_r8 * r8_pi * l ) ** 2 / 8.0D+00
    zeta = sqrt ( sqrt ( r8_pi ) * l ) * exp ( zeta_arg )

    if ( mod ( i, 2 ) == 0 ) then
      phi(1:n) = sin ( ihalf_r8 * r8_pi * x(1:n) / lp )
    else
      phi(1:n) = cos ( ihalf_r8 * r8_pi * x(1:n) / lp )
    end if

    dc_arg(1:n) = dc_arg(1:n) + zeta * phi(1:n) * omega(i)

  end do

  dc(1:n) = dc0 + exp ( dc_arg(1:n) )

  return
end
subroutine diffusivity_2d_pwc ( h, w, a, b, c, d, omega, n, x, y, rho )

!*****************************************************************************80
!
!! DIFFUSIVITY_2D_PWC: piecewise constant diffusivity function.
!
!  Discussion:
!
!    The 2D stochastic diffusion equation has the form
!
!      - Del ( RHO(X,Y;OMEGA) Del U(X,Y;OMEGA) ) = F(X,Y).
!
!    Here, the diffusivity RHO is assumed to be a piecewise constant function,
!    defined on rectangular grid over [A,B]x[C,D] that is H cells high
!    and W cells wide.  The parameters OMEGA(H*W) are assumed to be given 
!    positive values that are bounded away from 0.
!
!    The underlying grid is assumed to be equally spaced, so each cell
!    is (D-C)/M units high and (B-A)/N wide.
!
!       ^  ^  +-------+-------+-------+-------+
!       D  H  | (H,1) | (H,2) | ----- | (H,W) |
!       |  |  +-------+-------+-------+-------+
!             | ----- | ----- | (I,J) | ----- |
!       Y  I  +-------+-------+-------+-------+
!             | (2,1) | (2,2) | ----- | (2,W) |
!       |  |  +-------+-------+-------+-------+
!       |  |  | (1,1) | (1,2) | ----- | (1,W) |
!       C  1  +-------+-------+-------+-------+
!       |  |
!       +--+--1-------------- J --------------W->
!       +--+--A-------------- X --------------B->
!
!    Indexing is tricky, since our (X,Y) coordinate system indexing and
!    conventions for indexing an array and for ordering indices do not
!    match up.  An arbitrary choice must be made, and here we associate
!    "width" and "X" and "J", as well as "height" and "Y" and "I".
!    We start in the lower left corner, and proceed right first, and
!    then move up one row.
!
!    The following coordinate systems are available:
!
!    * cell index (I,J):  (1,1) <= (I,J) <= (H,W). 
!      I refers to "height" and J to "width".
!
!    * cell counter K:    1     <= K     <= H*W
!      Count from lower left to lower right, then
!      go up one row.
!
!    * physical cell center coordinates (XC,YC)
!      the (X,Y) coordinates of the centers of a cell
!
!    * physical coordinates (X,Y) in [A,B]x[C,D].
!      the (X,Y) coordinates of any point
!      
!    * normalized coordinates (X01,Y01) in [0,1]x[0,1].
!   
!    Transformations from one coordinate system to another:
!
!    * (I,J) --> (XC,YC)
!                X01 = (2*J-1)/2/W
!                Y01 = (2*I-1)/2/H
!                XC = A + (B-A) * X01
!                YC = C + (D-C) * Y01
!
!    * (I,J) --> K
!                K = (I-1)*W+J
!
!    * K     --> (I,J)
!                I = floor ( K - 1 ) / W ) + 1
!                J = mod ( K - 1, W ) + 1
!
!    * K     --> (XC,YC)
!                I = floor ( K - 1 ) / W ) + 1
!                J = mod ( K - 1, W ) + 1
!                X01 = (2*J-1)/2/W
!                Y01 = (2*I-1)/2/H
!                XC = A + (B-A) * X01
!                YC = C + (D-C) * Y01
!
!    * (X,Y) --> (I,J)
!                X01 = ( X - A ) / ( B - A )
!                Y01 = ( Y - C ) / ( D - C )
!                I = round ( ( 2 * h * Y01 + 1 ) / 2 );
!                J = round ( ( 2 * w * X01 + 1 ) / 2 );
!
!    * (X,Y) --> K
!                X01 = ( X - A ) / ( B - A )
!                Y01 = ( Y - C ) / ( D - C )
!                I = round ( ( 2 * h * Y01 + 1 ) / 2 );
!                J = round ( ( 2 * w * X01 + 1 ) / 2 );
!                K = (I-1)*W+J
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2019
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Ivo Babuska, Fabio Nobile, Raul Tempone,
!    A stochastic collocation method for elliptic partial differential equations
!    with random input data,
!    SIAM Journal on Numerical Analysis,
!    Volume 45, Number 3, 2007, pages 1005-1034.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) H, W, specifies the number of elements.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of X.
!
!    Input, real ( kind = 8 ) C, D, the lower and upper limits of Y.
!
!    Input, real ( kind = 8 ) OMEGA(H*W), the parameters.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the points where the diffusion coefficient is to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) RHO(N), the value of the diffusion coefficient at each point.
!
  implicit none

  integer ( kind = 4 ) h
  integer ( kind = 4 ) n
  integer ( kind = 4 ) w

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) omega(h*w)
  real ( kind = 8 ) rho(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x01
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y01

  do ij = 1, n

    x01 = ( x(ij) - a ) / ( b - a )
    y01 = ( y(ij) - c ) / ( d - c )

    i = nint ( ( 2.0D+00 * h * y01 + 1.0D+00 ) / 2.0D+00 )
    i = max ( i, 1 )
    i = min ( i, h )

    j = nint ( ( 2.0D+00 * w * x01 + 1.0D+00 ) / 2.0D+00 )
    j = max ( j, 1 )
    j = min ( j, w )

    k = ( i - 1 ) * w + j

    rho(ij) = omega(k)

  end do

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
function r8_epsilon ( )

!*****************************************************************************80
!
!! R8_EPSILON returns the R8 roundoff unit.
!
!  Discussion:
!
!    The roundoff unit is a number R which is a power of 2 with the
!    property that, to the precision of the computer's arithmetic,
!      1 < 1 + R
!    but
!      1 = ( 1 + R / 2 )
!
!    FORTRAN90 provides the superior library routine
!
!      EPSILON ( X )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_EPSILON, the round-off unit.
!
  implicit none

  real ( kind = 8 ) r8_epsilon

  r8_epsilon = 2.220446049250313D-016

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
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
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8mat_max ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_MAX returns the maximum entry of an R8MAT.
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
!    11 December 2004
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
!    Input, real ( kind = 8 ) A(M,N), the M by N matrix.
!
!    Output, real ( kind = 8 ) R8MAT_MAX, the maximum entry of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_max

  r8mat_max = maxval ( a(1:m,1:n) )

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
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
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
subroutine r8vec_max ( n, a, amax )

!*****************************************************************************80
!
!! R8VEC_MAX returns the maximum value in an R8VEC.
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
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), the array.
!
!    Output, real ( kind = 8 ) AMAX, the value of the largest entry.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) amax

  amax = maxval ( a(1:n) )

  return
end
subroutine r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat )

!*****************************************************************************80
!
!! R8VEC_MESH_2D creates a 2D mesh from X and Y vectors.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    NX = 2
!    XVEC = ( 1, 2, 3 )
!    NY = 3
!    YVEC = ( 4, 5 )
!
!    XMAT = (
!      1, 2, 3
!      1, 2, 3 )
!
!    YMAT = (
!      4, 4, 4
!      5, 5, 5 ) 
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2013
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, the number of X and Y values.
!
!    Input, real ( kind = 8 ) XVEC(NX), YVEC(NY), the X and Y coordinate
!    values.
!
!    Output, real ( kind = 8 ) XMAT(NX,NY), YMAT(NX,NY), the coordinate
!    values of points on an NX by NY mesh.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  integer ( kind = 4 ) j
  real ( kind = 8 ) xmat(nx,ny)
  real ( kind = 8 ) xvec(nx)
  real ( kind = 8 ) ymat(nx,ny)
  real ( kind = 8 ) yvec(ny)

  do j = 1, ny
    xmat(1:nx,j) = xvec(1:nx)
  end do

  do j = 1, ny
    ymat(1:nx,j) = yvec(j)
  end do

 return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    This routine can generate a vector of values on one call.  It
!    has the feature that it should provide the same results
!    in the same order no matter how we break up the task.
!
!    Before calling this routine, the user may call RANDOM_SEED
!    in order to set the seed of the random number generator.
!
!    The Box-Muller method is used, which is efficient, but
!    generates an even number of values each time.  On any call
!    to this routine, an even number of new values are generated.
!    Depending on the situation, one value may be left over.
!    
!    However, using saved memory means that the code is not easy to
!    manipulate by changing the seed.  In order to restore that ability,
!    I am suppressing the use of the second computed value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.  If N is
!    negative,then the code will flush its internal memory; in particular,
!    if there is a saved value to be used on the next call, it is
!    instead discarded.  This is useful if the user has reset the
!    random number seed, for instance.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
!    X that we need to compute.  This starts off as 1:N, but is adjusted
!    if we have a saved value that can be immediately stored in X(1),
!    and so on.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r(n+1)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n

  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    call r8vec_uniform_01 ( 2, seed, r )

    x(x_hi_index) = &
             sqrt ( -2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) == 0 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( -2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( -2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

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
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
subroutine r8vec_uniform_ab ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    Each dimension ranges from A to B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_AB - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine theta_solve ( a, cl, m, theta )

!*****************************************************************************80
!
!! THETA_SOLVE solves a pair of transcendental equations.
!
!  Discussion:
!
!    The vector THETA returned by this function is needed in order to define
!    the terms in a Karhunen-Loeve expansion of a diffusion coefficient.
!
!    The two equations are:
!
!      1/CL - THETA * TAN ( A * THETA ) = 0
!      THETA - 1/CL * TAN ( A * THETA ) = 0
!
!    A and CL are taken to be positive.  Over each open interval 
!
!      ( n - 1/2 pi, n + 1/2 pi ) / A, for N = 0, 1, ...
!
!    the function TAN ( A * THETA ) monotonically rises from -oo to +00; 
!    therefore, it can be shown that there is one root of each equation 
!    in every interval of this form.  Moreover, because of the positivity
!    of A and CL, we can restrict our search to the interval 
!
!      [ n pi, n + 1/2 pi ) / A, for N = 0, 1, ...
!
!    This function computes K such roots, starting in the first interval,
!    finding those two roots, moving to the next interval, and so on, until
!    the requested number of roots have been found.  Odd index roots will
!    correspond to the first equation, and even index roots to the second.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Howard Elman, Darran Furnival,
!    Solving the Stochastic Steady-State Diffusion Problem Using Multigrid,
!    University of Maryland Department of Computer Science,
!    Technical Report TR-4786.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the "radius" of the domain, D = (-A,A)x(-A,A).
!    0 < A.
!
!    Input, real ( kind = 8 ) CL, the correlation length.
!    0 < CL.
!
!    Input, integer ( kind = 4 ) M, the number of values to compute.
!
!    Output, real ( kind = 8 ) THETA(M), the values of Theta.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a
  real ( kind = 8 ) bmatol
  real ( kind = 8 ) cl
  real ( kind = 8 ) eps
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) ftol
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_epsilon
  real ( kind = 8 ) theta(m)
  real ( kind = 8 ) xa
  real ( kind = 8 ) xa_init
  real ( kind = 8 ) xb
  real ( kind = 8 ) xb_init
  real ( kind = 8 ) xc

  k = 0
  theta(1:m) = 0.0D+00
!
!  [ XA_INIT, XB_INIT] = [ n * pi, n+1/2 pi ] / a, n = 0, 1, 2, ...
!
  xa_init = 0.0D+00
  xb_init = ( r8_pi / 2.0D+00 ) / a
  eps = r8_epsilon ( )

  do
!
!  Seek root of equation 1 in interval.
!
    if ( m <= k ) then
      exit
    end if

    k = k + 1
    xa = xa_init
    fa = 1.0D+00 / cl - xa * tan ( a * xa )
    ftol = eps * ( abs ( fa ) + 1.0D+00 )
    xb = xb_init
    fb = - fa
    fc = fa
    bmatol = 100.0D+00 * eps * ( abs ( xa ) + abs ( xb ) )

    do while ( bmatol < xb - xa )

      xc = ( xa + xb ) / 2.0D+00
      fc = 1.0D+00 / cl - xc * tan ( a * xc )

      if ( abs ( fc ) <= ftol ) then
        exit
      else if ( 0.0D+00 < fc ) then
        xa = xc
      else
        xb = xc
      end if

    end do

    theta(k) = xc
!
!  Seek root of equation 2 in interval.
!
    if ( m <= k ) then
      exit
    end if

    k = k + 1
!
!  In the first interval, we need to skip the zero root of equation 2.
!
    if ( k == 2 ) then

      k = k - 1

    else

      xa = xa_init
      fa = xa - tan ( a * xa ) / cl
      ftol = eps * ( abs ( fa ) + 1.0D+00 )
      xb = xb_init
      fb = - fa

      do while ( bmatol < xb - xa )

        xc = ( xa + xb ) / 2.0D+00
        fc = xc - tan ( a * xc ) / cl

        if ( abs ( fc ) <= ftol ) then
          exit
        else if ( 0.0D+00 < fc ) then
          xa = xc
        else
          xb = xc
        end if

      end do

      theta(k) = xc

    end if
!
!  Advance the interval.
!
    xa_init = xa_init + r8_pi / a
    xb_init = xb_init + r8_pi / a

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
