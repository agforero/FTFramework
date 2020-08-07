program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIANGLE.
!
!  Discussion:
!
!    This driver computes the interpolation of the Franke function
!    on the triangle T(U,V,W) with vertices U=(U1,U2)=(0,0), 
!    V=(V1,V2)=(1,0) and W=(W1,W2)=(0,1) (unit triangle) 
!    at the first family of Padua points. 
!
!    The degree of interpolation is DEG = 60 and the number of target 
!    points is NTG = NTG1 ** 2 - NTG1 + 1, NTG1 = 100.
!
!    The maps from the reference square [-1,1]^2 to the triangle
!    are SIGMA1 and SIGMA2 with inverses ISIGM1 and ISIGM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Local, integer ( kind = 4 ) DEGMAX, the maximum degree of interpolation.
!
!    Local, integer ( kind = 4 ) NPDMAX, the maximum number of Padua points
!    = (DEGMAX + 1) * (DEGMAX + 2) / 2.
!
!    Local, integer ( kind = 4 ) NTG1MX, the maximum value of the parameter determining 
!    the number of target points.
!
!    Local, integer ( kind = 4 ) NTGMAX, the maximum number of target points,
!    dependent on NTG1MX.
!
!    Local, integer ( kind = 4 ) DEG, the degree of interpolation.
!
!    Local, integer ( kind = 4 ) NTG1, the parameter determining the number 
!   of target points.
!
!    Local, integer ( kind = 4 ) NPD, the number of Padua points = (DEG + 1) * (DEG + 2) / 2.
!
!    Local, integer ( kind = 4 ) NTG, the number of target points, dependent on NTG1.
!
!    Local, real ( kind = 8 ) PD1(NPDMAX), the first coordinates of 
!    the Padua points.
!
!    Local, real ( kind = 8 ) PD2(NPDMAX), the second coordinates of the 
!    Padua points.
!
!    Local, real ( kind = 8 ) WPD(NPDMAX), the weights.
!
!    Local, real ( kind = 8 ) FPD(NPDMAX), the function at the Padua points.
!
!    Workspace, real ( kind = 8 ) RAUX1(DEGMAX+1)*(DEGMAX+2)).
!
!    Workspace, real ( kind = 8 ) RAUX2(DEGMAX+1)*(DEGMAX+2)).
!
!    Local, real ( kind = 8 ) C0(0:DEGMAX+1,0:DEGMAX+1), the coefficient matrix.
!
!    Local, real ( kind = 8 ) TG1(NTGMAX), the first coordinates of the 
!    target points.
!
!    Local, real ( kind = 8 ) TG2(NTGMAX), the second coordinates of the 
!    target points.
!
!    Local, real ( kind = 8 ) INTFTG(NTGMAX), the values of the 
!    interpolated function.
!
!    Local, real ( kind = 8 ) MAXERR, the maximum norm of the error at target 
!    points.
!
!    Local, real ( kind = 8 ) ESTERR, the estimated error.
!
  implicit none

  integer ( kind = 4 ), parameter :: degmax = 60
  integer ( kind = 4 ), parameter :: ntg1mx = 100

  integer ( kind = 4 ), parameter :: npdmax = ( degmax + 1 ) * ( degmax + 2 ) / 2
  integer ( kind = 4 ), parameter :: ntgmax = ntg1mx ** 2 - ntg1mx + 1

  real ( kind = 8 ) c0(0:degmax+1,0:degmax+1)
  integer ( kind = 4 ) deg
  real ( kind = 8 ) esterr
  integer ( kind = 4 ) family
  character * ( 255 ) filename
  real ( kind = 8 ) fmax
  real ( kind = 8 ) fmin
  real ( kind = 8 ) fpd(npdmax)
  real ( kind = 8 ) franke
  real ( kind = 8 ) fxy
  integer ( kind = 4 ) i
  real ( kind = 8 ) intftg(ntgmax)
  real ( kind = 8 ) isigm1
  real ( kind = 8 ) isigm2
  real ( kind = 8 ) ixy
  real ( kind = 8 ) maxdev
  real ( kind = 8 ) maxerr
  real ( kind = 8 ) mean
  integer ( kind = 4 ) npd
  integer ( kind = 4 ) ntg
  integer ( kind = 4 ) ntg1
  real ( kind = 8 ) pd1(npdmax)
  real ( kind = 8 ) pd2(npdmax)
  real ( kind = 8 ) pd2val
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) raux1((degmax+1)*(degmax+2))
  real ( kind = 8 ) raux2((degmax+1)*(degmax+2))
  real ( kind = 8 ) sigma1
  real ( kind = 8 ) sigma2
  real ( kind = 8 ) tg1(ntgmax)
  real ( kind = 8 ) tg2(ntgmax)
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) wpd(npdmax)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  u1 = 0.0D+00
  u2 = 0.0D+00
  v1 = 1.0D+00
  v2 = 0.0D+00
  w1 = 0.0D+00
  w2 = 1.0D+00
  family = 1
  deg = 60
  ntg1 = 100

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Interpolation of the Franke function'
  write ( *, '(a)' ) '  on the unit triangle T((0,0),(1,0),(0,1))'
  write ( *, '(a,i6)' ) '  at degree = ', deg

  if ( degmax < deg ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE - Fatal error!'
    write ( *, '(a)' ) '  DEGMAX < DEG.'
    write ( *, '(a,i6)' ) '  DEG =    ', deg
    write ( *, '(a,i6)' ) '  DEGMAX = ', degmax
    stop 1
  end if
!
!  Build the first family of Padua points in the square [-1,1]^2
!     
  call pdpts ( deg, pd1, pd2, wpd, npd )
!     
!  Compute the Franke function at Padua points mapped to T(U,V,W).
!   
  do i = 1, npd
    x = sigma1 ( pd1(i), pd2(i), u1, u2, v1, v2, w1, w2 )
    y = sigma2 ( pd1(i), pd2(i), u1, u2, v1, v2, w1, w2 )
    fpd(i) = franke ( x, y )
  end do
!
!  Write X, Y, F(X,Y) to a file.
!
  filename = 'triangle_fpd.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, npd
    x = sigma1 ( pd1(i), pd2(i), u1, u2, v1, v2, w1, w2 )
    y = sigma2 ( pd1(i), pd2(i), u1, u2, v1, v2, w1, w2 )
    write ( 10, '(3g14.6)' ) x, y, fpd(i)
  end do
  close ( unit = 10 )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Wrote F(x,y) at Padua points in "' &
    // trim ( filename ) // '".'
!     
!  Compute the matrix C0 of the coefficients in the bivariate
!  orthonormal Chebyshev basis
!     
  call padua2 ( deg, degmax, npd, wpd, fpd, raux1, raux2, c0, esterr )
!     
!  Build the set of target points on T(U,V,W)
!     
  call target ( u1, u2, v1, v2, w1, w2, ntg1, ntgmax, tg1, tg2, ntg )
!     
!  Evaluate the interpolant at the target points.
!    
  do i = 1, ntg
    x = isigm1 ( tg1(i), tg2(i), u1, u2, v1, v2, w1, w2 )
    y = isigm2 ( tg1(i), tg2(i), u1, u2, v1, v2, w1, w2 )
    intftg(i) = pd2val ( deg, degmax, c0, x, y )
  end do
!
!  Write the function value at target points to a file.
!
  filename = 'triangle_ftg.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, ntg
    write ( 10, '(3g14.6)' ) tg1(i), tg2(i), franke ( tg1(i), tg2(i) )
  end do
  close ( unit = 10 )
  write ( *, '(a)' ) '  Wrote F(x,y) at target points in "' &
    // trim ( filename ) // '".'
!
!  Write the interpolated function value at target points to a file.
!
  filename = 'triangle_itg.txt'
  open ( unit = 10, file = filename, status = 'replace' )
  do i = 1, ntg
    write ( 10, '(3g14.6)' ) tg1(i), tg2(i), intftg(i)
  end do
  close ( unit = 10 )
  write ( *, '(a)' ) '  Wrote I(F)(x,y) at target points in "' &
    // trim ( filename ) // '".'
!
!  Compute the error relative to the max deviation from the mean.
!     
  maxerr = 0.0D+00
  mean = 0.0D+00
  fmax = - r8_huge ( )
  fmin = + r8_huge ( )

  do i = 1, ntg
    fxy = franke ( tg1(i), tg2(i) )
    ixy = intftg(i)
    maxerr = max ( maxerr, abs ( fxy - ixy ) )
    mean = mean + fxy
    fmax = max ( fmax, fxy )
    fmin = min ( fmin, fxy )
  end do
 
  if ( fmax == fmin ) then
    maxdev = 1.0D+00
  else
    mean = mean / real ( ntg, kind = 8 )
    maxdev = max ( fmax - mean, mean - fmin )
  end if
!
!  Print error ratios.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,e10.4)' ) '  Estimated error:  ', esterr / maxdev
  write ( *, '(a,e10.4)' ) '  Actual error:     ', maxerr / maxdev
  write ( *, '(a,e10.4)' ) '  Expected error:   ', 0.1226D-09
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
function sigma1 ( t1, t2, u1, u2, v1, v2, w1, w2 )

!*****************************************************************************80
!
!! SIGMA1 maps first coordinate from square to triangle.
!
!  Discussion:
!
!    This function returns the first component of the map
!    from the square [-1,1]^2 to the triangle T(U,V,W). 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T1, T2, the coordinates of a point in the square.
!
!    Input, real ( kind = 8 ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Output, real ( kind = 8 ) SIGMA1, the X coordinate of the corresponding
!    point in the triangle.
!
  implicit none

  real ( kind = 8 ) sigma1
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2

  sigma1 = ( v1 - u1 ) * ( 1.0D+00 + t1 ) &
    * ( 1.0D+00 - t2 ) / 4.0D+00 &
    + ( w1 - u1 ) * ( 1.0D+00 + t2 ) / 2.0D+00 + u1

  return
end
function isigm1 ( sigma1, sigma2, u1, u2, v1, v2, w1, w2 )

!*****************************************************************************80
!
!! ISIGM1 maps first coordinate from triangle to the square.
!
!  Discussion:
!
!    This functions returns the first component of the map
!    from the the triangle T(U,V,W) to the square [-1,1]^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SIGMA1, SIGMA2, the coordinates of a point 
!    in the triangle.
!
!    Input, real ( kind = 8 ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Output, real ( kind = 8 ) ISIGM1, the X coordinate of the corresponding
!    point in the square.
!
  implicit none

  real ( kind = 8 ) isigm1
  real ( kind = 8 ) rho1
  real ( kind = 8 ) rho2
  real ( kind = 8 ) sigma1
  real ( kind = 8 ) sigma2
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2

  rho1 = ( sigma1 * ( w2 - u2 ) - sigma2 * ( w1 - u1 ) &
    + ( w1 - u1 ) * u2 - ( w2 - u2 ) * u1 ) / &
    ( ( v1 - u1 ) * ( w2 - u2 ) - ( v2 - u2 ) * ( w1 - u1 ) )

  rho2 = ( sigma1 * ( v2 - u2 ) - sigma2 * ( v1 - u1 ) &
    + ( v1 - u1 ) * u2 - ( v2 - u2 ) * u1 ) / &
    ( ( w1 - u1 ) * ( v2 - u2 ) - ( w2 - u2 ) * ( v1 - u1 ) )

  if ( rho2 == 1.0D+00 ) then
    isigm1 = 0.0D+00
  else
    isigm1 = 2.0D+00 * rho1 / ( 1.0D+00 - rho2 ) - 1.0D+00
  end if  
 
  return
end
function sigma2 ( t1, t2, u1, u2, v1, v2, w1, w2 )

!*****************************************************************************80
!
!! SIGMA2 maps the second coordinate from square to triangle.
!
!  Discussion:
!
!    This functions returns the second component of the map
!    from the square [-1,1]^2 to the triangle T(U,V,W).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T1, T2, the coordinates of a point in the square.
!
!    Input, real ( kind = 8 ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Output, real ( kind = 8 ) SIGMA2, the Y coordinate of the corresponding
!    point in the triangle.
!
  implicit none

  real ( kind = 8 ) sigma2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2

  sigma2 = ( v2 - u2 ) * ( 1.0D+00 + t1 ) &
    * ( 1.0D+00 - t2 ) / 4.0D+00 + ( w2 - u2 ) &
    * ( 1.0D+00 + t2 ) / 2.0D+00 + u2

  return
end
function isigm2 ( sigma1, sigma2, u1, u2, v1, v2, w1, w2 )

!*****************************************************************************80
!
!! ISIGM2 maps second coordinate from triangle to the square.
!
!  Discussion:
!
!    This functions returns the second component of the map
!    from the the triangle T(U,V,W) to the square [-1,1]^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) SIGMA1, SIGMA2, the coordinates of a point 
!    in the ellipse.
!
!    Input, real ( kind = 8 ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Output, real ( kind = 8 ) ISIGM2, the Y coordinate of the corresponding
!    point in the triangle.
!
  implicit none

  real ( kind = 8 ) isigm2
  real ( kind = 8 ) rho2
  real ( kind = 8 ) sigma1
  real ( kind = 8 ) sigma2
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2

  rho2 = ( sigma1 * ( v2 - u2 ) - &
    sigma2 * ( v1 - u1) + ( v1 - u1 ) * u2 - ( v2 - u2 ) * u1 ) / &
    ( ( w1 - u1 ) * ( v2 - u2 ) - ( w2 - u2 ) * ( v1 - u1 ) )

  if ( rho2 == 1.0D+00 ) then
    isigm2 = 1.0D+00
  else
    isigm2 = 2.0D+00 * rho2 - 1.0D+00
  end if 
  
  return
end
subroutine target ( u1, u2, v1, v2, w1, w2, ntg1, ntgmax, tg1, tg2, ntg )

!*****************************************************************************80
!
!! TARGET returns the target points on the triangle.
!
!  Discussion:
!
!    Target points on the triangle T(U,V,W).
!    The number of target points is NTG = NTG1^2 - NGT1 + 1.
!
!  Licensing:
!
!    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
!    Marco Vianello.
!    This FORTRAN90 version by John Burkardt.
!  
!  Modified:
!
!    12 February 2014
!
!  Author:
!
!    Marco Caliari, Stefano De Marchi, Marco Vianello
!
!  Reference:
!
!    Marco Caliari, Stefano de Marchi, Marco Vianello,
!    Algorithm 886:
!    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
!    ACM Transactions on Mathematical Software,
!    Volume 35, Number 3, October 2008, Article 21, 11 pages.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U1, U2, V1, V2, W1, W2, the coordinates of the 
!    vertices of the triangle.
!
!    Input, integer ( kind = 4 ) NTG1, a parameter determining the number 
!    of target points
!
!    Input, integer ( kind = 4 ) NTGMAX, the maximum number of target points.
!
!    Output, real ( kind = 8 ) TG1(NTG), TG2(NTG), the X and Y coordinates
!    of the target points.
!
!    Output, integer ( kind = 4 ) NTG, the number of target points computed.
!
  implicit none

  integer ( kind = 4 ) ntgmax

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntg
  integer ( kind = 4 ) ntg1
  real ( kind = 8 ) tg1(ntgmax)
  real ( kind = 8 ) tg2(ntgmax)
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2

  if ( ntg1 < 2 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TARGET - Fatal error!'
    write ( *, '(a)' ) '  NTG1 < 2'
    write ( *, '(a,i4)' ) '  NTG1 = ', ntg1
    stop 1
  end if  
   
  if ( ntgmax < ntg1 ** 2 - ntg1 + 1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'TARGET - Fatal error!'
    write ( *, '(a)' ) '  NTGMAX < NTG1 * NTG1 - NTG1 + 1.'
    write ( *, '(a,i4)' ) '  NTG1 = ', ntg1
    write ( *, '(a,i4)' ) '  NTGMAX = ', ntgmax
    stop 1
  end if    

  ntg = 0
  do i = 1, ntg1 - 1
    do j = 1, ntg1

      ntg = ntg + 1

      tg1(ntg) = ( v1 - u1 ) * real ( i - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) &
        + ( w1 - u1 ) * ( real ( j - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) ) &
        * ( 1.0D+00 - real ( i - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) ) + u1

      tg2(ntg) = ( v2 - u2 ) * real ( i - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) &
        + ( w2 - u2 ) * ( real ( j - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) ) &
        * ( 1.0D+00 - real ( i - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) ) + u2

    end do
  end do

  i = ntg1
  j = 1
  ntg = ntg + 1

  tg1(ntg) = ( v1 - u1 ) * real ( i - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) &
    + ( w1 - u1 ) * ( real ( j - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) ) &
    * ( 1.0D+00 - real ( i - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) ) + u1

  tg2(ntg) = ( v2 - u2 ) * real ( i - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) &
    + ( w2 - u2 ) * ( real ( j - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) ) &
    * ( 1.0D+00 - real ( i - 1, kind = 8 ) / real ( ntg1 - 1, kind = 8 ) ) + u2

  return
end
