MODULE poly_zeroes
IMPLICIT NONE
PRIVATE :: aberth, newton, start, cnvex, cmerge, left, right, ctest
PUBLIC  :: polzeros
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)

CONTAINS

!************************************************************************
!    NUMERICAL COMPUTATION OF THE ROOTS OF A POLYNOMIAL HAVING          *
!        COMPLEX COEFFICIENTS, BASED ON ABERTH'S METHOD.                *
!                      Version 1.4, June   1996                         *
!    (D. Bini, Dipartimento di Matematica, Universita' di Pisa)         *
!                         (bini@dm.unipi.it)                            *
!
!  ***************************************************************************
!  * All the software  contained in this library  is protected by copyright. *
!  * Permission  to use, copy, modify, and  distribute this software for any *
!  * purpose without fee is hereby granted, provided that this entire notice *
!  * is included  in all copies  of any software which is or includes a copy *
!  * or modification  of this software  and in all copies  of the supporting *
!  * documentation for such software.                                        *
!  ***************************************************************************
!  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
!  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
!  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
!  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
!  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
!  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
!  ***************************************************************************
!  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
!  * ABOVE STATEMENT.                                                        *
!  ***************************************************************************
!
!   AUTHOR:
!
!       DARIO ANDREA BINI
!       UNIVERSITY OF PISA, ITALY
!       E-MAIL: bini@dm.unipi.it
!
!   REFERENCE:
!
!    -  NUMERICAL COMPUTATION OF POLYNOMIAL ZEROS BY MEANS OF
!       ABERTH'S METHOD
!       NUMERICAL ALGORITHMS, 13 (1996), PP. 179-200
!
!   This version, which is compatible with Lahey's free ELF90 compiler
!   by Alan Miller, CSIRO Mathematical & Information Sciences,
!   Private Bag 10, Clayton South MDC, Victoria, Australia 3169.
!   amiller @ bigpond.net.au     http://users.bigpond.net.au/amiller
!   Latest revision of ELF90 version - 5 May 1997
!
!************************************************************************
! Work performed under the support of the ESPRIT BRA project 6846 POSSO *
!************************************************************************
!**********         SUBROUTINES AND FUNCTIONS                 ***********
!************************************************************************
!  The following routines/functions are listed:                         *
!  POLZEROS  :  computes polynomial roots by means of Aberth's method   *
!    ABERTH  :  computes the Aberth correction                          *
!    NEWTON  :  computes p(x)/p'(x) by means of Ruffini-Horner's rule   *
!    START   :  Selects N starting points by means of Rouche's theorem  *
!    CNVEX   :  Computes the convex hull, used by START                 *
!    CMERGE  :  Used by CNVEX                                           *
!    LEFT    :  Used by CMERGE                                          *
!    RIGHT   :  Used by CMERGE                                          *
!    CTEST   :  Convexity test, Used by CMERGE                          *
!************************************************************************
!                                                                       *
!                                                                       *
!************************************************************************
!********************** SUBROUTINE POLZEROS *****************************
!************************************************************************
!                        GENERAL COMMENTS                               *
!************************************************************************
!  This routine approximates the roots of   the  polynomial             *
!  p(x)=a(n+1)x^n+a(n)x^(n-1)+...+a(1), a(j)=cr(j)+I ci(j), I**2=-1,    *
!  where a(1) and a(n+1) are nonzero.                                   *
!  The coefficients are complex*16 numbers. The routine is fast, robust *
!  against overflow, and allows to deal with polynomials of any degree. *
!  Overflow situations are very unlikely and may occurr if there exist  *
!  simultaneously coefficients of moduli close to BIG and close to      *
!  SMALL, i.e., the greatest and the smallest positive real*8 numbers,  *
!  respectively. In this limit situation the program outputs a warning  *
!  message. The computation can be speeded up by performing some side   *
!  computations in single precision, thus slightly reducing the         *
!  robustness of the program (see the comments in the routine ABERTH).  *
!  Besides a set of approximations to the roots, the program delivers a *
!  set of a-posteriori error bounds which are guaranteed in the most    *
!  part of cases. In the situation where underflow does not allow to    *
!  compute a guaranteed bound, the program outputs a warning message    *
!  and sets the bound to 0. In the situation where the root cannot be   *
!  represented as a complex*16 number the error bound is set to -1.     *
!************************************************************************
!  The computation is performed by means of Aberth's method             *
!  according to the formula                                             *
!           x(i)=x(i)-newt/(1-newt*abcorr), i=1,...,n             (1)   *
!  where newt=p(x(i))/p'(x(i)) is the Newton correction and abcorr=     *
!  =1/(x(i)-x(1))+...+1/(x(i)-x(i-1))+1/(x(i)-x(i+1))+...+1/(x(i)-x(n)) *
!  is the Aberth correction to the Newton method.                       *
!************************************************************************
!  The value of the Newton correction is computed by means of the       *
!  synthetic division algorithm (Ruffini-Horner's rule) if |x|<=1,      *
!  otherwise the following more robust (with respect to overflow)       *
!  formula is applied:                                                  *
!                    newt=1/(n*y-y**2 R'(y)/R(y))                 (2)   *
!  where                                                                *
!                    y=1/x                                              *
!                    R(y)=a(1)*y**n+...+a(n)*y+a(n+1)            (2')   *
!  This computation is performed by the routine NEWTON.                 *
!************************************************************************
!  The starting approximations are complex numbers that are             *
!  equispaced on circles of suitable radii. The radius of each          *
!  circle, as well as the number of roots on each circle and the        *
!  number of circles, is determined by applying Rouche's theorem        *
!  to the functions a(k+1)*x**k and p(x)-a(k+1)*x**k, k=0,...,n.        *
!  This computation is performed by the routine START.                  *
!************************************************************************
!                              STOP CONDITION                           *
!************************************************************************
! If the condition                                                      *
!                     |p(x(j))|<EPS s(|x(j)|)                      (3)  *
! is satisfied,    where      s(x)=s(1)+x*s(2)+...+x**n * s(n+1),       *
! s(i)=|a(i)|*(1+3.8*(i-1)),  EPS is the machine precision (EPS=2**-53  *
! for the IEEE arithmetic), then the approximation x(j) is not updated  *
! and the subsequent iterations (1)  for i=j are skipped.               *
! The program stops if the condition (3) is satisfied for j=1,...,n,    *
! or if the maximum number NITMAX of  iterations   has   been reached.  *
! The condition (3) is motivated by a backward rounding error analysis  *
! of the Ruffini-Horner rule, moreover the condition (3) guarantees     *
! that the computed approximation x(j) is an exact root of a slightly   *
! perturbed polynomial.                                                 *
!************************************************************************
!             INCLUSION DISKS, A-POSTERIORI ERROR BOUNDS                *
!************************************************************************
! For each approximation x of a root, an a-posteriori absolute error    *
! bound r is computed according to the formula                          *
!                   r=n(|p(x)|+EPS s(|x|))/|p'(x)|                 (4)  *
! This provides an inclusion disk of center x and radius r containing a *
! root.                                                                 *
!************************************************************************
!************************************************************************
!*************       MEANING OF THE INPUT VARIABLES         *************
!************************************************************************
!************************************************************************
!                                                                       *
!  -- N     : degree of the polynomial.                                 *
!  -- POLY  : complex vector of N+1 components, POLY(i) is the          *
!           coefficient of x**(i-1), i=1,...,N+1 of the polynomial p(x) *
!  -- EPS   : machine precision of the floating point arithmetic used   *
!            by the computer, EPS=2**(-53)  for the IEEE standard.      *
!  -- BIG   : the max real*8, BIG=2**1023 for the IEEE standard.        *
!  -- SMALL : the min positive real*8, SMALL=2**(-1074) for the IEEE.   *
!  -- NITMAX: the max number of allowed iterations.                     *
!************************************************************************
!************************************************************************
!*************      MEANING OF THE OUTPUT VARIABLES         *************
!************************************************************************
!************************************************************************
!  ROOT   : complex vector of N components, containing the              *
!           approximations to the roots of p(x).                        *
!  RADIUS : real vector of N components, containing the error bounds to *
!           the approximations of the roots, i.e. the disk of center    *
!           ROOT(i) and radius RADIUS(i) contains a root of p(x), for   *
!           i=1,...,N. RADIUS(i) is set to -1 if the corresponding root *
!           cannot be represented as floating point due to overflow or  *
!           underflow.                                                  *
!  ERR    : vector of (N+1) components detecting an error condition;    *
!           ERR(j)=.TRUE. if after NITMAX iterations the stop condition *
!                         (3) is not satisfied for x(j)=ROOT(j);        *
!           ERR(j)=.FALSE.  otherwise, i.e., the root is reliable,      *
!                         i.e., it can be viewed as an exact root of a  *
!                         slightly perturbed polynomial.                *
!           The vector ERR is used also in the routine convex hull for  *
!           storing the abscissae of the vertices of the convex hull.   *
!           ERR(N+1) is only used for the convex hull.                  *
!  ITER   : number of iterations peformed.                              *
!************************************************************************
!************************************************************************
!************    MEANING OF THE AUXILIARY VARIABLES         *************
!************************************************************************
!************************************************************************
!  APOLY  : real vector of N+1 components used to store the moduli of   *
!           the coefficients of p(x) and the coefficients of s(x) used  *
!           to test the stop condition (3).                             *
!  APOLYR : real vector of N+1 components used to test the stop         *
!           condition                                                   *
!************************************************************************
!*****         WARNING:   2 is the output unit                    *******
!************************************************************************
SUBROUTINE polzeros (n, poly, eps, big, small, nitmax, root, radius, err, &
                     iter)
IMPLICIT NONE
INTEGER, INTENT(IN)            :: n, nitmax
INTEGER, INTENT(OUT)           :: iter
COMPLEX (KIND=dp), INTENT(IN)  :: poly(:)
COMPLEX (KIND=dp), INTENT(OUT) :: root(:)
REAL (KIND=dp), INTENT(OUT)    :: radius(:)
REAL (KIND=dp), INTENT(IN)     :: eps, small, big
LOGICAL, INTENT(OUT)           :: err(:)

! Local variables
INTEGER                   :: i, nzeros
COMPLEX (KIND=dp)         :: corr, abcorr
REAL (KIND=dp)            :: amax, apoly(n+1), apolyr(n+1)
REAL (KIND=dp), PARAMETER :: zero = 0.0_dp

! Check consistency of data
IF (ABS(poly(n+1)) == zero)THEN
  WRITE(*,*)'Inconsistent data: the leading coefficient is zero'
  STOP
END IF
IF (ABS(poly(1)) == zero) THEN
  WRITE(*,*)'The constant term is zero: deflate the polynomial'
  STOP
END IF

! Compute the moduli of the coefficients
amax = 0
DO i = 1, n+1
  apoly(i) = ABS(poly(i))
  amax = MAX(amax, apoly(i))
  apolyr(i) = apoly(i)
END DO
IF((amax) >= (big/(n+1))) THEN
  WRITE(*,*)'WARNING: COEFFICIENTS TOO BIG, OVERFLOW IS LIKELY'
  WRITE(2,*)'WARNING: COEFFICIENTS TOO BIG, OVERFLOW IS LIKELY'
END IF
! Initialize
DO i = 1, n
  radius(i) = zero
  err(i) = .true.
END DO

! Select the starting points
CALL start(n, apolyr, root, radius, nzeros, small, big, err)

! Compute the coefficients of the backward-error polynomial
DO i = 1, n+1
  apolyr(n-i+2) = eps*apoly(i)*(3.8*(n-i+1) + 1)
  apoly(i) = eps*apoly(i)*(3.8*(i-1) + 1)
END DO
IF((apoly(1) == 0).OR.(apoly(n+1) == 0)) THEN
  WRITE(*,*)'WARNING: THE COMPUTATION OF SOME INCLUSION RADIUS'
  WRITE(*,*)'MAY FAIL. THIS IS REPORTED BY RADIUS = 0'
  WRITE(2,*)'WARNING: THE COMPUTATION OF SOME INCLUSION RADIUS'
  WRITE(2,*)'MAY FAIL. THIS IS REPORTED BY RADIUS = 0'
END IF
DO i = 1,n
  err(i) = .true.
  IF(radius(i) == -1) err(i) = .false.
END DO

! Starts Aberth's iterations
DO iter = 1, nitmax
  DO i = 1, n
    IF (err(i)) THEN
      CALL newton(n, poly, apoly, apolyr, root(i), small, radius(i), corr, &
                  err(i))
      IF (err(i)) THEN
        CALL aberth(n, i, root, abcorr)
        root(i) = root(i) - corr/(1 - corr*abcorr)
      ELSE
        nzeros = nzeros + 1
        IF (nzeros == n) RETURN
      END IF
    END IF
  END DO
END DO
RETURN
END SUBROUTINE polzeros




!************************************************************************
!                             SUBROUTINE NEWTON                         *
!************************************************************************
! Compute  the Newton's correction, the inclusion radius (4) and checks *
! the stop condition (3)                                                *
!************************************************************************
! Input variables:                                                      *
!     N     : degree of the polynomial p(x)                             *
!     POLY  : coefficients of the polynomial p(x)                       *
!     APOLY : upper bounds on the backward perturbations on the         *
!             coefficients of p(x) when applying Ruffini-Horner's rule  *
!     APOLYR: upper bounds on the backward perturbations on the         *
!             coefficients of p(x) when applying (2), (2')              *
!     Z     : value at which the Newton correction is computed          *
!     SMALL : the min positive REAL (KIND=dp) ::, SMALL=2**(-1074) for the IEEE.   *
!************************************************************************
! Output variables:                                                     *
!     RADIUS: upper bound to the distance of Z from the closest root of *
!             the polynomial computed according to (4).                 *
!     CORR  : Newton's correction                                       *
!     AGAIN : this variable is .true. if the computed value p(z) is     *
!             reliable, i.e., (3) is not satisfied in Z. AGAIN is       *
!             .false., otherwise.                                       *
!************************************************************************
SUBROUTINE newton(n, poly, apoly, apolyr, z, small, radius, corr, again)
IMPLICIT NONE
INTEGER, INTENT(IN)            :: n
COMPLEX (KIND=dp), INTENT(IN)  :: poly(:), z
COMPLEX (KIND=dp), INTENT(OUT) :: corr
REAL (KIND=dp), INTENT(IN)     :: apoly(:), apolyr(:), small
REAL (KIND=dp), INTENT(OUT)    :: radius
LOGICAL, INTENT(OUT)           :: again

! Local variables
INTEGER           :: i
COMPLEX (KIND=dp) :: p, p1, zi, den, ppsp
REAL (KIND=dp)    :: ap, az, azi, absp

az = ABS(z)
! If |z|<=1 then apply Ruffini-Horner's rule for p(z)/p'(z)
! and for the computation of the inclusion radius
IF(az <= 1)THEN
  p = poly(n+1)
  ap = apoly(n+1)
  p1 = p
  DO i = n, 2, -1
    p = p*z + poly(i)
    p1 = p1*z + p
    ap = ap*az + apoly(i)
  END DO
  p = p*z + poly(1)
  ap = ap*az + apoly(1)
  corr = p/p1
  absp = ABS(p)
  ap = ap
  again = (absp > (small + ap))
  IF(.NOT.again) radius = n*(absp + ap)/ABS(p1)
  RETURN
ELSE
! If |z| > 1 then apply Ruffini-Horner's rule to the reversed polynomial
! and use formula (2) for p(z)/p'(z). Analogously do for the inclusion radius.
  zi = 1/z
  azi = 1/az
  p = poly(1)
  p1 = p
  ap = apolyr(n+1)
  DO i = n, 2, -1
    p = p*zi + poly(n-i+2)
    p1 = p1*zi + p
    ap = ap*azi + apolyr(i)
  END DO
  p = p*zi + poly(n+1)
  ap = ap*azi + apolyr(1)
  absp = ABS(p)
  again = (absp > (small+ap))
  ppsp = (p*z)/p1
  den = n*ppsp - 1
  corr = z*(ppsp/den)
  IF(again)RETURN
  radius = ABS(ppsp) + (ap*az)/ABS(p1)
  radius = n*radius/ABS(den)
  radius = radius*az
END IF
RETURN
END SUBROUTINE newton



!************************************************************************
!                             SUBROUTINE ABERTH                         *
!************************************************************************
! Compute  the Aberth correction. To save time, the reciprocation of    *
! ROOT(J)-ROOT(I) could be performed in single precision (complex*8)    *
! In principle this might cause overflow if both ROOT(J) and ROOT(I)    *
! have too small moduli.                                                *
!************************************************************************
! Input variables:                                                      *
!     N     : degree of the polynomial                                  *
!     ROOT  : vector containing the current approximations to the roots *
!     J     : index of the component of ROOT with respect to which the  *
!             Aberth correction is computed                             *
!************************************************************************
! Output variable:                                                      *
!     ABCORR: Aberth's correction (compare (1))                         *
!************************************************************************
SUBROUTINE aberth(n, j, root, abcorr)
IMPLICIT NONE
INTEGER, INTENT(IN)            :: n, j
COMPLEX (KIND=dp), INTENT(IN)  :: root(:)
COMPLEX (KIND=dp), INTENT(OUT) :: abcorr

! Local variables
INTEGER                   :: i
COMPLEX (KIND=dp)         :: z, zj
REAL (KIND=dp), PARAMETER :: zero = 0.0_dp

abcorr = CMPLX(zero, zero, KIND=dp)
zj = root(j)
DO i = 1, j-1
  z = zj - root(i)
  abcorr = abcorr + 1/z
END DO
DO i = j+1, n
  z = zj - root(i)
  abcorr = abcorr + 1/z
END DO
RETURN
END SUBROUTINE aberth



!************************************************************************
!                             SUBROUTINE START                          *
!************************************************************************
! Compute  the starting approximations of the roots                     *
!************************************************************************
! Input variables:                                                      *
!     N     :  number of the coefficients of the polynomial             *
!     A     :  moduli of the coefficients of the polynomial             *
!     SMALL : the min positive REAL (KIND=dp) ::, SMALL=2**(-1074) for the IEEE.   *
!     BIG   : the max REAL (KIND=dp) ::, BIG=2**1023 for the IEEE standard.        *
! Output variables:                                                     *
!     Y     :  starting approximations                                  *
!     RADIUS:  if a component is -1 then the corresponding root has a   *
!              too big or too small modulus in order to be represented  *
!              as double float with no overflow/underflow               *
!     NZ    :  number of roots which cannot be represented without      *
!              overflow/underflow                                       *
! Auxiliary variables:                                                  *
!     H     :  needed for the computation of the convex hull            *
!************************************************************************
! This routines selects starting approximations along circles center at *
! 0 and having suitable radii. The computation of the number of circles *
! and of the corresponding radii is performed by computing the upper    *
! convex hull of the set (i,log(A(i))), i=1,...,n+1.                    *
!************************************************************************
SUBROUTINE start(n, a, y, radius, nz, small, big, h)
IMPLICIT NONE
INTEGER, INTENT(IN)            :: n
INTEGER, INTENT(OUT)           :: nz
LOGICAL, INTENT(OUT)           :: h(:)
COMPLEX (KIND=dp), INTENT(OUT) :: y(:)
REAL (KIND=dp), INTENT(IN)     :: small, big
REAL (KIND=dp), INTENT(IN OUT) :: a(:)
REAL (KIND=dp), INTENT(OUT)    :: radius(:)

! Local variables
INTEGER                   :: i, iold, nzeros, j, jj
REAL (KIND=dp)            :: r, th, ang, temp, xsmall, xbig
REAL (KIND=dp), PARAMETER :: pi2 = 6.2831853071796, sigma = 0.7

xsmall = LOG(small)
xbig = LOG(big)
nz = 0
! Compute the logarithm A(I) of the moduli of the coefficients of
! the polynomial and then the upper covex hull of the set (A(I),I)
DO i = 1, n+1
  IF(a(i) /= 0) THEN
    a(i) = LOG(a(i))
  ELSE
    a(i) = -1.d30
  END IF
END DO
CALL cnvex(n+1, a, h)
! Given the upper convex hull of the set (A(I),I) compute the moduli
! of the starting approximations by means of Rouche's theorem
iold = 1
th = pi2/n
DO i = 2, n+1
  IF (h(i)) THEN
    nzeros = i - iold
    temp = (a(iold) - a(i))/nzeros
! Check if the modulus is too small
    IF((temp < -xbig).AND.(temp >= xsmall))THEN
      WRITE(*,*)'WARNING:',nzeros,' ZERO(S) ARE TOO SMALL TO'
      WRITE(*,*)'REPRESENT THEIR INVERSES AS COMPLEX (KIND=dp) ::, THEY'
      WRITE(*,*)'ARE REPLACED BY SMALL NUMBERS, THE CORRESPONDING'
      WRITE(*,*)'RADII ARE SET TO -1'
      WRITE(2,*)'WARNING:',nzeros,' ZERO(S) ARE TOO SMALL TO '
      WRITE(2,*)'REPRESENT THEIR INVERSES AS COMPLEX (KIND=dp) ::, THEY'
      WRITE(2,*)'ARE REPLACED BY SMALL NUMBERS, THE CORRESPONDING'
      WRITE(2,*)'RADII ARE SET TO -1'
      nz = nz + nzeros
      r = 1.0D0/big
    END IF
    IF(temp < xsmall)THEN
      nz = nz + nzeros
      WRITE(*,*)'WARNING: ',nzeros,' ZERO(S) ARE TOO SMALL TO BE'
      WRITE(*,*)'REPRESENTED AS COMPLEX (KIND=dp) ::, THEY ARE SET TO 0'
      WRITE(*,*)'THE CORRESPONDING RADII ARE SET TO -1'
      WRITE(2,*)'WARNING: ',nzeros,' ZERO(S) ARE TOO SMALL TO BE'
      WRITE(2,*)'REPRESENTED AS COMPLEX (KIND=dp) ::, THEY ARE SET 0'
      WRITE(2,*)'THE CORRESPONDING RADII ARE SET TO -1'
    END IF
! Check if the modulus is too big
    IF(temp > xbig)THEN
      r = big
      nz = nz + nzeros
      WRITE(*,*)'WARNING: ', nzeros, ' ZEROS(S) ARE TOO BIG TO BE'
      WRITE(*,*)'REPRESENTED AS COMPLEX (KIND=dp) ::,'
      WRITE(*,*)'THE CORRESPONDING RADII ARE SET TO -1'
      WRITE(2,*)'WARNING: ',nzeros, ' ZERO(S) ARE TOO BIG TO BE'
      WRITE(2,*)'REPRESENTED AS COMPLEX (KIND=dp) ::,'
      WRITE(2,*)'THE CORRESPONDING RADII ARE SET TO -1'
    END IF
    IF((temp <= xbig).AND.(temp > MAX(-xbig, xsmall)))THEN
      r = EXP(temp)
    END IF
! Compute NZEROS approximations equally distributed in the disk of
! radius R
    ang = pi2/nzeros
    DO j = iold, i-1
      jj = j-iold+1
      IF((r <= (1.0D0/big)).OR.(r == big)) radius(j) = -1
      y(j) = r*(COS(ang*jj + th*i + sigma) + (0,1)*SIN(ang*jj + th*i + sigma))
    END DO
    iold = i
  END IF
END DO
RETURN
END SUBROUTINE start


!************************************************************************
!                             SUBROUTINE CNVEX                          *
!************************************************************************
! Compute  the upper convex hull of the set (i,a(i)), i.e., the set of  *
! vertices (i_k,a(i_k)), k=1,2,...,m, such that the points (i,a(i)) lie *
! below the straight lines passing through two consecutive vertices.    *
! The abscissae of the vertices of the convex hull equal the indices of *
! the TRUE  components of the logical output vector H.                  *
! The used method requires O(nlog n) comparisons and is based on a      *
! divide-and-conquer technique. Once the upper convex hull of two       *
! contiguous sets  (say, {(1,a(1)),(2,a(2)),...,(k,a(k))} and           *
! {(k,a(k)), (k+1,a(k+1)),...,(q,a(q))}) have been computed, then       *
! the upper convex hull of their union is provided by the subroutine    *
! CMERGE. The program starts with sets made up by two consecutive       *
! points, which trivially constitute a convex hull, then obtains sets   *
! of 3,5,9... points,  up to  arrive at the entire set.                 *
! The program uses the subroutine  CMERGE; the subroutine CMERGE uses   *
! the subroutines LEFT, RIGHT and CTEST. The latter tests the convexity *
! of the angle formed by the points (i,a(i)), (j,a(j)), (k,a(k)) in the *
! vertex (j,a(j)) up to within a given tolerance TOLER, where i<j<k.    *
!************************************************************************
SUBROUTINE cnvex(n, a, h)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: n
LOGICAL, INTENT(OUT)       :: h(:)
REAL (KIND=dp), INTENT(IN) :: a(:)

! Local variables
INTEGER :: i, j, k, m, nj, jc

h(1:n) = .true.

! compute K such that N-2 <= 2**K < N-1
k = INT(LOG(n-2.0D0)/LOG(2.0D0))
IF(2**(k+1) <= (n-2)) k = k+1

! For each M=1,2,4,8,...,2**K, consider the NJ pairs of consecutive
! sets made up by M+1 points having the common vertex
! (JC,A(JC)), where JC=M*(2*J+1)+1 and J=0,...,NJ,
! NJ = MAX(0, INT((N-2-M)/(M+M))).
! Compute the upper convex hull of their union by means of subroutine CMERGE
m = 1
DO i = 0, k
  nj = MAX(0, INT((n-2-m)/(m+m)))
  DO j = 0, nj
    jc = (j+j+1)*m+1
    CALL cmerge(n, a, jc, m, h)
  END DO
  m = m+m
END DO
RETURN
END SUBROUTINE cnvex



!************************************************************************
!                             SUBROUTINE LEFT                           *
!************************************************************************
! Given as input the integer I and the vector H of logical, compute the *
! the maximum integer IL such that IL<I and H(IL) is TRUE.              *
!************************************************************************
! Input variables:                                                      *
!     H   : vector of logical                                           *
!     I   : integer                                                     *
!************************************************************************
! Output variable:                                                      *
!     IL  : maximum integer such that IL<I, H(IL)=.TRUE.                *
!************************************************************************
SUBROUTINE left(h, i, il)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: i
INTEGER, INTENT(OUT) :: il
LOGICAL, INTENT(IN)  :: h(:)

DO il = i-1, 0, -1
  IF (h(il)) RETURN
END DO
RETURN
END SUBROUTINE left



!************************************************************************
!                             SUBROUTINE RIGHT                          *
!************************************************************************
!************************************************************************
! Given as input the integer I and the vector H of logical, compute the *
! the minimum integer IR such that IR>I and H(IL) is TRUE.              *
!************************************************************************
!************************************************************************
! Input variables:                                                      *
!     N   : length of the vector H                                      *
!     H   : vector of logical                                           *
!     I   : integer                                                     *
!************************************************************************
! Output variable:                                                      *
!     IR  : minimum integer such that IR>I, H(IR)=.TRUE.                *
!************************************************************************
SUBROUTINE right(n, h, i, ir)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: n, i
INTEGER, INTENT(OUT) :: ir
LOGICAL, INTENT(IN)  :: h(:)

DO ir = i+1, n
  IF (h(ir)) RETURN
END DO
RETURN
END SUBROUTINE right



!************************************************************************
!                             SUBROUTINE CMERGE                         *
!************************************************************************
! Given the upper convex hulls of two consecutive sets of pairs         *
! (j,A(j)), compute the upper convex hull of their union                *
!************************************************************************
! Input variables:                                                      *
!     N    : length of the vector A                                     *
!     A    : vector defining the points (j,A(j))                        *
!     I    : abscissa of the common vertex of the two sets              *
!     M    : the number of elements of each set is M+1                  *
!************************************************************************
! Input/Output variable:                                                *
!     H    : vector defining the vertices of the convex hull, i.e.,     *
!            H(j) is .TRUE. if (j,A(j)) is a vertex of the convex hull  *
!            This vector is used also as output.                        *
!************************************************************************
SUBROUTINE cmerge(n, a, i, m, h)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: n, m, i
LOGICAL, INTENT(IN OUT)    :: h(:)
REAL (KIND=dp), INTENT(IN) :: a(:)

! Local variables
INTEGER :: ir, il, irr, ill
LOGICAL :: tstl, tstr

! at the left and the right of the common vertex (I,A(I)) determine
! the abscissae IL,IR, of the closest vertices of the upper convex
! hull of the left and right sets, respectively
CALL left(h, i, il)
CALL right(n, h, i, ir)

! check the convexity of the angle formed by IL,I,IR
IF (ctest(a, il, i, ir)) THEN
  RETURN
ELSE
! continue the search of a pair of vertices in the left and right
! sets which yield the upper convex hull
  h(i) = .false.
  DO
    IF (il == (i-m)) THEN
      tstl = .true.
    ELSE
      CALL left(h, il, ill)
      tstl = ctest(a, ill, il, ir)
    END IF
    IF (ir == MIN(n, i+m)) THEN
      tstr = .true.
    ELSE
      CALL right(n, h, ir, irr)
      tstr = ctest(a, il, ir, irr)
    END IF
    h(il) = tstl
    h(ir) = tstr
    IF (tstl.AND.tstr) RETURN
    IF(.NOT.tstl) il = ill
    IF(.NOT.tstr) ir = irr
  END DO
END IF

RETURN
END SUBROUTINE cmerge



!************************************************************************
!                             FUNCTION CTEST                            *
!************************************************************************
! Test the convexity of the angle formed by (IL,A(IL)), (I,A(I)),       *
! (IR,A(IR)) at the vertex (I,A(I)), up to within the tolerance         *
! TOLER. If convexity holds then the function is set to .TRUE.,         *
! otherwise CTEST=.FALSE. The parameter TOLER is set to 0.4 by default. *
!************************************************************************
! Input variables:                                                      *
!     A       : vector of double                                        *
!     IL,I,IR : integers such that IL < I < IR                          *
!************************************************************************
! Output:                                                               *
!     .TRUE. if the angle formed by (IL,A(IL)), (I,A(I)), (IR,A(IR)) at *
!            the vertex (I,A(I)), is convex up to within the tolerance  *
!            TOLER, i.e., if                                            *
!            (A(I)-A(IL))*(IR-I)-(A(IR)-A(I))*(I-IL)>TOLER.             *
!     .FALSE.,  otherwise.                                              *
!************************************************************************
FUNCTION ctest(a, il, i, ir) RESULT(OK)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: i, il, ir
REAL (KIND=dp), INTENT(IN) :: a(:)
LOGICAL                    :: OK

! Local variables
REAL (KIND=dp)            :: s1, s2
REAL (KIND=dp), PARAMETER :: toler = 0.4D0

s1 = a(i) - a(il)
s2 = a(ir) - a(i)
s1 = s1*(ir-i)
s2 = s2*(i-il)
OK = .false.
IF(s1 > (s2+toler)) OK = .true.
RETURN
END FUNCTION ctest


END MODULE poly_zeroes
