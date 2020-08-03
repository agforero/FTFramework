MODULE Bivariate_Trivariate_Normal

! Appendix: Fortran Source Code
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-02-25  Time: 11:57:13
! Latest revision - 25 February 2001
! amiller @ bigpond.net.au
 
!     Code from Alan Genz's web site for bivariate & trivariate
!     normal integrals, adapted by Alan Miller.

!     http://www.sci.wsu.edu/math/faculty/genz/papers/

!     Alan C Genz
!     2001-01-22


IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /trvbkd/b2p, b3p, rho21, rho31, rho
REAL (dp), SAVE  :: b2p, b3p, rho21, rho31, rho


CONTAINS


FUNCTION bvn( sh, sk, r ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.
!     This function uses an algorithm given in the paper
!        "Numerical Computation of Bivariate and
!             Trivariate Normal Probabilities", by
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu

! BVN - calculate the probability that X is larger than SH, and Y is
!       larger than SK.

! Parameters

!   SH  REAL, integration limit
!   SK  REAL, integration limit
!   R   REAL, correlation coefficient
!   LG  INTEGER, number of Gauss Rule Points and Weights


REAL, INTENT(IN)  :: sh, sk, r
REAL              :: fn_val

REAL, PARAMETER     :: zero = 0, twopi = 6.283185
INTEGER, PARAMETER  :: lg = 5
REAL, PARAMETER     :: x(lg) =  &
                       (/ -0.9061798, -0.5384693, 0.0, 0.5384693, 0.9061798 /)
!          X(I) = +-SQRT((5+-SQRT(10/7))/9), 0
REAL, PARAMETER     :: w(lg) =  &
                       (/ 0.2369269,  0.4786287, 0.5688889, 0.4786287, 0.2369269 /)

INTEGER  :: i
REAL     :: as, a, b, c, rs, xs
REAL     :: sn, asr, h, k, bs, hs, hk

h = sh
k = sk
hk = h*k
fn_val = 0

!     Absolute value of the correlation coefficient is less than 0.8

IF ( ABS(r) < 0.8 ) THEN
  hs = ( h*h + k*k )/2
  asr = ASIN( r )/2
  DO i = 1, lg
    sn = SIN( asr*(x(i)+1) )
    fn_val = fn_val + w(i)*EXP( ( sn*hk - hs )/( 1 - sn*sn ) )
  END DO
  fn_val = asr*fn_val/twopi + phi(-h)*phi(-k)
ELSE
  
!     For larger correlation coefficient
  
  IF ( r < 0 ) THEN
    k = -k
    hk = -hk
  END IF
  IF ( ABS(r) < 1 ) THEN
    as = ( 1 - r )*( 1 + r )
    a = SQRT( as )
    bs = ( h - k )**2
    c = ( 4 - hk )/8
    fn_val = a*EXP( -(bs/as + hk)/2 )*( 1 - c*(bs - as)/3 )
    IF ( -hk < 100 ) THEN
      b = SQRT(bs)
      fn_val = fn_val - EXP(-hk/2)*SQRT(twopi)*phi(-b/a)*b*(1-c*bs/3)
    END IF
    a = a/2
    DO i = 1, lg
      xs = ( a*( x(i) + 1 ) )**2
      rs = SQRT( 1 - xs )
      asr = -( bs/xs + hk )/2
      IF ( asr > -100 ) fn_val = fn_val + a*w(i)*   &
                        ( EXP(-bs/(2*xs) - hk/(1+rs))/rs - EXP(asr)*(1+c*xs) )
    END DO
    fn_val = -fn_val/twopi
  END IF
  IF ( r > 0 ) fn_val =  fn_val + phi( -MAX( h, k ) )
  IF ( r < 0 ) fn_val = -fn_val + MAX( zero, phi(-h) - phi(-k) )
END IF

RETURN
END FUNCTION bvn



FUNCTION phi(x) RESULT(fn_val)

REAL, INTENT(IN)  :: x
REAL              :: fn_val

fn_val = phid( REAL(x, KIND=dp) )
RETURN
END FUNCTION phi



FUNCTION bvnd( dh, dk, r ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.

!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu

!    This function is based on the method described by
!        Drezner, Z and G.O. Wesolowsky, (1989),
!        On the computation of the bivariate normal inegral,
!        Journal of Statist. Comput. Simul. 35, pp. 101-107,
!    with major modifications for REAL (dp), and for |R| close to 1.

! BVND - calculate the probability that X is larger than DH, and Y is
!       larger than DK.

! Parameters

!   DH  REAL (dp), integration limit
!   DK  REAL (dp), integration limit
!   R   REAL (dp), correlation coefficient


REAL (dp), INTENT(IN)  :: dh, dk, r
REAL (dp)              :: fn_val

REAL (dp), PARAMETER  :: zero = 0.0_dp, twopi = 6.283185307179586_dp

REAL (dp), PARAMETER  :: x(10,3) = RESHAPE(  &
                             !     Gauss Legendre Points, N =  6
     (/ -0.9324695142031522_dp, -0.6612093864662647_dp,  &
        -0.2386191860831970_dp,  0.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
                             !     Gauss Legendre Points, N = 12
        -0.9815606342467191_dp, -0.9041172563704750_dp,  &
        -0.7699026741943050_dp, -0.5873179542866171_dp,  &
        -0.3678314989981802_dp, -0.1252334085114692_dp,  &
         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
                             !     Gauss Legendre Points, N = 20
        -0.9931285991850949_dp, -0.9639719272779138_dp,  &
        -0.9122344282513259_dp, -0.8391169718222188_dp,  &
        -0.7463319064601508_dp, -0.6360536807265150_dp,  &
        -0.5108670019508271_dp, -0.3737060887154196_dp,  &
        -0.2277858511416451_dp, -0.07652652113349733_dp /), (/ 10, 3 /) )
REAL (dp), PARAMETER  :: w(10,3) = RESHAPE(  &
                             !     Gauss Legendre Weights, N =  6
  (/ 0.1713244923791705_dp, 0.3607615730481384_dp, 0.4679139345726904_dp, &
     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
                          !     Gauss Legendre Weights, N = 12
     0.04717533638651177_dp, 0.1069393259953183_dp, 0.1600783285433464_dp, &
     0.2031674267230659_dp,  0.2334925365383547_dp, 0.2491470458134029_dp, &
     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
                          !     Gauss Legendre Weights, N = 20
     0.01761400713915212_dp, 0.04060142980038694_dp, 0.06267204833410906_dp, &
     0.08327674157670475_dp, 0.1019301198172404_dp, 0.1181945319615184_dp, &
     0.1316886384491766_dp,  0.1420961093183821_dp, 0.1491729864726037_dp, &
     0.1527533871307259_dp /), (/ 10, 3 /) )

INTEGER    :: i, is, lg, ng
REAL (dp)  :: as, a, b, c, d, rs, xs, bvn
REAL (dp)  :: sn, asr, h, k, bs, hs, hk

IF ( ABS(r) < 0.3_dp ) THEN
  ng = 1
  lg = 3
ELSE IF ( ABS(r) < 0.75_dp ) THEN
  ng = 2
  lg = 6
ELSE
  ng = 3
  lg = 10
END IF
h = dh
k = dk
hk = h*k
bvn = 0
IF ( ABS(r) < 0.925_dp ) THEN
  hs = ( h*h + k*k )/2
  asr = ASIN(r)
  DO i = 1, lg
    DO is = -1, 1, 2
      sn = SIN( asr*(  is*x(i,ng) + 1 )/2 )
      bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/( 1 - sn*sn ) )
    END DO
  END DO
  bvn = bvn*asr/( 2*twopi ) + phid(-h)*phid(-k)
ELSE
  IF ( r < 0 ) THEN
    k = -k
    hk = -hk
  END IF
  IF ( ABS(r) < 1 ) THEN
    as = ( 1 - r )*( 1 + r )
    a = SQRT(as)
    bs = ( h - k )**2
    c = ( 4 - hk )/8
    d = ( 12 - hk )/16
    asr = -( bs/as + hk )/2
    IF ( asr > -100 ) bvn = a*EXP(asr)  &
                       *( 1 - c*( bs - as )*( 1 - d*bs/5 )/3 + c*d*as*as/5 )
    IF ( -hk < 100 ) THEN
      b = SQRT(bs)
      bvn = bvn - EXP( -hk/2 )*SQRT(twopi)*phid(-b/a)*b  &
            *( 1 - c*bs*( 1 - d*bs/5 )/3 )
    END IF
    a = a/2
    DO i = 1, lg
      DO is = -1, 1, 2
        xs = ( a*(  is*x(i,ng) + 1 ) )**2
        rs = SQRT( 1 - xs )
        asr = -( bs/xs + hk )/2
        IF ( asr > -100 ) THEN
          bvn = bvn + a*w(i,ng)*EXP( asr )  &
                *( EXP( -hk*( 1 - rs )/( 2*( 1 + rs ) ) )/rs  &
                - ( 1 + c*xs*( 1 + d*xs ) ) )
        END IF
      END DO
    END DO
    bvn = -bvn/twopi
  END IF
  IF ( r > 0 ) bvn =  bvn + phid( -MAX( h, k ) )
  IF ( r < 0 ) bvn = -bvn + MAX( zero, phid(-h) - phid(-k) )
END IF

fn_val = bvn
RETURN
END FUNCTION bvnd



FUNCTION tvnd( limit, sigma ) RESULT(fn_val)

!     A function for computing trivariate normal probabilities.
!     This function uses an algorithm given in the paper
!        "Numerical Computation of Bivariate and
!             Trivariate Normal Probabilities",
!     by
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu

! TVND calculates the probability that X(I) < LIMIT(I), I = 1, 2, 3.

!     Bivariate normal distribution function BVND is required.

! Parameters

!   LIMIT  REAL (dp) array of three upper integration limits.
!   SIGMA  REAL (dp) array of three correlation coefficients,
!          SIGMA should contain the lower left portion of the
!          correlation matrix R.
!          SIGMA(1) = R(2,1), SIGMA(2) = R(3,1), SIGMA(3) = R(3,2).

!    TVND cuts the outer integral over -infinity to B1 to an integral
!      from -8.5 to B1 and then uses an adaptive integration method to
!      compute the integral of a bivariate normal distribution function.


REAL (dp), INTENT(IN)  :: limit(3), sigma(3)
REAL (dp)              :: fn_val

! EXTERNAL trvfnd
! COMMON /trvbkd/b2p, b3p, rho21, rho31, rho

LOGICAL   :: tail
REAL (dp) :: sq21, sq31
REAL (dp) :: b1, b2, b3, rho32
REAL (dp), PARAMETER :: sqtwpi = 2.506628274631000_dp, xcut = -8.5_dp,  &
                        eps = 5D-16

b1 = limit(1)
b2 = limit(2)
b3 = limit(3)
rho21 = sigma(1)
rho31 = sigma(2)
rho32 = sigma(3)
IF ( ABS(b2) >= MAX( ABS(b1), ABS(b3) ) ) THEN
  b1 = b2
  b2 = limit(1)
  rho31 = rho32
  rho32 = sigma(2)
ELSE IF ( ABS(b3) >= MAX( ABS(b1), ABS(b2) ) ) THEN
  b1 = b3
  b3 = limit(1)
  rho21 = rho32
  rho32 = sigma(1)
END IF
tail = .false.
IF ( b1 > 0 ) THEN
  tail = .true.
  b1 = -b1
  rho21 = -rho21
  rho31 = -rho31
END IF
IF ( b1 > xcut ) THEN
  IF ( 2*ABS(rho21) < 1 ) THEN
    sq21 = SQRT( 1 - rho21**2 )
  ELSE
    sq21 = SQRT( ( 1 - rho21 )*( 1 + rho21 ) )
  END IF
  IF ( 2*ABS(rho31) < 1 ) THEN
    sq31 = SQRT( 1 - rho31**2 )
  ELSE
    sq31 = SQRT( ( 1 - rho31 )*( 1 + rho31 ) )
  END IF
  rho = ( rho32 - rho21*rho31 )/(sq21*sq31)
  b2p = b2/sq21
  rho21 = rho21/sq21
  b3p = b3/sq31
  rho31 = rho31/sq31
  fn_val = adoned( trvfnd, xcut, b1, eps )/sqtwpi
ELSE
  fn_val = 0
END IF
IF ( tail ) fn_val = bvnd( -b2, -b3, rho32 ) - fn_val

RETURN
END FUNCTION tvnd



FUNCTION trvfnd(t) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: t
REAL (dp)              :: fn_val

! COMMON /trvbkd/b2, b3, rho21, rho31, rho

fn_val = EXP( -t*t/2 )*bvnd(t*rho21 - b2p, t*rho31 - b3p, rho)
RETURN
END FUNCTION trvfnd



FUNCTION adoned( f, a, b, tol ) RESULT(fin)

!     One Dimensional Globally Adaptive Integration Function

REAL (dp), INTENT(IN)  :: a, b, tol
REAL (dp)              :: fin

! EXTERNAL f

INTERFACE
  FUNCTION f(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

INTEGER, PARAMETER  :: nl = 100
REAL (dp)  :: ei(nl), ai(nl), bi(nl), fi(nl), ERR
INTEGER    :: i, im, ip

ip = 1
ai(ip) = a
bi(ip) = b
CALL krnrdd( ai(ip), bi(ip), f, ei(ip), fi(ip) )
im = 1
10   im = im + 1
bi(im) = bi(ip)
ai(im) = (ai(ip) + bi(ip))/2
bi(ip) = ai(im)
fin = fi(ip)
CALL krnrdd( ai(ip), bi(ip), f, ei(ip), fi(ip) )
CALL krnrdd( ai(im), bi(im), f, ei(im), fi(im) )
ERR = ABS( fin - fi(ip) - fi(im) )/2
ei(ip) = ei(ip) + ERR
ei(im) = ei(im) + ERR
ip = 1
ERR = 0
fin = 0
DO i = 1, im
  IF ( ei(i) > ei(ip) ) ip = i
  fin = fin + fi(i)
  ERR = ERR + ei(i)
END DO
IF ( ERR > tol .AND. im < nl ) GO TO 10

RETURN
END FUNCTION adoned



SUBROUTINE krnrdd( a, b, f, abserr, fn_val)

!     Kronrod Rule

REAL (dp), INTENT(IN)   :: a, b
REAL (dp), INTENT(OUT)  :: abserr, fn_val

! EXTERNAL f

INTERFACE
  FUNCTION f(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

REAL (dp)  :: abscis, center, fc, funsum, hflgth, resltg, resltk

!     The abscissae and weights are given for the interval (-1,1)
!     because of symmetry only the positive abscisae and their
!     corresponding weights are given.

!     XGK    - abscissae of the 2N+1-point Kronrod rule:
!              XGK(2), XGK(4), ...  N-point Gauss rule abscissae;
!              XGK(1), XGK(3), ...  abscissae optimally added
!              to the N-point Gauss rule.

!     WGK    - weights of the 2N+1-point Kronrod rule.

!     WG     - weights of the N-point Gauss rule.

INTEGER  :: j
INTEGER, PARAMETER :: n = 11

REAL (dp), PARAMETER  :: wg(0:n/2) = (/  &
     0.2729250867779007_dp, 0.05566856711617449_dp, 0.1255803694649048_dp,  &
     0.1862902109277352_dp, 0.2331937645919914_dp, 0.2628045445102478_dp  /)
REAL (dp), PARAMETER  :: xgk(0:n) = (/  &
     0.0000000000000000_dp, 0.9963696138895427_dp, 0.9782286581460570_dp,   &
     0.9416771085780681_dp, 0.8870625997680953_dp, 0.8160574566562211_dp,   &
     0.7301520055740492_dp, 0.6305995201619651_dp, 0.5190961292068118_dp,   &
     0.3979441409523776_dp, 0.2695431559523450_dp, 0.1361130007993617_dp /)
REAL (dp), PARAMETER  :: wgk(0:n) = (/  &
     0.1365777947111183_dp, 0.9765441045961290D-02, 0.2715655468210443D-01, &
     0.4582937856442671D-01, 0.6309742475037484D-01, 0.7866457193222764D-01, &
     0.9295309859690074D-01, 0.1058720744813894_dp, 0.1167395024610472_dp,  &
     0.1251587991003195_dp, 0.1312806842298057_dp, 0.1351935727998845_dp /)

!     List of major variables

!     CENTER  - mid point of the interval
!     HFLGTH  - half-length of the interval
!     ABSCIS   - abscissae
!     RESLTG   - result of the N-point Gauss formula
!     RESLTK   - result of the 2N+1-point Kronrod formula

hflgth = ( b - a )/2
center = ( b + a )/2

!           compute the 2N+1-point Kronrod approximation to
!           the integral, and estimate the absolute error.

fc = f(center)
resltg = fc*wg(0)
resltk = fc*wgk(0)
DO j = 1, n
  abscis = hflgth*xgk(j)
  funsum = f( center - abscis ) + f( center + abscis )
  resltk = resltk + wgk(j)*funsum
  IF( MOD( j, 2 ) == 0 ) resltg = resltg + wg(j/2)*funsum
END DO
fn_val = resltk*hflgth
abserr = 3*ABS( ( resltk - resltg )*hflgth )

RETURN
END SUBROUTINE krnrdd



FUNCTION phid(z) RESULT(fn_val)

! Normal distribution probabilities accurate to 1.e-15.
! Z = no. of standard deviations from the mean.
! P, Q = probabilities to the left & right of Z.   P + Q = 1.
!       PDF = the probability density.

!       Based upon algorithm 5666 for the error function, from:
!       Hart, J.F. et al, 'Computer Approximations', Wiley 1968

!       Programmer: Alan Miller, modified by Alan Genz

! Latest revision of Fortran 77 version - 30 March 1986
! Latest revision of Fortran 90 version - 23 June 1997

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN)  :: z
REAL (dp)              :: fn_val

! Local variables
REAL (dp), PARAMETER :: p0 = 220.2068679123761_dp, p1 = 221.2135961699311_dp,  &
                        p2 = 112.0792914978709_dp, p3 = 33.91286607838300_dp,  &
                        p4 = 6.373962203531650_dp, p5 = .7003830644436881_dp,  &
                        p6 = 0.03526249659989109_dp,  &
                        q0 = 440.4137358247522_dp, q1 = 793.8265125199484_dp,  &
                        q2 = 637.3336333788311_dp, q3 = 296.5642487796737_dp,  &
                        q4 = 86.78073220294608_dp, q5 = 16.06417757920695_dp,  &
                        q6 = 1.755667163182642_dp, q7 = .8838834764831844D-1, &
                        cutoff = 7.071067811865475_dp

REAL (dp)            :: zabs, expntl

zabs = ABS(z)

! |Z| > 37.

IF (zabs > 37._dp) THEN
    fn_val = 0._dp
  ELSE

! |Z| <= 37.

    expntl = EXP(-0.5_dp*zabs**2)

! |Z| < CUTOFF = 10/sqrt(2).

  IF (zabs < cutoff) THEN
      fn_val = expntl*((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs + p2)*zabs     &
                 + p1)*zabs + p0) / (((((((q7*zabs + q6)*zabs + q5)*zabs &
                 + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs +q0)
  
! |Z| >= CUTOFF.
  
  ELSE
    fn_val = expntl/(zabs + 1._dp/(zabs + 2._dp/(zabs + 3._dp/(zabs +  &
                     4._dp/(zabs + 0.65_dp)))))
  END IF
END IF

IF (z > 0.0_dp) fn_val = 1.0_dp - fn_val

RETURN
END FUNCTION phid

END MODULE Bivariate_Trivariate_Normal
