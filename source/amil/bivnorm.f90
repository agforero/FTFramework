FUNCTION bvn ( lower, upper, infin, correl ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.
!     Extracted from Alan Genz's package for multivariate normal integration.

!  Parameters

!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, correlation coefficient.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN) :: lower(:), upper(:), correl
INTEGER, INTENT(IN)   :: infin(:)
REAL (dp)             :: fn_val

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

IF ( infin(1) == 2  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )    &
            - bvnu ( upper(1), lower(2), correl )  &
            - bvnu ( lower(1), upper(2), correl )  &
            + bvnu ( upper(1), upper(2), correl )
ELSE IF ( infin(1) == 2  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )  &
            - bvnu ( upper(1), lower(2), correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )  &
            - bvnu ( lower(1), upper(2), correl )
ELSE IF ( infin(1) == 2  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )  &
            - bvnu ( -lower(1), -upper(2), correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )  &
            - bvnu ( -upper(1), -lower(2), correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( lower(1), -upper(2), -correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( -upper(1), lower(2), -correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )
END IF

RETURN


CONTAINS


FUNCTION bvnu( sh, sk, r ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.

!       Yihong Ge
!       Department of Computer Science and Electrical Engineering
!       Washington State University
!       Pullman, WA 99164-2752
!       Email : yge@eecs.wsu.edu
!     and
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu

! BVN - calculate the probability that X is larger than SH and Y is
!       larger than SK.

! Parameters

!   SH  REAL, integration limit
!   SK  REAL, integration limit
!   R   REAL, correlation coefficient
!   LG  INTEGER, number of Gauss Rule Points and Weights

REAL (dp), INTENT(IN) :: sh, sk, r
REAL (dp)             :: fn_val

! Local variables
INTEGER              :: i, lg, ng
REAL (dp), PARAMETER :: twopi = 6.283185307179586
REAL (dp)            :: as, a, b, c, d, rs, xs
REAL (dp)            :: bvn, sn, asr, h, k, bs, hs, hk
!     Gauss Legendre Points and Weights, N =  6
! DATA ( w(i,1), x(i,1), i = 1,3) /  &
! 0.1713244923791705D+00,-0.9324695142031522D+00,  &
! 0.3607615730481384D+00,-0.6612093864662647D+00,  &
! 0.4679139345726904D+00,-0.2386191860831970D+00/
!     Gauss Legendre Points and Weights, N = 12
! DATA ( w(i,2), x(i,2), i = 1,6) /  &
! 0.4717533638651177D-01,-0.9815606342467191D+00,  &
! 0.1069393259953183D+00,-0.9041172563704750D+00,  &
! 0.1600783285433464D+00,-0.7699026741943050D+00,  &
! 0.2031674267230659D+00,-0.5873179542866171D+00,  &
! 0.2334925365383547D+00,-0.3678314989981802D+00,  &
! 0.2491470458134029D+00,-0.1252334085114692D+00/
!     Gauss Legendre Points and Weights, N = 20
! DATA ( w(i,3), x(i,3), i = 1,10) /  &
! 0.1761400713915212D-01,-0.9931285991850949D+00,  &
! 0.4060142980038694D-01,-0.9639719272779138D+00,  &
! 0.6267204833410906D-01,-0.9122344282513259D+00,  &
! 0.8327674157670475D-01,-0.8391169718222188D+00,  &
! 0.1019301198172404D+00,-0.7463319064601508D+00,  &
! 0.1181945319615184D+00,-0.6360536807265150D+00,  &
! 0.1316886384491766D+00,-0.5108670019508271D+00,  &
! 0.1420961093183821D+00,-0.3737060887154196D+00,  &
! 0.1491729864726037D+00,-0.2277858511416451D+00,  &
! 0.1527533871307259D+00,-0.7652652113349733D-01/
REAL (dp), PARAMETER :: w(10,3) = RESHAPE( (/      &
      0.1713244923791705D+00, 0.3607615730481384D+00, 0.4679139345726904D+00, &
        0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,                  &
      0.4717533638651177D-01, 0.1069393259953183D+00, 0.1600783285433464D+00, &
      0.2031674267230659D+00, 0.2334925365383547D+00, 0.2491470458134029D+00, &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      0.1761400713915212D-01, 0.4060142980038694D-01, 0.6267204833410906D-01, &
      0.8327674157670475D-01, 0.1019301198172404D+00, 0.1181945319615184D+00, &
      0.1316886384491766D+00, 0.1420961093183821D+00, 0.1491729864726037D+00, &
      0.1527533871307259D+00 /), (/ 10, 3 /) )
REAL (dp), PARAMETER :: x(10,3) = RESHAPE( (/      &
      -0.9324695142031522D+00, -0.6612093864662647D+00,   &
      -0.2386191860831970D+00,  0.0D0, 0.0D0, 0.0D0,      &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      -0.9815606342467191D+00, -0.9041172563704750D+00,   &
      -0.7699026741943050D+00, -0.5873179542866171D+00,   &
      -0.3678314989981802D+00, -0.1252334085114692D+00,   &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      -0.9931285991850949D+00, -0.9639719272779138D+00,   &
      -0.9122344282513259D+00, -0.8391169718222188D+00,   &
      -0.7463319064601508D+00, -0.6360536807265150D+00,   &
      -0.5108670019508271D+00, -0.3737060887154196D+00,   &
      -0.2277858511416451D+00, -0.7652652113349733D-01 /), (/ 10, 3 /) )

IF ( ABS(r) < 0.3 ) THEN
  ng = 1
  lg = 3
ELSE IF ( ABS(r) < 0.75 ) THEN
  ng = 2
  lg = 6
ELSE
  ng = 3
  lg = 10
END IF
h = sh
k = sk
hk = h*k
bvn = zero
IF ( ABS(r) < 0.925 ) THEN
  hs = ( h*h + k*k )/2
  asr = ASIN(r)
  DO  i = 1, lg
    sn = SIN(asr*( x(i,ng)+1 )/2)
    bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/(one - sn*sn ) )
    sn = SIN(asr*(-x(i,ng)+1 )/2)
    bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/(one - sn*sn ) )
  END DO
  bvn = bvn*asr/(2*twopi) + phi(-h)*phi(-k)
ELSE
  IF ( r < zero ) THEN
    k = -k
    hk = -hk
  END IF
  IF ( ABS(r) < one ) THEN
    as = ( one - r )*( one + r )
    a = SQRT(as)
    bs = ( h - k )**2
    c = ( 4. - hk )/8
    d = ( 12. - hk )/16.
    bvn = a*EXP( -(bs/as + hk)/2. )  &
    *( one - c*(bs - as)*(one - d*bs/5.)/3. + c*d*as*as/5. )
    IF ( hk > -160. ) THEN
      b = SQRT(bs)
      bvn = bvn - EXP(-hk/2)*SQRT(twopi)*phi(-b/a)*b  &
                      *( one - c*bs*( one - d*bs/5. )/3. )
    END IF
    a = a/2
    DO i = 1, lg
      xs = ( a*(x(i,ng) + one) )**2
      rs = SQRT( one - xs )
      bvn = bvn + a*w(i,ng)*( EXP( -bs/(2*xs) - hk/(1+rs) )/rs  &
                - EXP( -(bs/xs+hk)/2. )*( one + c*xs*( one + d*xs ) ) )
      xs = as*(-x(i,ng) + one)**2/4.
      rs = SQRT( 1 - xs )
      bvn = bvn + a*w(i,ng)*EXP( -(bs/xs + hk)/2 )  &
                * ( EXP( -hk*(one - rs)/(2*(one + rs)) )/rs - &
                       ( one + c*xs*( one + d*xs ) ) )
    END DO
    bvn = -bvn/twopi
  END IF
  IF ( r > 0 ) bvn =  bvn + phi( -MAX( h, k ) )
  IF ( r < 0 ) bvn = -bvn + MAX( zero, phi(-h) - phi(-k) )
END IF
fn_val = bvn

RETURN
END FUNCTION bvnu



SUBROUTINE normp(z, p, q, pdf)

! Normal distribution probabilities accurate to 1.e-15.
! Z = no. of standard deviations from the mean.
! P, Q = probabilities to the left & right of Z.   P + Q = 1.
! PDF = the probability density.

! Based upon algorithm 5666 for the error function, from:
! Hart, J.F. et al, 'Computer Approximations', Wiley 1968

! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 30 March 1986
! Latest revision of Fortran 90 version - 12 August 1997

IMPLICIT NONE
REAL (dp), INTENT(IN)            :: z
REAL (dp), INTENT(OUT), OPTIONAL :: p, q, pdf

! Local variables
REAL (dp), PARAMETER :: p0 = 220.2068679123761D0, p1 = 221.2135961699311D0,  &
                        p2 = 112.0792914978709D0, p3 = 33.91286607838300D0,  &
                        p4 = 6.373962203531650D0, p5 = .7003830644436881D0,  &
                        p6 = .3526249659989109D-01,  &
                        q0 = 440.4137358247522D0, q1 = 793.8265125199484D0,  &
                        q2 = 637.3336333788311D0, q3 = 296.5642487796737D0,  &
                        q4 = 86.78073220294608D0, q5 = 16.06417757920695D0,  &
                        q6 = 1.755667163182642D0, q7 = .8838834764831844D-1, &
                        cutoff = 7.071D0, root2pi = 2.506628274631001D0
REAL (dp)            :: zabs, expntl, pp, qq, ppdf

zabs = ABS(z)

! |Z| > 37.

IF (zabs > 37.d0) THEN
  IF (PRESENT(pdf)) pdf = zero
  IF (z > zero) THEN
    IF (PRESENT(p)) p = one
    IF (PRESENT(q)) q = zero
  ELSE
    IF (PRESENT(p)) p = zero
    IF (PRESENT(q)) q = one
  END IF
  RETURN
END IF

! |Z| <= 37.

expntl = EXP(-0.5D0*zabs**2)
ppdf = expntl/root2pi
IF (PRESENT(pdf)) pdf = ppdf

! |Z| < CUTOFF = 10/sqrt(2).

IF (zabs < cutoff) THEN
  pp = expntl*((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs + p2)*zabs     &
                   + p1)*zabs + p0) / (((((((q7*zabs + q6)*zabs + q5)*zabs &
                   + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs +q0)

! |Z| >= CUTOFF.

ELSE
  pp = ppdf/(zabs + one/(zabs + 2.d0/(zabs + 3.d0/(zabs + 4.d0/(zabs + 0.65D0)))))
END IF

IF (z < zero) THEN
  qq = one - pp
ELSE
  qq = pp
  pp = one - qq
END IF

IF (PRESENT(p)) p = pp
IF (PRESENT(q)) q = qq

RETURN
END SUBROUTINE normp



FUNCTION phi(z) RESULT(p)
REAL (dp), INTENT(IN) :: z
REAL (dp)             :: p

CALL normp(z, p)

RETURN
END FUNCTION phi

END FUNCTION bvn

