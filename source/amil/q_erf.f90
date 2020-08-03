SUBROUTINE q_erf(x, erf, erfc)

!   Calculate the error function & its complement in quadruple precision.
!   0 <= x <= 1 uses series 7.1.5 from Abramowitz & Stegun
!   1 < x < 4 uses Chebyshev series
!   4 <= x uses a continued fraction
!
!   If x < 0 then erfc(x) = 1 + erf(|x|) and erf(x) = - erf(|x|)
!
!   WARNING: x must NOT occupy the same locations as any of
!            the other arguments if x < 0.
!
!   Programmer: Alan Miller    Alan.Miller @ vic.cmis.csiro.au
!               www.ozemail.com.au/~milleraj   FAX: (+61) 3-9545-8080
!   Latest revision - (Fortran 77 version) 11 May 1988
!   Latest revision - 4 February 1998

USE quadruple_precision
IMPLICIT NONE
TYPE (quad), INTENT(IN)            :: x
TYPE (quad), INTENT(OUT), OPTIONAL :: erf, erfc

!     Chebyshev coefficients
TYPE (quad), PARAMETER :: coeff(40) = (/ quad(0.4888515056984716D+00, -.1253349490474187D-15),  &
         quad(-.1358986289159460D+00, 0.1828677302762691D-16),  &
         quad(0.3563511414853497D-01, -.1026978994630672D-17),  &
         quad(-.8884098618222492D-02, 0.1869288268188562D-17),  &
         quad(0.2118423997553177D-02, 0.1709987556862464D-18),  &
         quad(-.4853987640800370D-03, -.6960364726525852D-19),  &
         quad(0.1072737545377177D-03, 0.1534550435363124D-19),  &
         quad(-.2293663448190487D-04, 0.6928883131080608D-20),  &
         quad(0.4756874691575556D-05, -.1627977086918653D-20),  &
         quad(-.9589944832662770D-06, -.2249363452839194D-22),  &
         quad(0.1882900688295897D-06, -.1968876883111176D-22),  &
         quad(-.3606320662143092D-07, 0.7483138202418590D-23),  &
         quad(0.6747602051937456D-08, -.1005646109123554D-23),  &
         quad(-.1234907042333986D-08, -.8553919151789978D-25),  &
         quad(0.2213150119337267D-09, -.8378636504284192D-27),  &
         quad(-.3887955444797815D-10, 0.4599438293044550D-27),  &
         quad(0.6701388253274512D-11, -.4533239266031100D-27),  &
         quad(-.1134237426103359D-11, -.2592968654684961D-27),  &
         quad(0.1886555920176928D-12, -.9535048043067878D-29),  &
         quad(-.3085780489354876D-13, 0.7801557018723816D-29),  &
         quad(0.4966710292058390D-14, -.1626370800837245D-29),  &
         quad(-.7871154824211008D-15, 0.8389350837750800D-31),  &
         quad(0.1228887884514987D-15, -.5265492467954702D-31),  &
         quad(-.1891087633190591D-16, -.2388153131040173D-32),  &
         quad(0.2869742042877464D-17, 0.7703719777548944D-34),  &
         quad(-.4296338289538412D-18, 0.9629649721936182D-35),  &
         quad(0.6348319360536490D-19, -.7222237291452136D-35),  &
         quad(-.9261746469308854D-20, -.1805559322863034D-35),  &
         quad(0.1334626827242897D-20, 0.2633107345841924D-36),  &
         quad(-.1900244302571750D-21, 0.3761581922631320D-37),  &
         quad(0.2674133030089899D-22, -.4701977403289150D-38),  &
         quad(-.3720610299437181D-23, -.2938735877055719D-39),  &
         quad(0.5119522849758822D-24, -.1578156674465586D-56),  &
         quad(-.6968654298053806D-25, -.2295887403949780D-41),  &
         quad(0.9386117844916958D-26, -.2008901478456058D-41),  &
         quad(-.1251273661158375D-26, 0.5022253696140144D-42),  &
         quad(0.1651286941713772D-27, -.3587324068671532D-43),  &
         quad(-.2157623168556909D-28, -.5605193857299268D-45),  &
         quad(0.2771990968956549D-29, -.3503246160812043D-45),  &
         quad(-.3468792422835851D-30, 0.7006492321624086D-46) /)

!       Local variables

TYPE (quad) :: d, dd, term, arg, y2, sv, xx, y, temp, erff, erfcc
REAL (dp)   :: a
INTEGER     :: m, j
LOGICAL     :: small
TYPE (quad), PARAMETER :: qone = quad(1._dp, 0._dp), qzero = quad(0._dp, 0._dp)

IF (x%hi < 0.d0) THEN
  arg = -x
ELSE
  arg = x
END IF

IF (arg%hi > 1._dp) GO TO 20
!----------------------------------------------------------------------
!
!                  2x          x^2   x^4    x^6
!       erf(x) = --------.(1 - --- + ---- - ---- + .. )
!                sqrt(pi)      1.3   2!.5   3!.7
!
!
small = .false.
term = two_on_rtpi * arg
erff = term
xx = arg * arg
m = 1
DO
  IF (small) THEN
    term%hi = - term%hi * xx%hi / m
    temp%hi = term%hi / (2*m + 1)
    temp%lo = 0._dp
  ELSE
    term = - term * xx
    SELECT CASE (m)
      CASE (1)
        term = term
      CASE (2)
        term = SCALE(term, -1)
      CASE (4)
        term = SCALE(term, -2)
      CASE (8)
        term = SCALE(term, -3)
      CASE DEFAULT
        term = term / DBLE(m)
    END SELECT
    temp = term / DBLE(2*m+1)
  END IF
  erff = erff + temp
  m = m + 1
  IF (.NOT. small .AND. ABS(temp%hi) < erff%hi * 1.d-16) small = .true.
  IF (ABS(temp%hi) < erff%hi * 1.d-30) EXIT
END DO

erfcc = qone - erff

GO TO 60
!
!----------------------------------------------------------------------
!
!    Use Chebyshev series for 1 < x < 4.
!    Adapted from the routine CHEBEV from 'Numerical Recipes' by Press,
!    Flannery, Teukolsky & Vetterling, Cambridge Uni Press, 1986.
!
20 IF (arg%hi >= 4._dp) GO TO 40
d = qzero
dd = qzero
y = (arg - 2.5_dp) / 1.5_dp
y2 = SCALE(y, 1)
DO j = 40, 2, -1
  sv = d
  d = y2 * d - dd + coeff(j)
  dd = sv
END DO
term = y * d - dd + SCALE(coeff(1), -1)

!     Finally multiply by exp(-x^2)

erfcc = term * EXP(-x*x)
erff = qone - erfcc

GO TO 60

!----------------------------------------------------------------------
!
!       Use a continued fraction for x >= 4.0

40 a = 0.5D0 * INT(230._dp/arg%hi) + 2._dp
erfcc = qzero
DO
  erfcc = erfcc + arg
  erfcc = quad(a, 0._dp) / erfcc
  a = a - 0.5D0
  IF (a < 0.1D0) EXIT
END DO
erfcc = erfcc + arg
erfcc = qone / erfcc

!       Multiply result by exp(-x**2) / sqrt(pi)

erfcc = erfcc * EXP(-x*x) / sqrtpi
erff = qone - erfcc

!       Replace erf & erfc if argument was negative.

60 IF (x%hi < 0._dp) THEN
  erfcc = qone + erff
  erff = -erff
END IF

IF (PRESENT(erf)) erf = erff
IF (PRESENT(erfc)) erfc = erfcc

RETURN
END SUBROUTINE q_erf



PROGRAM test_q_erf
! Test quadruple-precision erf function by approximating its first derivative
! using divided differences.
! f'(x) = [f(x-2h) - 8.f(x-h) + 8.f(x+h) - f(x-2h)] / (12.h)  approx.

USE quadruple_precision
IMPLICIT NONE

TYPE (quad) :: x, h, arg, erf_2, erf_1, erf1, erf2, erfc, av_erfc,  &
               d_est, deriv, error
TYPE (quad), PARAMETER :: qone = quad(1._dp, 0._dp)

INTERFACE
  SUBROUTINE q_erf(x, erf, erfc)
    USE quadruple_precision
    IMPLICIT NONE
    TYPE (quad), INTENT(IN)            :: x
    TYPE (quad), INTENT(OUT), OPTIONAL :: erf, erfc
  END SUBROUTINE q_erf
END INTERFACE

DO
  WRITE(*, '(a)', ADVANCE='NO') ' Enter x: '
  READ(*, *) x%hi
  x%lo = 0._dp
  h = SCALE(qone, -24)
  IF (ABS(x%hi) < 0.477_dp) THEN       ! Use erf if erf(x) > 0.5
    arg = x - 2._dp * h
    CALL q_erf(arg, erf_2, av_erfc)
    arg = x - h
    CALL q_erf(arg, erf_1, erfc)
    av_erfc = av_erfc + erfc
    arg = x + h
    CALL q_erf(arg, erf1, erfc)
    av_erfc = av_erfc + erfc
    arg = x + 2._dp * h
    CALL q_erf(arg, erf2, erfc)
    av_erfc = (av_erfc + erfc) / 4._dp
    d_est = (erf_2 - erf2 + SCALE( erf1 - erf_1, 3)) / (12._dp * h)

  ELSE                                 ! Else use erfc(x)
    arg = x - 2._dp * h
    CALL q_erf(arg, erfc=av_erfc)
    erf_2 = av_erfc
    arg = x - h
    CALL q_erf(arg, erfc=erf_1)
    av_erfc = av_erfc + erf_1
    arg = x + h
    CALL q_erf(arg, erfc=erf1)
    av_erfc = av_erfc + erf1
    arg = x + 2._dp * h
    CALL q_erf(arg, erfc=erf2)
    av_erfc = (av_erfc + erf2) / 4._dp
    d_est = (erf2 - erf_2 + SCALE( erf_1 - erf1, 3)) / (12._dp * h)
  END IF

  deriv = two_on_rtpi * EXP(-x*x)
  error = d_est - deriv
  WRITE(*, '(a)') '      Estimate        Derivative         Error          erfc'
  WRITE(*, '(2f20.16, g12.4, f20.16)') d_est%hi, deriv%hi, error%hi, av_erfc%hi
END DO

STOP
END PROGRAM test_q_erf
