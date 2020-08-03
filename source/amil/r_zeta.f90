MODULE Riemann_zeta

! From CERNLIB, with permission.
! Code converted using TO_F90 by Alan Miller
! Date: 2002-02-27  Time: 16:11:48

! Author: K.S. Kolbig, 7 June 1992
! Revision 1.1.1.1  1996/04/01 15:02:00  mclareni
! Mathlib gen/C (C315)

! The Riemann Zeta function
! Reference:
! Cody, W.J., Hillstrom, K.E. & Thather, H.C., `Chebyshev approximations
! for the Riemann zeta function', Math. Comp., vol.25 (1971), 537-547.

! N.B. The value returned by the CERNLIB routine was (zeta - 1) when the
!      argument was > 1.  This version returns the value of zeta.


IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


FUNCTION r_zeta(x) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

REAL (dp), PARAMETER :: delta = 1.0D-13, z1 = 1.0_dp, hf = z1/2, th = z1/3,  &
                        pi = 3.14159265358979324_dp, pih = pi/2, pi2 = 2*pi

REAL (dp), PARAMETER :: p1(0:8) = (/ 1.28716812148244639D+10,  &
  1.37539693203702511D+10, 5.10665591836440610D+09, 8.56147100243331486D+08, &
  7.48361812438023298D+07, 4.86010658546188251D+06, 2.73957499022140609D+05, &
  4.63171084318342712D+03, 5.78758100409666066D+01 /)
REAL (dp), PARAMETER :: q1(0:8) = (/ 2.57433624296484624D+10,  &
  5.93816564867959016D+09, 9.00633037326123344D+08, 8.04253663428328989D+07, &
  5.60971175954192006D+06, 2.24743120289913752D+05, 7.57457890934153756D+03, &
 -2.37383578137377262D+01, 1.0_dp /)

REAL (dp), PARAMETER :: p2(0:8) = (/ -6.88197293216348954D+06,  &
  7.48218916305315972D+06, -2.07584504810211014D+06, 3.55302557096214295D+05, &
 -4.06706449551854889D+04, 3.19804864027146911D+03, -1.69820937033722853D+02, &
  5.61485842394289048D+00, -8.93888705926154944D-02 /)
REAL (dp), PARAMETER :: q2(0:8) = (/ -1.29725624934891554D+09,  &
 -9.48715407579907817D+08, -1.05496193474005203D+08, 4.67774488211993048D+06, &
  3.12936040573813534D+06, 4.59581803839305070D+05, 3.88176109610396834D+04,  &
  1.92561544834491423D+03, 5.12578125000000000D+01 /)

REAL (dp), PARAMETER :: p3(0:9) = (/ 1.66156480515774676D-11,  &
 -4.68068827660654529D-09, 5.83519727319147047D-07, -4.17644012643145602D-05, &
  1.85468422843597959D-03, -5.11288800220490241D-02, 8.10450231751100353D-01, &
 -5.69951948768478923D+00, 0.0_dp, 0.0_dp /)
REAL (dp), PARAMETER :: q3(0:9) = (/ -6.99562633519191655D-10,  &
 -1.77757961895149257D-08, -9.82231825734078036D-07, -2.84927282759096488D-05, &
 -5.81727909388048094D-04, -1.15848749169766586D-02, -1.28149124051978196D-01, &
 -1.11913057349097709D+00, -7.67928761604628813D-01,  1.0_dp /)

REAL (dp), PARAMETER :: p4(0:8) = (/ 1.03144877188859712D-15,  &
 -5.12584613964688241D-13, 1.12948794194873548D-10, -1.44234665373130952D-08, &
  1.16824676984458098D-06, -6.14975167990314806D-05, 2.05594677988830328D-03, &
 -3.99339429394668869D-02, 3.45234976736178457D-01 /)
REAL (dp), PARAMETER :: q4(0:8) = (/ 5.93959417288419050D-11,  &
 -6.04755359079991806D-09, 3.64680208668388563D-07, -1.29456905568011812D-05, &
  3.20189498470229250D-04, -5.07801557099994077D-03, 5.49628907881587266D-02, &
 -3.24517611155972419D-01, 1.0_dp /)

REAL (dp)  :: alfa, ap, aq, b0, b1, b2, f, h, t, v
INTEGER    :: ix, j

v = x
f = 1.0_dp
IF (x /= 0 .AND. x < hf) THEN
  ix = x - delta
  IF (ABS(ix-x) <= delta) THEN
    IF (MOD(-ix,2) == 0) THEN
      h = 0.0_dp
      GO TO 70
    ELSE
      v = 1.0_dp - x
      f = 2 * (-z1) ** ((1-ix)/2) * dgamma(v) / pi2 ** v
    END IF
  ELSE
    v = 1.0_dp - x
    f = 2 * SIN(pih*x) * dgamma(v) / pi2 ** v
  END IF
END IF

IF (x == 0) THEN
  h = -3 * hf
ELSE IF (x == 1) THEN
  fn_val = 0.0_dp
  WRITE(*, *) 'Riemanns ZETA(X) HAS POLE AT X = 1'
  RETURN
ELSE IF (v <= 5) THEN
  ap = p1(8)
  aq = q1(8)
  DO  j = 7, 0, -1
    ap = p1(j) + v * ap
    aq = q1(j) + v * aq
  END DO
  h = ap / (aq*(v-1)) - 1
ELSE IF (v <= 11) THEN
  t = th * (v-8)
  alfa = t + t
  b1 = 0
  b2 = 0
  DO  j = 8, 0, -1
    b0 = p2(j) + alfa * b1 - b2
    b2 = b1
    b1 = b0
  END DO
  h = b0 - t * b2
  b1 = 0
  b2 = 0
  DO  j = 8, 0, -1
    b0 = q2(j) + alfa * b1 - b2
    b2 = b1
    b1 = b0
  END DO
  h = h / (b0 - t*b2)
ELSE IF (v <= 25) THEN
  t = 1 / v
  ap = p3(7)
  DO  j = 6, 0, -1
    ap = p3(j) + t * ap
  END DO
  aq = q3(9)
  DO  j = 8, 0, -1
    aq = q3(j) + t * aq
  END DO
  h = hf ** (v - t*ap/aq)
ELSE IF (v <= 55) THEN
  t = 1 / v
  ap = p4(8)
  aq = q4(8)
  DO  j = 7, 0, -1
    ap = p4(j) + t * ap
    aq = q4(j) + t * aq
  END DO
  h = hf ** (v-t*ap/aq)
ELSE IF (v <= 90) THEN
  h = hf ** v + th ** v
ELSE
  h = hf ** v
END IF
IF (x < 1) h = f * (1 + h)

70 IF (x > 1.0_dp) THEN
   fn_val = 1.0_dp + h
ELSE
   fn_val = h
END IF

RETURN
END FUNCTION r_zeta



! $Log: gamma64.F,v $
! Revision 1.1.1.1  1996/04/01 15:01:54  mclareni
! Mathlib gen

FUNCTION dgamma(x) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

REAL (dp), PARAMETER  :: c(0:15) = (/ 3.65738772508338244_dp,  &
   1.95754345666126827_dp, 0.33829711382616039_dp, 0.04208951276557549_dp,  &
   0.00428765048212909_dp, 0.00036521216929462_dp, 0.00002740064222642_dp,  &
   0.00000181240233365_dp, 0.00000010965775866_dp, 0.00000000598718405_dp,  &
   0.00000000030769081_dp, 0.00000000001431793_dp, 0.00000000000065109_dp,  &
   0.00000000000002596_dp, 0.00000000000000111_dp, 0.00000000000000004_dp /)
REAL (dp)  :: alfa, b0, b1, b2, f, h, u
INTEGER    :: i

u = x
IF (u <= 0) THEN
  WRITE(*, '(a, g13.5, a)') ' Argument value: ', u, ' <= 0 for routine DGAMMA'
  h = 0
  RETURN
END IF

f = 1
IF (u < 3) THEN
  DO  i = 1, INT(4-u)
    f = f / u
    u = u + 1
  END DO
ELSE
  DO  i = 1, INT(u-3)
    u = u - 1
    f = f * u
  END DO
END IF
h = u + u - 7
alfa = h + h
b1 = 0
b2 = 0
DO  i = 15, 0, -1
  b0 = c(i) + alfa * b1 - b2
  b2 = b1
  b1 = b0
END DO

fn_val = f * (b0 - h*b2)
RETURN
END FUNCTION dgamma

END MODULE Riemann_zeta



PROGRAM Test_Riemann_zeta
USE Riemann_zeta
IMPLICIT NONE

! Output the LHS of Table 23.3 (page 811) of Abramowitz & Stegun

INTEGER             :: n
REAL (dp)           :: zeta
CHARACTER (LEN=21)  :: text

DO n = 2, 42
  zeta = r_zeta( DBLE(n) )
  WRITE(text, '(f18.16)') zeta
  text = text(:7) // ' ' // text(8:12) // ' ' // text(13:17) // ' ' // text(18:18)
  WRITE(*, '(i3, "  ", a)') n, text
END DO

STOP
END PROGRAM Test_Riemann_zeta
