RECURSIVE FUNCTION prtrng(q, v, r) RESULT(prob)

! N.B. Argument IFAULT has been removed.

! Code converted using TO_F90 by Alan Miller
! Date: 1999-04-02  Time: 20:07:46

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

REAL (dp), INTENT(IN)  :: q
REAL (dp), INTENT(IN)  :: v
REAL (dp), INTENT(IN)  :: r
REAL (dp)              :: prob

!   Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2

!   Evaluates the probability from 0 to q for a studentized range
!   having v degrees of freedom and r samples.

!   Uses subroutine ALNORM = algorithm AS66.

!   Arrays vw and qw store transient values used in the quadrature
!   summation.  Node spacing is controlled by step.  pcutj and pcutk control
!   truncation.  Minimum and maximum number of steps are controlled by
!   jmin, jmax, kmin and kmax.  Accuracy can be increased by use of a finer
!   grid - Increase sizes of arrays vw and qw, and jmin, jmax, kmin, kmax and
!   1/step proportionally.

REAL (dp) :: vw(30), qw(30)
REAL (dp) :: g, gmid, r1, c, h, v2, gstep, pk1, pk2, gk, pk
REAL (dp) :: w0, pz, x, hj, ehj, pj
INTEGER   :: ifault, j, jj, jump, k
REAL (dp), PARAMETER :: pcutj = 0.00003_dp, pcutk = 0.0001_dp, step = 0.45_dp, &
                        vmax = 120.0_dp, zero = 0.0_dp, fifth = 0.2_dp,  &
                        half = 0.5_dp, one = 1.0_dp, two = 2.0_dp,   &
                        cv1 = 0.193064705_dp, cv2 = 0.293525326_dp,  &
                        cvmax = 0.39894228_dp, cv(4) = (/ 0.318309886_dp,  &
                        -0.268132716D-2, 0.347222222D-2, 0.833333333D-1 /)
INTEGER, PARAMETER   :: jmin = 3, jmax = 15, kmin = 7, kmax = 15

!        Check initial values

prob = zero
ifault = 0
IF (v < one .OR. r < two) ifault = 1
IF (q <= zero .OR. ifault == 1) GO TO 99

!        Computing constants, local midpoint, adjusting steps.

g = step * r ** (-fifth)
gmid = half * LOG(r)
r1 = r - one
c = LOG(r * g * cvmax)
IF(v > vmax) GO TO 20

h = step * v ** (-half)
v2 = v * half
IF (v == one) c = cv1
IF (v == two) c = cv2
IF (.NOT. (v == one .OR. v == two)) c = SQRT(v2)  &
    * cv(1) / (one + ((cv(2) / v2 + cv(3)) / v2 + cv(4)) / v2)
c = LOG(c * r * g * h)

!        Computing integral
!        Given a row k, the procedure starts at the midpoint and works
!        outward (index j) in calculating the probability at nodes
!        symmetric about the midpoint.  The rows (index k) are also
!        processed outwards symmetrically about the midpoint.  The
!        centre row is unpaired.

20 gstep = g
qw(1) = -one
qw(jmax + 1) = -one
pk1 = one
pk2 = one
DO  k = 1, kmax
  gstep = gstep - g
  21 gstep = -gstep
  gk = gmid + gstep
  pk = zero
  IF (pk2 <= pcutk .AND. k > kmin) GO TO 26
  w0 = c - gk * gk * half
  pz = alnorm(gk, .true.)
  x = alnorm(gk - q, .true.) - pz
  IF (x > zero) pk = EXP(w0 + r1 * LOG(x))
  IF (v > vmax) GO TO 26
  
  jump = -jmax
  22 jump = jump + jmax
  DO  j = 1, jmax
    jj = j + jump
    IF (qw(jj) > zero) GO TO 23
    hj = h * j
    IF (j < jmax) qw(jj + 1) = -one
    ehj = EXP(hj)
    qw(jj) = q * ehj
    vw(jj) = v * (hj + half - ehj * ehj * half)
    
    23 pj = zero
    x = alnorm(gk - qw(jj), .true.) - pz
    IF (x > zero) pj = EXP(w0 + vw(jj) + r1 * LOG(x))
    pk = pk + pj
    IF (pj > pcutj) CYCLE
    IF (jj > jmin .OR. k > kmin) EXIT
  END DO
  h = -h
  IF (h < zero) GO TO 22
  
  26 prob = prob + pk
  IF (k > kmin .AND. pk <= pcutk .AND. pk1 <= pcutk) GO TO 99
  pk2 = pk1
  pk1 = pk
  IF (gstep > zero) GO TO 21
END DO

99 IF (ifault /= 0) WRITE(*, '(a, i3)') ' IFAULT = ', ifault
RETURN

CONTAINS


FUNCTION qtrng(p, v, r) RESULT(quantile)

! N.B. Argument IFAULT has been removed.

REAL (dp), INTENT(IN)  :: p
REAL (dp), INTENT(IN)  :: v
REAL (dp), INTENT(IN)  :: r
REAL (dp)              :: quantile

!   Algorithm AS 190.1  Appl. Statist. (1983) Vol.32, No. 2

!   Approximates the quantile p for a studentized range distribution
!   having v degrees of freedom and r samples for probability 0.9 < p < 0.99.

!   Uses functions  alnorm, ppnd, prtrng and qtrng0 -
!   Algorithms AS 66, AS 241, AS 190 and AS 190.2

REAL (dp) :: q1, p1, q2, p2, d, e1, e2
INTEGER   :: ifault, j, nfault
INTEGER, PARAMETER    :: jmax = 8
REAL (dp), PARAMETER  :: pcut = 0.001_dp, p75 = 0.75_dp, p80 = 0.80_dp,  &
                         p90 = 0.9_dp, p99 = 0.99_dp, p995 = 0.995_dp,   &
                         p175 = 1.75_dp, one = 1.0_dp, two = 2.0_dp,    &
                         five = 5.0_dp, eps = 1.0D-04

!        Check input parameters

ifault = 0
nfault = 0
IF (v < one .OR. r < two) ifault = 1
IF (p < p90 .OR. p > p99) ifault = 2
IF (ifault /= 0) GO TO 99

!        Obtain initial values

q1 = qtrng0(p, v, r)
p1 = prtrng(q1, v, r)
IF (nfault /= 0) GO TO 99
quantile = q1
IF (ABS(p1-p) < pcut) GO TO 99
IF (p1 > p) p1 = p175 * p - p75 * p1
IF (p1 < p) p2 = p + (p - p1) * (one - p) / (one - p1) * p75
IF (p2 < p80) p2 = p80
IF (p2 > p995) p2 = p995
q2 = qtrng0(p2, v, r)
IF (nfault /= 0) GO TO 99

!        Refine approximation

DO  j = 2, jmax
  p2 = prtrng(q2, v, r)
  IF (nfault /= 0) GO TO 99
  e1 = p1 - p
  e2 = p2 - p
  quantile = (q1 + q2) / two
  d = e2 - e1
  IF (ABS(d) > eps) quantile = (e2 * q1 - e1 * q2) / d
  IF(ABS(e1) < ABS(e2)) GO TO 12
  q1 = q2
  p1 = p2
  12   IF (ABS(p1 - p) < pcut * five) GO TO 99
  q2 = quantile
END DO

99 IF (nfault /= 0) ifault = 9
IF (ifault /= 0) WRITE(*, '(a, i4, a)') ' IFAULT = ', ifault, ' from QTRNG'
RETURN
END FUNCTION qtrng



FUNCTION qtrng0(p, v, r) RESULT(initq)

! N.B. Argument IFAULT has been removed.

REAL (dp), INTENT(IN)  :: p
REAL (dp), INTENT(IN)  :: v
REAL (dp), INTENT(IN)  :: r
REAL (dp)              :: initq

!  Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2

!  Calculates an initial quantile p for a studentized range distribution
!  having v degrees of freedom and r samples for probability p, 0.8 < p < 0.995

!  Uses function ppnd - Algorithm AS 241

INTEGER    :: ifault
REAL (dp)  :: q, t
REAL (dp), PARAMETER  :: vmax = 120.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, c1 = 0.8843_dp, c2 = 0.2368_dp, &
                         c3 = 1.214_dp, c4 = 1.208_dp, c5 = 1.4142_dp

CALL ppnd16(half + half * p, t, ifault)
IF (v < vmax) t = t + (t * t* t + t) / v / four
q = c1 - c2 * t
IF (v < vmax) q = q - c3 / v + c4 * t / v
initq = t * (q * LOG(r - one) + c5)

RETURN
END FUNCTION qtrng0



FUNCTION alnorm(x, upper) RESULT(fn_val)

!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

REAL (dp), INTENT(IN) :: x
LOGICAL, INTENT(IN)   :: upper
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp, &
                        con = 1.28_dp
REAL (dp) :: z, y
LOGICAL   :: up

!*** machine dependent constants
REAL (dp), PARAMETER :: ltone = 7.0_dp, utzero = 18.66_dp

REAL (dp), PARAMETER :: p = 0.398942280444_dp, q = 0.39990348504_dp,   &
                        r = 0.398942280385_dp, a1 = 5.75885480458_dp,  &
                        a2 = 2.62433121679_dp, a3 = 5.92885724438_dp,  &
                        b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
                        c1 = -3.8052E-8_dp, c2 = 3.98064794E-4_dp,     &
                        c3 = -0.151679116635_dp, c4 = 4.8385912808_dp, &
                        c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
                        d1 = 1.00000615302_dp, d2 = 1.98615381364_dp,  &
                        d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, &
                        d5 = 30.789933034_dp

up = upper
z = x
IF(z >=  zero) GO TO 10
up = .NOT. up
z = -z
10 IF(z <= ltone .OR. up .AND. z <= utzero) GO TO 20
fn_val = zero
GO TO 40
20 y = half*z*z
IF(z > con) GO TO 30

fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
GO TO 40
30 fn_val = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
40 IF(.NOT. up) fn_val = one - fn_val

RETURN
END FUNCTION alnorm



SUBROUTINE ppnd16 (p, normal_dev, ifault)

! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.

! The hash sums below are the sums of the mantissas of the
! coefficients.   They are included for use in checking
! transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

REAL (dp), INTENT(IN)   :: p
INTEGER, INTENT(OUT)    :: ifault
REAL (dp), INTENT(OUT)  :: normal_dev

! Local variables

REAL (dp) :: zero = 0.d0, one = 1.d0, half = 0.5d0, split1 = 0.425d0,  &
             split2 = 5.d0, const1 = 0.180625d0, const2 = 1.6d0, q, r

! Coefficients for P close to 0.5

REAL (dp) :: a0 = 3.3871328727963666080D0, &
             a1 = 1.3314166789178437745D+2, &
             a2 = 1.9715909503065514427D+3, &
             a3 = 1.3731693765509461125D+4, &
             a4 = 4.5921953931549871457D+4, &
             a5 = 6.7265770927008700853D+4, &
             a6 = 3.3430575583588128105D+4, &
             a7 = 2.5090809287301226727D+3, &
             b1 = 4.2313330701600911252D+1, &
             b2 = 6.8718700749205790830D+2, &
             b3 = 5.3941960214247511077D+3, &
             b4 = 2.1213794301586595867D+4, &
             b5 = 3.9307895800092710610D+4, &
             b6 = 2.8729085735721942674D+4, &
             b7 = 5.2264952788528545610D+3
! HASH SUM AB    55.8831928806149014439

! Coefficients for P not close to 0, 0.5 or 1.

REAL (dp) :: c0 = 1.42343711074968357734D0, &
             c1 = 4.63033784615654529590D0, &
             c2 = 5.76949722146069140550D0, &
             c3 = 3.64784832476320460504D0, &
             c4 = 1.27045825245236838258D0, &
             c5 = 2.41780725177450611770D-1, &
             c6 = 2.27238449892691845833D-2, &
             c7 = 7.74545014278341407640D-4, &
             d1 = 2.05319162663775882187D0, &
             d2 = 1.67638483018380384940D0, &
             d3 = 6.89767334985100004550D-1, &
             d4 = 1.48103976427480074590D-1, &
             d5 = 1.51986665636164571966D-2, &
             d6 = 5.47593808499534494600D-4, &
             d7 = 1.05075007164441684324D-9
! HASH SUM CD    49.33206503301610289036

! Coefficients for P near 0 or 1.

REAL (dp) :: e0 = 6.65790464350110377720D0, &
             e1 = 5.46378491116411436990D0, &
             e2 = 1.78482653991729133580D0, &
             e3 = 2.96560571828504891230D-1, &
             e4 = 2.65321895265761230930D-2, &
             e5 = 1.24266094738807843860D-3, &
             e6 = 2.71155556874348757815D-5, &
             e7 = 2.01033439929228813265D-7, &
             f1 = 5.99832206555887937690D-1, &
             f2 = 1.36929880922735805310D-1, &
             f3 = 1.48753612908506148525D-2, &
             f4 = 7.86869131145613259100D-4, &
             f5 = 1.84631831751005468180D-5, &
             f6 = 1.42151175831644588870D-7, &
             f7 = 2.04426310338993978564D-15
! HASH SUM EF    47.52583317549289671629

ifault = 0
q = p - half
IF (ABS(q) <= split1) THEN
  r = const1 - q * q
  normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / &
           (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
  RETURN
ELSE
  IF (q < zero) THEN
    r = p
  ELSE
    r = one - p
  END IF
  IF (r <= zero) THEN
    ifault = 1
    normal_dev = zero
    RETURN
  END IF
  r = SQRT(-LOG(r))
  IF (r <= split2) THEN
    r = r - const2
    normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / &
             (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
  ELSE
    r = r - split2
    normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / &
             (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
  END IF
  IF (q < zero) normal_dev = - normal_dev
  RETURN
END IF
RETURN
END SUBROUTINE ppnd16

END FUNCTION prtrng
