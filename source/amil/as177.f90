SUBROUTINE nscor1 (s, n, n2, work, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-17  Time: 12:13:20

!     Algorithm AS 177   Appl. Statist. (1982) Vol. 31, No. 2

!     Exact calculation of Normal Scores

IMPLICIT NONE
INTEGER, PARAMETER   :: dp = SELECTED_REAL_KIND(15, 100)

INTEGER, INTENT(IN)        :: n2
REAL (dp), INTENT(OUT)     :: s(n2)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: work(4,721)
INTEGER, INTENT(OUT)       :: ifault

INTERFACE
  SUBROUTINE init(work)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15, 100)
    REAL (dp), INTENT(OUT)  :: work(4,721)
  END SUBROUTINE init
END INTERFACE

REAL (dp) :: c, scor, ai1, ani, an

INTEGER   :: i, i1, j, ni
REAL (dp), PARAMETER  :: one = 1.0D0, zero = 0.0D0, h = 0.025D0
INTEGER, PARAMETER    :: nstep = 721

ifault=3
IF (n2 /= n/2) RETURN
ifault=1
IF (n <= 1) RETURN
ifault=0
IF (n > 2000) ifault=2

an=n
c=LOG(an)

!        Accumulate ordinates for calculation of integral for rankits

DO  i=1, n2
  i1=i-1
  ni=n-i
  ai1=i1
  ani=ni
  scor=zero
  DO  j=1,nstep
    scor = scor + EXP(work(2,j) + ai1*work(3,j) + ani*work(4,j) + c) * work(1,j)
  END DO
  s(i)=scor * h
  c=c + LOG(ani/i)
END DO
RETURN

END SUBROUTINE nscor1



SUBROUTINE init(work)

!        Algorithm AS 177.1   Appl. Statist. (1982) Vol. 31, No. 2

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15, 100)

INTERFACE
  FUNCTION alnorm( x, upper ) RESULT( fn_val )
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15, 100)
    REAL(dp), INTENT(IN)   ::  x
    LOGICAL, INTENT(IN)    ::  upper
    REAL(dp)               ::  fn_val
  END FUNCTION alnorm
END INTERFACE

REAL (dp), INTENT(OUT)  :: work(4,721)

REAL (dp)  :: xx
REAL (dp), PARAMETER  :: xstart = -9.0D0, h = 0.025D0, pi2 = -0.918938533D0, &
                         half = 0.5D0
INTEGER, PARAMETER    :: nstep = 721
INTEGER  :: i

xx=xstart

!        Set up arrays for calculation of integral

DO  i=1,nstep
  work(1,i)=xx
  work(2,i)=pi2 - xx * xx * half
  work(3,i)=LOG(alnorm(xx, .true.))
  work(4,i)=LOG(alnorm(xx, .false.))
  xx=xstart + i * h
END DO
RETURN
END SUBROUTINE init



!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 2001

FUNCTION alnorm( x, upper ) RESULT( fn_val )
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(15, 100)

REAL(DP), INTENT(IN)   ::  x
LOGICAL, INTENT(IN)    ::  upper
REAL(DP)               ::  fn_val

!  Local variables
REAL(DP), PARAMETER   ::  zero=0.0_DP, one=1.0_DP, half=0.5_DP, con=1.28_DP
REAL(DP)              ::  z, y
LOGICAL               ::  up

!  Machine dependent constants
REAL(DP), PARAMETER  ::  ltone = 7.0_DP, utzero = 18.66_DP
REAL(DP), PARAMETER  ::  p = 0.398942280444_DP, q = 0.39990348504_DP,   &
                         r = 0.398942280385_DP, a1 = 5.75885480458_DP,  &
                         a2 = 2.62433121679_DP, a3 = 5.92885724438_DP,  &
                         b1 = -29.8213557807_DP, b2 = 48.6959930692_DP, &
                         c1 = -3.8052E-8_DP, c2 = 3.98064794E-4_DP,     &
                         c3 = -0.151679116635_DP, c4 = 4.8385912808_DP, &
                         c5 = 0.742380924027_DP, c6 = 3.99019417011_DP, &
                         d1 = 1.00000615302_DP, d2 = 1.98615381364_DP,  &
                         d3 = 5.29330324926_DP, d4 = -15.1508972451_DP, &
                         d5 = 30.789933034_DP

up = upper
z = x
IF( z < zero ) THEN
   up = .NOT. up
   z = -z
END IF
IF( z <= ltone .OR. (up .AND. z <= utzero) ) THEN
   y = half*z*z
   IF( z > con ) THEN
      fn_val = r*EXP( -y )/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
   ELSE
      fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
   END IF
ELSE
   fn_val = zero
END IF

IF( .NOT. up ) fn_val = one - fn_val
RETURN
END FUNCTION alnorm



!----------------------------------------------------------------------------


SUBROUTINE nscor2(s, n, n2, ier)

!     Algorithm as 177.3, applied statistics, v.31, 161-165, 1982.

!     Calculates approximate expected values of normal order statistics.
!     claimed accuracy is 0.0001, though usually accurate to 5-6 dec.

!     Arguments:
!     s(n2)   = output, the first n2 expected values.
!     n       = input, the sample size.
!     n2      = input, the number of order statistics required; must
!                      be <= n/2.
!     ier     = output, error indicator
!                   = 0 if no error detected
!                   = 1 if n <= 1.
!                   = 2 if n > 2000, in which case the order statistics
!                          are still calculated, but may be inaccurate.
!                   = 3 if n2 > n/2 (n.b. this differs from the
!                          published algorithm which returns an error
!                          if n2 is not equal to n/2.)

!     Calls ppnd = Applied Statistics algorithm 111.
!     An alternative is ppnd7 in algorithm AS 241.

! N.B. This routine is not in double precision.

IMPLICIT NONE

INTEGER, INTENT(IN)   :: n2
REAL, INTENT(OUT)     :: s(n2)
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(OUT)  :: ier

INTERFACE
  FUNCTION correc(i, n) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: i
    INTEGER, INTENT(IN)  :: n
    REAL                 :: fn_val
  END FUNCTION correc

  SUBROUTINE ppnd7 (p, normal_dev, ifault)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: sp = SELECTED_REAL_KIND(6, 30)
    REAL (sp), INTENT(IN)   :: p
    INTEGER, INTENT(OUT)    :: ifault
    REAL (sp), INTENT(OUT)  :: normal_dev
  END SUBROUTINE ppnd7
END INTERFACE

REAL     :: an, ai, e1, e2, l1
INTEGER  :: i, k
REAL , PARAMETER  ::   &
    eps(4) = (/ 0.419885, 0.450536, 0.456936, 0.468488 /),  &
    dl1(4) = (/ 0.112063, 0.121770, 0.239299, 0.215159 /),  &
    dl2(4) = (/ 0.080122, 0.111348, -0.211867, -0.115049 /),  &
    gam(4) = (/ 0.474798, 0.469051, 0.208597, 0.259784 /),  &
    lam(4) = (/ 0.282765, 0.304856, 0.407708, 0.414093 /),  &
    bb = -0.283833, d = -0.106136, b1 = 0.5641896

!     input parameter checks.

ier = 3
IF(n2 > n/2) RETURN
ier = 1
IF(n <= 1) RETURN
ier = 0
IF(n > 2000) ier = 2
s(1) = b1
IF(n == 2) RETURN

!     calculate normal tail areas for first 3 order statistics.

an = n
k = 3
IF(n2 < k) k = n2
DO  i = 1,k
  ai = i
  e1 = (ai - eps(i))/(an + gam(i))
  e2 = e1**lam(i)
  s(i) = e1 + e2*(dl1(i) + e2*dl2(i))/an - correc(i,n)
END DO
IF(n2 == k) GO TO 20

!     calculate normal areas for other cases.

DO  i = 4,n2
  ai = i
  l1 = lam(4) + bb/(ai + d)
  e1 = (ai - eps(4))/(an + gam(4))
  e2 = e1**l1
  s(i) = e1 + e2*(dl1(4) + e2*dl2(4))/an - correc(i,n)
END DO

!     convert tail areas to normal deviates.

20 DO  i = 1,n2
  CALL ppnd7(s(i), s(i), ier)
END DO

RETURN
END SUBROUTINE nscor2



FUNCTION correc(i, n) RESULT(fn_val)

! Calculates correction for tail area of the i-th largest of n order statistics.

IMPLICIT NONE

INTEGER, INTENT(IN)  :: i
INTEGER, INTENT(IN)  :: n
REAL                 :: fn_val

REAL  :: an
REAL, PARAMETER  :: c1(7) = (/ 9.5, 28.7, 1.9, 0., -7.0, -6.2, -1.6 /),  &
    c2(7) = (/ -6195., -9569., -6728., -17614., -8278., -3570., 1075. /),  &
    c3(7) = (/ 9.338E4, 1.7516E5, 4.1040E5, 2.1576E6, 2.376E6, 2.065E6,  &
    2.065E6 /), mic = 1.e-6, c14 = 1.9E-5

fn_val = c14
IF(i*n == 4) RETURN
fn_val = 0.0
IF(i < 1 .OR. i > 7) RETURN
IF(i /= 4 .AND. n > 20) RETURN
IF(i == 4 .AND. n > 40) RETURN
an = n
an = 1.0/(an*an)
fn_val = (c1(i) + an*(c2(i) + an*c3(i)))*mic
RETURN
END FUNCTION correc



SUBROUTINE ppnd7 (p, normal_dev, ifault)

! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477- 484.

! Produces the normal deviate Z corresponding to a given lower tail area of P;
! Z is accurate to about 1 part in 10**7.

! The hash sums below are the sums of the mantissas of the coefficients.
! They are included for use in checking transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

IMPLICIT NONE
INTEGER, PARAMETER  :: sp = SELECTED_REAL_KIND(6, 30)

REAL (sp), INTENT(IN)   :: p
INTEGER, INTENT(OUT)    :: ifault
REAL (sp), INTENT(OUT)  :: normal_dev

! Local variables

REAL (sp) :: zero = 0.0, one = 1.0, half = 0.5, split1 = 0.425,  &
             split2 = 5.0, const1 = 0.180625, const2 = 1.6, q, r

! Coefficients for P close to 0.5

REAL (sp) :: a0 = 3.3871327179E+00, a1 = 5.0434271938E+01, &
                   a2 = 1.5929113202E+02, a3 = 5.9109374720E+01, &
                   b1 = 1.7895169469E+01, b2 = 7.8757757664E+01, &
                   b3 = 6.7187563600E+01
! HASH SUM AB          32.3184577772

! Coefficients for P not close to 0, 0.5 or 1.

REAL (sp) :: c0 = 1.4234372777E+00, c1 = 2.7568153900E+00, &
             c2 = 1.3067284816E+00, c3 = 1.7023821103E-01, &
             d1 = 7.3700164250E-01, d2 = 1.2021132975E-01
! HASH SUM CD    15.7614929821

! Coefficients for P near 0 or 1.

REAL (sp) :: e0 = 6.6579051150E+00, e1 = 3.0812263860E+00, &
             e2 = 4.2868294337E-01, e3 = 1.7337203997E-02, &
             f1 = 2.4197894225E-01, f2 = 1.2258202635E-02
! HASH SUM EF    19.4052910204

ifault = 0
q = p - half
IF (ABS(q) <= split1) THEN
  r = const1 - q * q
  normal_dev = q * (((a3 * r + a2) * r + a1) * r + a0) / &
               (((b3 * r + b2) * r + b1) * r + one)
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
    normal_dev = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one)
  ELSE
    r = r - split2
    normal_dev = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one)
  END IF
  IF (q < zero) normal_dev = - normal_dev
  RETURN
END IF
RETURN
END SUBROUTINE ppnd7
