FUNCTION random_Poisson(mu, first) RESULT(ran_Poisson)
!**********************************************************************

!     Translated to Fortran 90 by Alan Miller from:
!                                    RANLIB
!
!           Library of Fortran Routines for Random Number Generation
!
!                           Compiled and Written by:
!
!                                Barry W. Brown
!                                 James Lovato
!
!                    Department of Biomathematics, Box 237
!                    The University of Texas, M.D. Anderson Cancer Center
!                    1515 Holcombe Boulevard
!                    Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                 GENerate POIsson random deviate Function

!  Generates a single random deviate from a Poisson
!  distribution with mean mu.

!                           Arguments

!  mu --> The mean of the Poisson distribution from which
!         a random deviate is to be generated.
!                           REAL mu

!                           Method

!  For details see:

!      Ahrens, J.H. and Dieter, U.
!      Computer Generation of Poisson Deviates
!      From Modified Normal Distributions.
!      ACM Trans. Math. Software, 8, 2 (June 1982), 163-179


!  MUPREV = PREVIOUS MU, MUOLD = MU AT LAST EXECUTION OF STEP P OR B.
!  TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
!  COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K) - DEL

!  SEPARATION OF CASES A AND B

IMPLICIT NONE

!     .. Scalar Arguments ..
REAL, INTENT(IN)    :: mu
LOGICAL, INTENT(IN) :: first
INTEGER             :: ran_Poisson

INTERFACE
  FUNCTION random_normal() RESULT (ran_norm)
    IMPLICIT NONE
    REAL :: ran_norm
  END FUNCTION random_normal
END INTERFACE

INTERFACE
  FUNCTION random_exponential() RESULT (ran_exp)
    IMPLICIT NONE
    REAL :: ran_exp
  END FUNCTION random_exponential
END INTERFACE

!     .. Local Scalars ..
REAL    :: a0 = -.5, a1 = .3333333, a2 = -.2500068, a3 = .2000118,        &
           a4 = -.1661269, a5 = .1421878, a6 = -.1384794, a7 = .1250060,  &
           b1, b2, c, c0, c1, c2, c3, d, del, difmuk, e, fk, fx, fy, g,   &
           omega, p, p0, px, py, q, s, t, u, v, x, xx
INTEGER :: j, k, kflag, l, m
!     ..
!     .. Local Arrays ..
REAL    :: fact(10) = (/ 1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880. /)
REAL    :: pp(35)
!     ..
!     ..
!     .. Executable Statements ..
IF (.NOT. first) GO TO 10
IF (mu < 10.0) GO TO 120

!     C A S E  A. (RECALCULATION OF S, D, L IF MU HAS CHANGED)

s = SQRT(mu)
d = 6.0*mu*mu

!     THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL PROBABILITIES FK
!     WHENEVER K >= M(MU).
!     L = INT(MU-1.1484) IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .

l = INT(mu-1.1484)

!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE

10 g = mu + s*random_normal()
IF (g < 0.0) GO TO 20
ran_Poisson = INT(g)

!     STEP I. IMMEDIATE ACCEPTANCE IF ran_Poisson IS LARGE ENOUGH

IF (ran_Poisson.GE.l) RETURN

!     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U

fk = REAL(ran_Poisson)
difmuk = mu - fk
CALL RANDOM_NUMBER(u)
IF (d*u >= difmuk*difmuk*difmuk) RETURN

!     STEP P. PREPARATIONS FOR STEPS Q AND H.
!             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
!             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
!             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
!             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
!             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.

20 IF (.NOT. first) GO TO 30
omega = .3989423/s
b1 = .4166667E-1 / mu
b2 = .3*b1*b1
c3 = .1428571*b1*b2
c2 = b2 - 15.*c3
c1 = b1 - 6.*b2 + 45.*c3
c0 = 1. - b1 + 3.*b2 - 15.*c3
c = .1069/mu
30 IF (g < 0.0) GO TO 50

!             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)

kflag = 0
GO TO 70

!     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)

40 IF (fy-u*fy <= py*EXP(px-fx)) RETURN

!     STEP E. EXPONENTIAL SAMPLE - random_expon() FOR STANDARD EXPONENTIAL
!             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
!             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)

50 e = random_exponential()
CALL RANDOM_NUMBER(u)
u = u + u - 1.0
t = 1.8 + SIGN(e, u)
IF (t <= (-.6744)) GO TO 50
ran_Poisson = INT(mu+s*t)
fk = REAL(ran_Poisson)
difmuk = mu - fk

!             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)

kflag = 1
GO TO 70

!     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)

60 IF (c*ABS(u) > py*EXP(px+e) - fy*EXP(fx+e)) GO TO 50
RETURN

!     STEP F. 'SUBROUTINE' F. CALCULATION OF PX, PY, FX, FY.
!             CASE ran_Poisson < 10 USES FACTORIALS FROM TABLE FACT

70 IF (ran_Poisson >= 10) GO TO 80
px = -mu
py = mu**ran_Poisson / fact(ran_Poisson + 1)
GO TO 110

!      CASE ran_Poisson >= 10 USES POLYNOMIAL APPROXIMATION
!      A0-A7 FOR ACCURACY WHEN ADVISABLE
!      .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)

80 del = .8333333E-1 / fk
del = del - 4.8*del*del*del
v = difmuk/fk
IF (ABS(v) <= 0.25) GO TO 90
px = fk*LOG(1.0+v) - difmuk - del
GO TO 100

90 px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
100 py = .3989423/SQRT(fk)
110 x = (0.5-difmuk)/s
xx = x*x
fx = -0.5*xx
fy = omega* (((c3*xx+c2)*xx+c1)*xx+c0)
IF (kflag > 0) THEN
  GO TO 60
ELSE
  GO TO 40
END IF

!     C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)

120 IF (first) THEN
  m = MAX(1, INT(mu))
  l = 0
  p = EXP(-mu)
  q = p
  p0 = p
END IF

!     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD

130 CALL RANDOM_NUMBER(u)
ran_Poisson = 0
IF (u <= p0) RETURN

!     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
!             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
!             (0.458 = PP(9) FOR MU = 10)

IF (l == 0) GO TO 150
j = 1
IF (u > 0.458) j = MIN(l, m)
DO k = j, l
  IF (u <= pp(k)) GO TO 180
END DO
IF (l == 35) GO TO 130

!     STEP C. CREATION OF NEW POISSON PROBABILITIES P
!             AND THEIR CUMULATIVES Q = PP(K)

150 l = l + 1
DO k = l, 35
  p = p*mu/REAL(k)
  q = q + p
  pp(k) = q
  IF (u <= q) GO TO 170
END DO
l = 35
GO TO 130

170 l = k
180 ran_Poisson = k
RETURN

END FUNCTION random_Poisson
