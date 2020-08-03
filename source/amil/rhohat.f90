SUBROUTINE rhohat(rho, y, tol)

!*********************************************************************

!     Find RHOHAT (ML estimator) of gamma density by Y

!     Requires subroutine POLG(X, PSI, POL, N, ICODE)

!     INPUT:  Y = LN(XBAR / Geometric mean) = LN(RHO) - PSI(RHO)
!             TOL = desired accuracy

!     OUTPUT: RHO = ML estimate of RHO

! N.B. The scale parameter for the gamma density is then obtained by
!      equating the population (arithmetic) mean to XBAR = sample mean.

! This version, which is compatible with Lahey's ELF90 compiler,
! is by Alan Miller: amiller @ bigpond.net.au
! URL: http://users.bigpond.net.au/amiller
! Latest revision - 27 August 1997

!*********************************************************************
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
REAL (dp), INTENT(IN)  :: y, tol
REAL (dp), INTENT(OUT) :: rho

! Local variables
REAL (dp)              :: r1, di
REAL (dp), PARAMETER   :: one = 1._dp, four = 4._dp, three = 3._dp
INTEGER                :: i

INTERFACE
  SUBROUTINE polg(x, psi, pol, n, icode)
    IMPLICIT NONE
    INTEGER, PARAMETER               :: dp = SELECTED_REAL_KIND(15, 60)
    REAL (dp), INTENT(IN)            :: x
    INTEGER, INTENT(IN), OPTIONAL    :: n
    REAL (dp), INTENT(OUT)           :: psi
    REAL (dp), INTENT(OUT), OPTIONAL :: pol(:)
    INTEGER, INTENT(OUT), OPTIONAL   :: icode
  END SUBROUTINE polg
END INTERFACE

IF (y <= EPSILON(one)) THEN
  WRITE(*, *) 'Error: argument y to routine RHOHAT must be > 0'
  RETURN
END IF

!     Initial RHO is Thom's estimator
rho = (one + SQRT(one + four * y / three)) / (four * y)

DO i = 1, 100
  r1 = rho
  CALL polg(r1, di)
  rho = r1 * (LOG(r1) - di) / y
  IF (ABS(r1 - rho) < tol) EXIT
END DO

RETURN
END SUBROUTINE rhohat



SUBROUTINE polg(x, psi, pol, n, icode)

!     This program computes PSI(X) & its first N derivatives.

!*********************************************************************

!     Input parameter X must be positive and < 4000
!                     N integer - highest derivative required >= 0
!     Output parameter PSI = 1st derivative of log (gamma (X) )
!                      POL dimension at least N contains the derivatives
!                      ICODE = 0 for normal exit
!                            = 1 if X <= 0
!                            = 2 if N < 0 or N > 20

!*********************************************************************

!     BER(S) = the Bernouilli number of order 2S (apart from signs)

IMPLICIT NONE
INTEGER, PARAMETER               :: dp = SELECTED_REAL_KIND(15, 60)
REAL (dp), INTENT(IN)            :: x
INTEGER, INTENT(IN), OPTIONAL    :: n
REAL (dp), INTENT(OUT)           :: psi
REAL (dp), INTENT(OUT), OPTIONAL :: pol(:)
INTEGER, INTENT(OUT), OPTIONAL   :: icode

! Local variables
REAL (dp)              :: x2, w(8), di, dj, dum, ds, sgn
REAL (dp), PARAMETER   :: ber(8) = (/ 0.16666666666666667D+00,  &
   0.33333333333333333D-01, 0.23809523809523810D-01, 0.33333333333333333D-01, &
   0.75757575757575758D-01, 0.25311355311355311D+00, 0.11666666666666667D+01, &
   0.70921568627450980D+01 /)
!     Values above are 1/6, 1/30, 1/42, 1/30, 5/66, 691/2730, 7/6 & 3617/510.

INTEGER                :: i, ic, in, indx, j, jj, nn

!     Check for valid input parameters X and N.

ic = 0
IF (x < 0.01D0 .OR. x > 4000.D0) ic = 1
IF (PRESENT(n)) THEN
  IF (n < 0 .OR. n > 20) ic = 2
  nn = n
ELSE
  nn = 0
END IF
IF (PRESENT(icode)) icode = ic
IF (ic > 0) RETURN

!     Increment X to a value > 40

indx = 0
x2 = x
DO i = 1, 40
  IF (x2 > 40.D0) EXIT
  indx = indx + 1
  x2 = x2 + 1.D0
END DO

IF (x2 > 500.D0) THEN
  in = 3
ELSE IF (x2 > 200.D0) THEN
  in = 4
ELSE IF (x2 > 100.D0) THEN
  in = 6
ELSE
  in = 8
END IF

!     Calculate polygamma function using table 6.4.11 on page 260 of
!     Bowman, K.O. & Shenton, L.R. (1988), Properties of Estimators for
!     the Gamma Distribution, Marcel Dekker.
!     X2 = 40 + fractional part of X

IF (nn == 0) GO TO 40
dum = 1.D0 / x2
sgn = 1.D0
DO i = 1, nn
  di = i
  pol(i) = dum * (1.D0 + di * 0.5D0 / x2)
  dum = - dum * di / x2
  DO j = 1, in
    dj = j
    IF (i == 1) THEN
      w(j) = sgn * ber(j) / x2 ** (2*j+1)
      sgn = - sgn
    ELSE
      w(j) = w(j) * (1.D0 - di - 2.D0 * dj) / x2
    END IF
    pol(i) = pol(i) + w(j)
  END DO
END DO

40 psi = LOG(x2) - 0.5D0 / x2
sgn = - 1.D0
DO i = 1, in
  di = 2 * i
  dum = ber(i) / x2 ** (2*i) / di
  psi = psi + sgn * dum
  sgn = - sgn
END DO

!     Calculate desirable polygamma function using 6.4.6 from page 260.

IF (indx == 0) RETURN
IF (PRESENT(pol)) THEN
  sgn = 1.D0
  DO i = 1, nn
    ds = 0.D0
    DO j = 1, indx
      dj = indx - j
      dum = 1.D0 / (x + dj)
      dj = dum
      DO jj = 1, i
        di = jj
        dj = dj * di * dum
      END DO
      ds = ds + dj
    END DO
    pol(i) = pol(i) + sgn * ds
    sgn = - sgn
  END DO
END IF

DO i = 1, indx
  dum = i - 1
  psi = psi - 1.D0 / (x + dum)
END DO

RETURN
END SUBROUTINE polg


PROGRAM drive
IMPLICIT NONE

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 60)
REAL (dp)          :: y, rho, tol = 1.D-04

INTERFACE
  SUBROUTINE rhohat(rho, y, tol)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
    REAL (dp), INTENT(IN)  :: y, tol
    REAL (dp), INTENT(OUT) :: rho
  END SUBROUTINE rhohat
END INTERFACE

DO
  WRITE(*, '(a)', ADVANCE='NO') ' Enter value of y: '
  READ(*, *) y
  CALL rhohat(rho, y, tol)
  WRITE(*, '(a, f10.6)') ' RHO = ', rho
  WRITE(*, *)
END DO

STOP
END PROGRAM drive
