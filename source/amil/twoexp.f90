PROGRAM twoexp

! Code converted using TO_F90 by Alan Miller
! Date: 1999-03-15  Time: 11:15:27

!           SAMPLE DRIVER PROGRAM FOR NONLINEAR LEAST SQUARES ROUTINE

!            SOLUTION TO SAMPLE PROBLEM (TO 3 SIGNIFICANT DIGITS):
!               LINEAR PARAMETERS  .375, 1.94, -1.46
!               NONLINEAR PARAMETERS  .0129, .0221

USE separable_leastsq
IMPLICIT NONE
REAL (dp) :: y(50), t(50,1), alf(12), beta(8), a(50,13), w(50), eta
INTEGER   :: i, ierr, im1, iprint, iv, j, l, lpnl, n, nl, nmax, p
INTEGER, PARAMETER :: INPUT = 5, output = 6

INTERFACE
  SUBROUTINE ada (lp1, nl, n, nmax, lpp2, iv, a, inc, t, alf, isel)
    USE separable_leastsq
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: lp1
    INTEGER, INTENT(IN)     :: nl
    INTEGER, INTENT(IN)     :: n
    INTEGER, INTENT(IN)     :: nmax
    INTEGER, INTENT(IN)     :: lpp2
    INTEGER, INTENT(IN)     :: iv
    REAL (dp), INTENT(OUT)  :: a(:,:)
    INTEGER, INTENT(OUT)    :: inc(:,:)
    REAL (dp), INTENT(IN)   :: t(:,:)
    REAL (dp), INTENT(IN)   :: alf(:)
    INTEGER, INTENT(IN)     :: isel
  END SUBROUTINE ada
END INTERFACE

!              READ AND LIST THE DATA

OPEN (INPUT, FILE='a.dat', STATUS='OLD')
OPEN (output, FILE='b.out', STATUS='NEW')
nmax = 50
iprint = 1
READ (INPUT,*) n, l, nl, p, iv
WRITE (output,110) n, l, nl, p, iv
DO  i=1,n
  w(i) = 1.0_dp
  READ (INPUT,*) t(i,1:iv), y(i)
END DO
DO  i=1,n
  WRITE (output,210) i, t(i,1:iv), y(i)
END DO

READ (INPUT,120) alf(1:nl)
WRITE (output,130) alf(1:nl)
WRITE (output,140)

CALL varpro (l, nl, n, nmax, l+p+2, iv, t, y, w, ada, a, iprint,  &
             alf, beta, ierr)
WRITE (6,30) ierr
30 FORMAT (t16, 'IERR =', i7)


!              PRINT LOWER TRIANGLE OF COVARIANCE MATRIX

IF (ierr <= -4) GO TO 90
lpnl = l + nl
WRITE (output,160)
DO  i=1,lpnl
  WRITE (output,170) a(i,1:i)
  a(i,i) = SQRT(a(i,i))
END DO

!              PRINT STANDARD DEVIATIONS OF PARAMETER ESTIMATES

WRITE (output,180)
WRITE (output,170) (a(j,j),j=1,lpnl)

!              CALCULATE AND PRINT LOWER TRIANG. OF CORRELATION MATRIX

IF (lpnl == 1) GO TO 60
DO  i=2,lpnl
  im1 = i-1
  DO  j=1,im1
    a(i,j) = a(i,j) / (a(i,i)*a(j,j))
  END DO
END DO

60 WRITE (output,190)
DO  i=1,lpnl
  a(i,i) = 1.0_dp
  WRITE (output,170) a(i,1:i)
END DO

!              PRINT RESIDUALS

WRITE (output,200)
DO  i=1,n
  eta = y(i) - a(i,lpnl+1) / SQRT(w(i))
  WRITE (output,210) i, w(i), t(i,1:iv), y(i), eta, a(i,lpnl+1)
END DO

90 STOP

110 FORMAT ('1   NON-LINEAR LEAST SQUARES PROBLEM'//  &
    ' NUMBER OF OBSERVATIONS =',i5, '   NUMBER OF LINEAR PARAMETERS =',i4//  &
    ' NUMBER OF NONLINEAR PARAMETERS =',i4,  &
    ' NUMBER OF NONVANISHING PARTIAL DERIVATIVES =',i4//  &
    ' NUMBER OF INDEPENDENT VARIABLES =',i4// '    I    T(I)            Y(I)'//)
120 FORMAT (4E20.7)
130 FORMAT ('0 INITIAL NONLINEAR PARAMETERS'//(4E20.7))
140 FORMAT ('0',50('*'))
160 FORMAT (/'       COVARIANCE MATRIX'/)
170 FORMAT (8E15.7)
180 FORMAT ('0    STANDARD DEVIATIONS OF PARAMETER ESTIMATES'/)
190 FORMAT (/'       CORRELATION MATRIX'/)
200 FORMAT (/'   I         W(I)            T(I)            Y(I)',  &
    '     PREDICTED Y     WEIGHTED RESIDUAL'//)
210 FORMAT (i5,7E16.7)
END PROGRAM twoexp



SUBROUTINE ada (lp1, nl, n, nmax, lpp2, iv, a, inc, t, alf, isel)

!        EVALUATE THE FUNCTIONS PHI AND THEIR PARTIAL DERIVATIVES
!        D PHI(J)/D ALF(K), AT THE SAMPLE POINTS T(I).
!        ISEL = 1  MEANS COMPUTE BOTH FUNCTIONS AND DERIVATIVES AND
!                  INITIALIZE INC AND CONSTANT PHI'S
!             = 2  MEANS COMPUTE ONLY THE NONCONSTANT FUNCTIONS PHI
!             = 3  MEANS COMPUTE ONLY THE DERIVATIVES.

!        THIS PARTICULAR ROUTINE IS FOR FITTING TWO EXPONENTIALS AND
!        ONE CONSTANT TERM:

!        ETA(C, ALF; T) = C  + C  * EXP(-ALF *T) + C * EXP(-ALF * T)
!                          1    2           1       3          2

!        WHERE C IS THE VECTOR OF LINEAR PARAMETERS TO BE DETERMINED,
!        AND ALF IS THE VECTOR OF NONLINEAR PARAMETERS TO BE DETERMINED.
!        HERE N = 33, L = 3, NL = P = 2, NCFUN = 1, PHI(1) = 1, PHI(2) =
!        EXP(-ALF(1)*T), AND PHI(3) = EXP(-ALF(2)*T).  THIS EXAMPLE AND
!        TEST DATA ARE TAKEN FROM REFERENCE 3 (SEE VARPRO).

!     ..................................................................

USE separable_leastsq
IMPLICIT NONE

INTEGER, INTENT(IN)     :: lp1
INTEGER, INTENT(IN)     :: nl
INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: nmax
INTEGER, INTENT(IN)     :: lpp2
INTEGER, INTENT(IN)     :: iv
REAL (dp), INTENT(OUT)  :: a(:,:)        ! a(nmax,lpp2)
INTEGER, INTENT(OUT)    :: inc(:,:)      ! inc(12,lp1)
REAL (dp), INTENT(IN)   :: t(:,:)        ! t(nmax)
REAL (dp), INTENT(IN)   :: alf(:)        ! alf(nl)
INTEGER, INTENT(IN)     :: isel

! Local variable

INTEGER :: i

SELECT CASE ( isel )
  CASE (1)
    GO TO 10
  CASE (2)
    GO TO 30
  CASE (3)
    GO TO 50
END SELECT

!              ISEL = 1, IN THIS CASE THE INCIDENCE MATRIX INC IS:

!               0  1  0  0     <- ALF(1)
!               0  0  1  0     <- ALF(2)

10 inc = 0
inc(1,2) = 1
inc(2,3) = 1

!              ISEL = 1, SET CONSTANT FUNCTION PHI(1) (REMOVE IF NO
!              CONSTANT FUNCTIONS)

a(1:n,1) = 1.0_dp

!              ISEL = 1 OR 2, COMPUTE NONCONSTANT FUNCTIONS PHI

30 DO  i=1,n
  a(i,2) = EXP(-alf(1)*t(i,1))
  a(i,3) = EXP(-alf(2)*t(i,1))
END DO
!                           COLUMN L+1 = 4 IS LEFT FOR WORKSPACE.
IF (isel == 2) GO TO 70

!              ISEL = 1 OR 3, COMPUTE DERIVATIVES DPHI(I) / D ALF(J)

50 DO  i=1,n
  a(i,5) = -t(i,1)*EXP(-alf(1)*t(i,1))
  a(i,6) = -t(i,1)*EXP(-alf(2)*t(i,1))
END DO

70 RETURN
END SUBROUTINE ada
