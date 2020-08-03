SUBROUTINE lfnorm(n, m, x, y, beta, z, ky, ifault)

!     ALGORITHM AS 135  APPL. STATIST. (1979) VOL.28, NO. 1

!     Min-Max estimates for a linear multiple regression problem

! Input arguments
! ---------------
! N      Number of cases
! M      Number of predictors
! X(,)   Values of the predictor variables, one case per row.
! Y()    Values of the dependent variable

! Output arguments
! ----------------
! BETA() Final estimates of the regression coefficients
! Z      The least maximum absolute deviation
! KY     The number of iterations
! IFAULT Error indicator = 0 for normal termination
!                        = 1 if the observation matrix has less than full rank
!
! N.B. Arguments NDIM & MDIM have been removed.

! Fortran 90 version by Alan.Miller @ vic.cmis.csiro.au
! Latest revision - 5 January 1999

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: x(:, :)   ! x(ndim, mdim)
REAL (dp), INTENT(IN)   :: y(:)
REAL (dp), INTENT(OUT)  :: beta(:)
REAL (dp), INTENT(OUT)  :: z
INTEGER, INTENT(OUT)    :: ky
INTEGER, INTENT(OUT)    :: ifault

REAL (dp) :: hilo(m), xrxf(m), xsxf(m), lu(m, m), dev1, sigr, sumxr,  &
             deviat, yest, ydev, sumxs, ratio, test, delta, div, swing,    &
             SAVE, sigs, top, acu, big
INTEGER   :: i, ibase(m), ii, ii1, irow, k, kk, kkk, k1, k2, m1, sss, rrr,  &
             indx(m)
LOGICAL   :: intl

REAL (dp), PARAMETER :: zero = 0._dp, one = 1._dp, two = 2._dp

acu = SQRT( EPSILON(1.0_dp) )
big = SQRT( HUGE(1.0_dp) )
ifault = 0
ky = 0
z = zero
m1 = m - 1

!     Set up initial LU decomposition

DO i = 1, m
  indx(i) = i
END DO
intl = .true.
kkk = 1
CALL update(kkk, x, lu, ibase, indx, intl, n, m, ifault)
IF (ifault /= 0) RETURN
intl = .false.
irow = kkk

!     Calculate beta value

k = indx(1)
k1 = ibase(1)
beta(k) = y(k1) / lu(k, 1)
DO ii = 2, m
  k = indx(ii)
  k1 = ibase(ii)
  beta(k) = y(k1)
  ii1 = ii - 1
  DO i = 1, ii1
    kk = indx(i)
    beta(k) = beta(k) - lu(kk, ii) * beta(kk)
  END DO
  beta(k) = beta(k) / lu(k, ii)
END DO
DO ii = 1, m1
  k1 = m - ii
  k = indx(k1)
  DO i = 1, ii
    kk = m - i + 1
    k2 = indx(kk)
    beta(k) = beta(k) - lu(k2, k1) * beta(k2)
  END DO
END DO

!     Search for and set first violated constraint as R-th constraint

50 irow = irow + 1
IF (irow > n) RETURN
dev1 = zero
DO i = 1, m
  dev1 = dev1 + x(irow, i) * beta(i)
END DO
dev1 = dev1 - y(irow)
IF (ABS(dev1) < acu) GO TO 50
sigr = SIGN(one, dev1)
rrr = irow

!     Adjust for the R-th constraint

k = indx(1)
xrxf(1) = x(rrr, k)
DO ii = 2, m
  k = indx(ii)
  xrxf(ii) = x(rrr, k)
  ii1 = ii - 1
  DO i = 1, ii1
    xrxf(ii) = xrxf(ii) - lu(k, i) * xrxf(i)
  END DO
END DO
k = indx(m)
xrxf(m) = xrxf(m) / lu(k, m)
hilo(m) = SIGN(one, -sigr * xrxf(m))
sumxr = sigr - hilo(m) * xrxf(m)
DO ii = 1, m1
  k1 = m - ii
  k = indx(k1)
  DO i = 1, ii
    k2 = m - i + 1
    xrxf(k1) = xrxf(k1) - lu(k, k2) * xrxf(k2)
  END DO
  xrxf(k1) = xrxf(k1) / lu(k, k1)
  hilo(k1) = SIGN(one, -sigr * xrxf(k1))
  sumxr = sumxr - hilo(k1) * xrxf(k1)
END DO
z = ABS(dev1 / sumxr)

!     Start of main iterative loop.
!     Search for the most violated S-th constraint

110 sss = 0
deviat = acu

!     Calculate beta value

k = indx(1)
k1 = ibase(1)
beta(k) = (y(k1) + z * hilo(1)) / lu(k, 1)
DO ii = 2, m
  k = indx(ii)
  k1 = ibase(ii)
  beta(k) = y(k1) + z * hilo(ii)
  ii1 = ii - 1
  DO i = 1, ii1
    kk = indx(i)
    beta(k) = beta(k) - lu(kk, ii) * beta(kk)
  END DO
  beta(k) = beta(k) / lu(k, ii)
END DO
DO ii = 1, m1
  k1 = m - ii
  k = indx(k1)
  DO i = 1, ii
    kk = m - i + 1
    k2 = indx(kk)
    beta(k) = beta(k) - lu(k2, k1) * beta(k2)
  END DO
END DO

!     Calculate residuals

DO i = 1, n
  yest = DOT_PRODUCT( x(i, 1:m), beta(1:m) )
  dev1 = ABS(y(i) - yest) - z
  IF (dev1 <= deviat) CYCLE
  ydev = yest - y(i)
  deviat = dev1
  sss = i
END DO

!     Check if at optimum

IF (sss == 0) RETURN

!     Set up information on the S-th constraint

sigs = SIGN(one, ydev)
k = indx(1)
xsxf(1) = x(sss, k)
DO ii = 2, m
  k = indx(ii)
  xsxf(ii) = x(sss, k)
  ii1 = ii - 1
  DO i = 1, ii1
    xsxf(ii) = xsxf(ii) - lu(k, i) * xsxf(i)
  END DO
END DO
k = indx(m)
xsxf(m) = xsxf(m) / lu(k, m)
sumxs = -sigs + hilo(m) * xsxf(m)
DO ii = 1, m1
  k1 = m - ii
  k = indx(k1)
  DO i = 1, ii
    k2 = m - i + 1
    xsxf(k1) = xsxf(k1) - lu(k, k2) * xsxf(k2)
  END DO
  xsxf(k1) = xsxf(k1) / lu(k, k1)
  sumxs = sumxs + hilo(k1) * xsxf(k1)
END DO

!     Search for minimum ratio

210 kkk = 0
ratio = big
DO i = 1, m
  IF (sigs * SIGN(one, xsxf(i)) /= hilo(i) .OR. ABS(xsxf(i)) < acu) CYCLE
  test = ABS(xrxf(i) / xsxf(i))
  IF (test >= ratio) CYCLE
  ratio = test
  kkk = i
END DO

!     Check if R-th constraint moves interior

IF (kkk /= 0) GO TO 260

!     Process the movement of the R-th constraint

delta = ABS(deviat / sumxs)

!     Calculate the largest tolerable delta

div = ABS(sumxr) - two
IF (div < acu) GO TO 240
swing = two * z / div
IF (swing >= delta) GO TO 240

!     Switch R and S constraint indicators

SAVE = sumxs
sumxs = -sumxr + sigr + sigr
sumxr = -SAVE
SAVE = sigr
sigr = sigs
sigs = -SAVE
deviat = ABS(sumxs * delta) - two * z
z = z + delta
DO i = 1, m
  SAVE = xsxf(i)
  xsxf(i) = xrxf(i)
  xrxf(i) = SAVE
END DO
i = rrr
rrr = sss
sss = i
GO TO 210

!     Replace the R-th constraint with the S-th constraint

240 sigr = sigs
xrxf(1:m) = xsxf(1:m)
sumxr = -sumxs
z = z + delta
rrr = sss
GO TO 110

!     Process the movement of the K-th constraint

260 delta = ABS(xrxf(kkk) * deviat / (xrxf(kkk) * sumxs + xsxf(kkk) * sumxr))
top = -two * z * xrxf(kkk)
div = xrxf(kkk) * xrxf(kkk) + hilo(kkk) * sumxr
IF (SIGN(one, top) /= SIGN(one, div)) GO TO 270
IF (ABS(div) < acu) GO TO 270
swing = top / div

!     Check to see if the K-th constraint swings across

IF (swing >= delta) GO TO 270
z = z + swing
deviat = deviat - swing * ABS(sumxs + xsxf(kkk) * sumxr / xrxf(kkk))
sumxr = sumxr + two * hilo(kkk) * xrxf(kkk)
sumxs = sumxs - two * hilo(kkk) * xsxf(kkk)
hilo(kkk) = -hilo(kkk)
GO TO 210

!     Update XRXF and the LU of the current basis

270 hilo(kkk) = sigs
sumxr = sigr
xrxf(kkk) = xrxf(kkk) / xsxf(kkk)
sumxr = sumxr - hilo(kkk) * xrxf(kkk)
DO i = 1, m
  IF (i == kkk) CYCLE
  xrxf(i) = xrxf(i) - xsxf(i) * xrxf(kkk)
  sumxr = sumxr - hilo(i) * xrxf(i)
END DO
ibase(kkk) = sss

!     Update LU decomposition

CALL update(kkk, x, lu, ibase, indx, intl, n, m, ifault)
IF (ifault /= 0) RETURN
z = z + delta
ky = ky + 1
GO TO 110

CONTAINS


SUBROUTINE update(kkk, x, lu, ibase, indx, intl, n, m, ifault)

!     ALGORITHM AS 135.1  APPL. STATIST. (1979) VOL.28, NO. 1

!     Update LU decomposition matrix

INTEGER, INTENT(IN OUT)   :: kkk
REAL (dp), INTENT(IN)     :: x(:, :)   ! x(ndim, mdim)
REAL (dp), INTENT(IN OUT) :: lu(:, :)
INTEGER, INTENT(IN OUT)   :: ibase(:)
INTEGER, INTENT(IN OUT)   :: indx(:)
LOGICAL, INTENT(IN OUT)   :: intl
INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(IN)       :: m
INTEGER, INTENT(OUT)      :: ifault

! Local variables

REAL (dp) :: subt, pivot
INTEGER   :: i, icol, ii, ii1, irow, isave, j, k, kk

irow = 0
DO ii = kkk, m
  IF (intl) GO TO 10
  irow = ibase(ii)
  GO TO 20

  10 irow = irow + 1
  IF (irow > n) THEN
    ifault = 1
    RETURN
  END IF
  20 lu(1:m, ii) = x(irow, 1:m)
  
!     Set up representation of incoming row
  
  IF (ii == 1) GO TO 60
  ii1 = ii - 1
  DO icol = 1, ii1
    k = indx(icol)
    subt = lu(k, ii)
    j = icol + 1
    DO i = j, m
      k = indx(i)
      lu(k, ii) = lu(k, ii) - subt * lu(k, icol)
    END DO
  END DO
  
!     Find maximum entry
  
  60 pivot = acu
  kk = 0
  DO i = ii, m
    k = indx(i)
    IF (ABS(lu(k, ii)) <= pivot) CYCLE
    pivot = ABS(lu(k, ii))
    kk = i
  END DO
  IF (kk == 0) GO TO 10
  
!     Switch order
  
  isave = indx(kk)
  indx(kk) = indx(ii)
  indx(ii) = isave
  
!     Put into columns of LU one at a time
  
  IF (intl) ibase(ii) = irow
  IF (ii == m) CYCLE
  j = ii + 1
  DO i = j, m
    k = indx(i)
    lu(k, ii) = lu(k, ii) / lu(isave, ii)
  END DO
END DO
kkk = irow

RETURN
END SUBROUTINE update

END SUBROUTINE lfnorm



PROGRAM test_lfnorm
IMPLICIT NONE

INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER                :: i, ier, niter, m, n
REAL (dp), ALLOCATABLE :: x(:,:), y(:), beta(:)
REAL (dp)              :: z

INTERFACE
  SUBROUTINE lfnorm(n, m, x, y, beta, z, ky, ifault)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n
    INTEGER, INTENT(IN)     :: m
    REAL (dp), INTENT(IN)   :: x(:, :)   ! x(ndim, mdim)
    REAL (dp), INTENT(IN)   :: y(:)
    REAL (dp), INTENT(OUT)  :: beta(:)
    REAL (dp), INTENT(OUT)  :: z
    INTEGER, INTENT(OUT)    :: ky
    INTEGER, INTENT(OUT)    :: ifault
  END SUBROUTINE lfnorm
END INTERFACE

DO
  WRITE(*, '(a)', ADVANCE='NO') ' Enter number of predictor variables: '
  READ(*, *) m
  IF (m < 1) THEN
    WRITE(*, *) '** Must be greater than zero **'
    CYCLE
  ELSE
    EXIT
  END IF
END DO

n = 3 * m

ALLOCATE( x(n,m), y(n), beta(m) )

! Generate artificial data with coefficients 1, 2, ..., m.

beta = (/ (DBLE(i),i=1,m) /)
DO i = 1, n
  CALL RANDOM_NUMBER( x(i,2:) )
  CALL RANDOM_NUMBER( z )
  x(i,1) = 1.0_dp
  y(i) = DOT_PRODUCT( x(i,:), beta ) + z - 0.5_dp
END DO

CALL lfnorm(n, m, x, y, beta, z, niter, ier)
IF (ier > 0) THEN
  WRITE(*, *) '** X-matrix has rank < m **'
ELSE
  WRITE(*, *) 'Regression coefficients'
  WRITE(*, '(t2, 10f7.3)') beta
  WRITE(*, '(a, i6)') ' No. of iterations = ', niter
  WRITE(*, '(a, f9.5)') ' Least max. abs. deviation = ', z
END IF

STOP
END PROGRAM test_lfnorm
