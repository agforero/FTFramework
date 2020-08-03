MODULE find_subsets

! A module of routines for finding and recording best-fitting subsets of
! regression variables

!   Version 1.10, 17 February 2004
!   Author: Alan Miller
!           formerly of CSIRO Division of Mathematical & Information Sciences
!   Phone: (+61) 3 9592-5085
!   e-mail: amiller @ bigpond.net.au
!   WWW-page: http://users.bigpond.net.au/amiller/

! 17 Feb 2004 Correction to subroutine EFROYM for the case in which all the
!             variables are selected. Thanks to David Jones.
! 12 Nov 1999 Made changes to routines exadd1 & add2 to prevent the calculation
!             of negative residual sums of squares which could occur in cases
!             in which the true RSS is zero.   Routine seq2 changed to avoid
!             cycling.
! 24 May 2000 Changed lsq_kind to dp (cosmetic change only)
! 4 June 2000 Added routine random_pick which picks a random set of variables.
! 29 Aug 2002 Set value of size in subroutine EFROYM when IER /= 0.

USE lsq

IMPLICIT NONE

INTEGER, SAVE                 :: max_size, nbest, lopt_dim1
REAL (dp), ALLOCATABLE, SAVE  :: bound(:), ress(:,:)
INTEGER, ALLOCATABLE, SAVE    :: lopt(:,:)


CONTAINS

SUBROUTINE init_subsets(nvar_max, fit_const, nvar)

INTEGER, INTENT(IN)           :: nvar_max
LOGICAL, INTENT(IN)           :: fit_const
INTEGER, INTENT(IN), OPTIONAL :: nvar

!     Local variables

INTEGER    :: i, ier
REAL (dp)  :: eps = 1.E-14
LOGICAL    :: lindep(ncol)

!     The LSQ module has probably already been initialized, but just in case ..

IF (.NOT. initialized) THEN
  IF (PRESENT(nvar)) CALL startup(nvar, fit_const)
END IF

IF (fit_const) THEN
  max_size = nvar_max + 1
ELSE
  max_size = nvar_max
END IF

lopt_dim1 = max_size * (max_size + 1) / 2
IF (ALLOCATED(bound)) DEALLOCATE(bound, ress, lopt)
ALLOCATE (bound(max_size), ress(max_size,nbest), lopt(lopt_dim1,nbest))

bound = HUGE(eps)
ress  = HUGE(eps)
lopt  = 0

CALL tolset(eps)
CALL sing(lindep, ier)

CALL ss()
DO i = 1, max_size
  CALL report(i, rss(i))
END DO

RETURN
END SUBROUTINE init_subsets


SUBROUTINE add1(first, last, ss, smax, jmax, ier)

! Calculate the reduction in residual sum of squares when one variable,
! selected from those in positions FIRST .. LAST, is added in position FIRST,
! given that the variables in positions 1 .. FIRST-1 (if any) are already
! included.

INTEGER, INTENT(IN)     :: first, last
INTEGER, INTENT(OUT)    :: jmax, ier
REAL (dp), INTENT(OUT)  :: ss(:), smax

!     Local variables

INTEGER    :: j, inc, pos, row, col
REAL (dp)  :: zero = 0.0_dp, diag, dy, ssqx, sxx(ncol), sxy(ncol)

!     Check call arguments

jmax = 0
smax = zero
ier = 0
IF (first > ncol) ier = 1
IF (last < first) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

!     Accumulate sums of squares & products from row FIRST

sxx(first:last) = zero
sxy(first:last) = zero
inc = ncol - last
pos = row_ptr(first)
DO row = first, last
  diag = d(row)
  dy = diag * rhs(row)
  sxx(row) = sxx(row) + diag
  sxy(row) = sxy(row) + dy
  DO col = row+1, last
    sxx(col) = sxx(col) + diag * r(pos)**2
    sxy(col) = sxy(col) + dy * r(pos)
    pos = pos + 1
  END DO
  pos = pos + inc
END DO

!     Incremental sum of squares for a variable = SXY * SXY / SXX.
!     Calculate whenever sqrt(SXX) > TOL for that variable.

DO j = first, last
  ssqx = sxx(j)
  IF (SQRT(ssqx) > tol(j)) THEN
    ss(j) = sxy(j)**2 / sxx(j)
    IF (ss(j) > smax) THEN
      smax = ss(j)
      jmax = j
    END IF
  ELSE
    ss(j) = zero
  END IF
END DO

RETURN
END SUBROUTINE add1


SUBROUTINE add2(first, last, smax, j1, j2, ier)

!     Calculate the maximum reduction in residual sum of squares when 2
!     variables, selected from those in positions FIRST .. LAST, are
!     added, given that the variables in positions 1 .. FIRST-1 (if
!     any) are already included.    J1, J2 are the positions of the two
!     best variables.   N.B. J2 < J1.

INTEGER, INTENT(IN)     :: first, last
INTEGER, INTENT(OUT)    :: j1, j2, ier
REAL (dp), INTENT(OUT)  :: smax

!     Local variables

INTEGER    :: start, i1, i2, row, pos1, pos2, inc
REAL (dp)  :: zero = 0.0_dp, temp, det, two = 2.0, sxx(ncol), sxy(ncol), sx1x2

!     Check call arguments

smax = zero
j1 = 0
j2 = 0
ier = 0
IF (first > ncol) ier = 1
IF (last <= first) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

start = row_ptr(first)

!     Cycle through all pairs of variables from those between FIRST & LAST.

DO i1 = first, last
  sxx(i1) = d(i1)
  sxy(i1) = d(i1) * rhs(i1)
  pos1 = start + i1 - first - 1
  DO row = first, i1-1
    temp = d(row) * r(pos1)
    sxx(i1) = sxx(i1) + temp*r(pos1)
    sxy(i1) = sxy(i1) + temp*rhs(row)
    pos1 = pos1 + ncol - row - 1
  END DO

  DO i2 = first, i1-1
    pos1 = start + i1 - first - 1
    pos2 = start + i2 - first - 1
    sx1x2 = zero
    DO row = first, i2-1
      sx1x2 = sx1x2 + d(row)*r(pos1)*r(pos2)
      inc = ncol - row - 1
      pos1 = pos1 + inc
      pos2 = pos2 + inc
    END DO
    sx1x2 = sx1x2 + d(i2)*r(pos1)

!     Calculate reduction in RSS for pair I1, I2.
!     The sum of squares & cross-products are in:
!              ( SXX(I1)  SX1X2   )      ( SXY(I1) )
!              ( SX1X2    SXX(I2) )      ( SXY(I2) )

    det = MAX( (sxx(i1) * sxx(i2) - sx1x2**2), zero)
    temp = SQRT(det)
    IF (temp < tol(i1)*SQRT(sxx(i2)) .OR.             &
        temp < tol(i2)*SQRT(sxx(i1))) CYCLE
    temp = ((sxx(i2)*sxy(i1) - two*sx1x2*sxy(i2))*sxy(i1) + sxx(i1)*sxy(i2)**2) &
           / det
    IF (temp > smax) THEN
      smax = temp
      j1 = i1
      j2 = i2
    END IF
  END DO ! i2 = first, i1-1
END DO   ! i1 = first, last

RETURN
END SUBROUTINE add2


SUBROUTINE bakwrd(first, last, ier)

!     Backward elimination from variables in positions FIRST .. LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.
!     On exit, the array VORDER contains the numbers of the variables
!     in the order in which they were deleted.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: pos, jmin, i
REAL (dp)  :: ss(last), smin

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

!     For POS = LAST, ..., FIRST+1 call DROP1 to find best variable to
!     find which variable to drop next.

DO pos = last, first+1, -1
  CALL drop1(first, pos, ss, smin, jmin, ier)
  CALL exdrop1(first, pos, ss, smin, jmin)
  IF (jmin > 0 .AND. jmin < pos) THEN
    CALL vmove(jmin, pos, ier)
    IF (nbest > 0) THEN
      DO i = jmin, pos-1
        CALL report(i, rss(i))
      END DO
    END IF
  END IF
END DO

RETURN
END SUBROUTINE bakwrd


SUBROUTINE drop1(first, last, ss, smin, jmin, ier)

! Calculate the increase in the residual sum of squares when the variable in
! position J is dropped from the model (i.e. moved to position LAST),
! for J = FIRST, ..., LAST-1.

INTEGER, INTENT(IN)     :: first, last
INTEGER, INTENT(OUT)    :: jmin, ier
REAL (dp), INTENT(OUT)  :: ss(:), smin

!     Local variables

INTEGER    :: j, pos1, inc, pos, row, col, i
REAL (dp)  :: large = HUGE(1.0_dp), zero = 0.0_dp, d1, rhs1, d2, x, wk(last), &
              vsmall = TINY(1.0_dp)

!     Check call arguments

jmin = 0
smin = large
ier = 0
IF (first > ncol) ier = 1
IF (last < first) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

!     POS1 = position of first element of row FIRST in r.

pos1 = row_ptr(first)
inc = ncol - last

!     Start of outer cycle for the variable to be dropped.

DO j = first, last
  d1 = d(j)
  IF (SQRT(d1) < tol(j)) THEN
    ss(j) = zero
    smin = zero
    jmin = j
    GO TO 50
  END IF
  rhs1 = rhs(j)
  IF (j == last) GO TO 40

!     Copy row J of R into WK.

  pos = pos1
  DO i = j+1, last
    wk(i) = r(pos)
    pos = pos + 1
  END DO
  pos = pos + inc

!     Lower the variable past each row.

  DO row = j+1, last
    x = wk(row)
    d2 = d(row)
    IF (ABS(x) * SQRT(d1) < tol(row) .OR. d2 < vsmall) THEN
      pos = pos + ncol - row
      CYCLE
    END IF
    d1 = d1 * d2 / (d2 + d1 * x**2)
    DO col = row+1, last
      wk(col) = wk(col) - x * r(pos)
      pos = pos + 1
    END DO
    rhs1 = rhs1 - x * rhs(row)
    pos = pos + inc
  END DO
  40 ss(j) = rhs1 * d1 * rhs1
  IF (ss(j) < smin) THEN
    jmin = j
    smin = ss(j)
  END IF

!     Update position of first element in row of r.

  50 IF (j < last) pos1 = pos1 + ncol - j
END DO

RETURN
END SUBROUTINE drop1


SUBROUTINE efroym(first, last, fin, fout, size, ier, lout)

!     Efroymson's stepwise regression from variables in positions FIRST,
!     ..., LAST.  If FIRST > 1, variables in positions prior to this are
!     forced in.  If LAST < ncol, variables in positions after this are
!     forced out.

!     A report is written to unit LOUT if LOUT >= 0.

INTEGER, INTENT(IN)    :: first, last, lout
INTEGER, INTENT(OUT)   :: size, ier
REAL (dp), INTENT(IN)  :: fin, fout

!     Local variables

INTEGER    :: jmax, jmin, i
REAL (dp)  :: one = 1.0, eps, zero = 0.0, ss(last), smax, base, var, f, smin

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (fin < fout .OR. fin <= zero) ier = ier + 256
IF (nobs <= ncol) ier = ier + 512
IF (ier /= 0) THEN
  size = 0
  RETURN
END IF

!     EPS approximates the smallest quantity such that the calculated value of
!     (1 + EPS) is > 1.   It is used to test for a perfect fit (RSS = 0).

eps = EPSILON(one)

!     SIZE = number of variables in the current subset

size = first - 1

!     Find the best variable to add next

20 CALL add1(size+1, last, ss, smax, jmax, ier)
IF (nbest > 0) CALL exadd1(size+1, smax, jmax, ss, last)

!     Calculate 'F-to-enter' value

IF (size > 0) THEN
  base = rss(size)
ELSE
  base = rss(1) + ss(1)
END IF
var = (base - smax) / (nobs - size - 1)
IF (var < eps*base) THEN
  ier = -1
  f = zero
ELSE
  f = smax / var
END IF
IF (lout >= 0) WRITE(lout, 900) vorder(jmax), f
900 FORMAT(' Best variable to add:  ', i4, '  F-to-enter = ', f10.2)

!     Exit if F < FIN or IER < 0 (perfect fit)

IF (f < fin .OR. ier < 0) RETURN

!     Add the variable to the subset (in position FIRST).

IF (lout >= 0) WRITE(lout, '(50x, "Variable added")')
size = size + 1
IF (jmax > first) CALL vmove(jmax, first, ier)
DO i = first, MIN(jmax-1, max_size)
  CALL report(i, rss(i))
END DO

!     See whether a variable entered earlier can be deleted now.

30 IF (size <= first) GO TO 20
CALL drop1(first+1, size, ss, smin, jmin, ier)
CALL exdrop1(first+1, size, ss, smin, jmin)
var = rss(size) / (nobs - size)
f = smin / var
IF (lout >= 0) WRITE(lout, 910) vorder(jmin), f
910 FORMAT(' Best variable to drop: ', i4, '  F-to-drop  = ', f10.2)

IF (f < fout) THEN
  IF (lout >= 0) WRITE(lout, '(50x, "Variable dropped")')
  CALL vmove(jmin, size, ier)
  IF (nbest > 0) THEN
    DO i = jmin, size-1
      CALL report(i, rss(i))
    END DO
  END IF
  size = size - 1
  GO TO 30
END IF

IF (size >= last) RETURN
GO TO 20
END SUBROUTINE efroym


SUBROUTINE exadd1(ivar, smax, jmax, ss, last)

!     Update the NBEST subsets of IVAR variables found from a call
!     to subroutine ADD1.

INTEGER, INTENT(IN)    :: ivar, jmax, last
REAL (dp), INTENT(IN)  :: smax, ss(:)

!     Local variables

REAL (dp)  :: zero = 0.0_dp, ssbase, sm, temp, wk(last)
INTEGER    :: i, j, ltemp, jm

IF (jmax == 0) RETURN
IF (ivar <= 0) RETURN
IF (ivar > max_size) RETURN
ltemp = vorder(ivar)
jm = jmax
sm = smax
IF (ivar > 1) ssbase = rss(ivar-1)
IF (ivar == 1) ssbase = rss(ivar) + ss(1)
wk(ivar:last) = ss(ivar:last)

DO i = 1, nbest
  temp = MAX(ssbase - sm, zero)
  IF (temp >= bound(ivar)) EXIT
  vorder(ivar) = vorder(jm)
  IF (jm == ivar) vorder(ivar) = ltemp
  CALL report(ivar, temp)
  IF (i >= nbest) EXIT
  wk(jm) = zero
  sm = zero
  jm = 0
  DO j = ivar, last
    IF (wk(j) <= sm) CYCLE
    jm = j
    sm = wk(j)
  END DO
  IF (jm == 0) EXIT
END DO

!     Restore VORDER(IVAR)

vorder(ivar) = ltemp

RETURN
END SUBROUTINE exadd1


SUBROUTINE exdrop1(first, last, ss, smin, jmin)
! Record any new subsets of (LAST-1) variables found from a call to DROP1

INTEGER, INTENT(IN)    :: first, last, jmin
REAL (dp), INTENT(IN)  :: ss(:), smin

! Local variables
INTEGER    :: list(1:last), i
REAL (dp)  :: rss_last, ssq

IF (jmin == 0 .OR. last < 1 .OR. last-1 > max_size) RETURN

rss_last = rss(last)
IF (rss_last + smin > bound(last-1)) RETURN

list = vorder(1:last)
DO i = first, last-1
  vorder(i:last-1) = list(i+1:last)
  ssq = rss_last + ss(i)
  CALL report(last-1, ssq)
  vorder(i) = list(i)
END DO

RETURN
END SUBROUTINE exdrop1


SUBROUTINE forwrd(first, last, ier)

!     Forward selection from variables in positions FIRST .. LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.
!     On exit, the array VORDER contains the numbers of the variables
!     in the order in which they were added.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: pos, jmax
REAL (dp)  :: ss(last), smax

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0) RETURN

!     For POS = FIRST .. max_size, call ADD1 to find best variable to put
!     into position POS.

DO pos = first, max_size
  CALL add1(pos, last, ss, smax, jmax, ier)
  IF (nbest > 0) CALL exadd1(pos, smax, jmax, ss, last)

!     Move the best variable to position POS.

  IF (jmax > pos) CALL vmove(jmax, pos, ier)
END DO

RETURN
END SUBROUTINE forwrd


SUBROUTINE report(nv, ssq)

!     Update record of the best NBEST subsets of NV variables, if
!     necessary, using SSQ.

INTEGER, INTENT(IN)    :: nv
REAL (dp), INTENT(IN)  :: ssq

!     Local variables

INTEGER    :: rank, pos1, j, list(nv)
REAL (dp)  :: under1 = 0.99999_dp, above1 = 1.00001_dp

!     If residual sum of squares (SSQ) for the new subset > the
!     appropriate bound, return.

IF(nv > max_size) RETURN
IF(ssq >= bound(nv)) RETURN
pos1 = (nv*(nv-1))/2 + 1

!     Find rank of the new subset

DO rank = 1, nbest
  IF(ssq < ress(nv,rank)*above1) THEN
    list = vorder(1:nv)
    CALL shell(list, nv)

!     Check list of variables if ssq is almost equal to ress(nv,rank) -
!     to avoid including the same subset twice.

    IF (ssq > ress(nv,rank)*under1) THEN
      IF (same_vars(list, lopt(pos1:,rank), nv)) RETURN
    END IF

!     Record the new subset, and move the others down one place.

    DO j = nbest-1, rank, -1
      ress(nv,j+1) = ress(nv,j)
      lopt(pos1:pos1+nv-1, j+1) = lopt(pos1:pos1+nv-1, j)
    END DO
    ress(nv,rank) = ssq
    lopt(pos1:pos1+nv-1, rank) = list(1:nv)
    bound(nv) = ress(nv,nbest)
    RETURN
  END IF
END DO

RETURN
END SUBROUTINE report



SUBROUTINE shell(l, n)

!      Perform a SHELL-sort on integer array L, sorting into increasing order.

!      Latest revision - 5 July 1995

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN OUT) :: l(:)

!     Local variables
INTEGER   :: start, finish, temp, new, i1, i2, incr, it

incr = n
DO
  incr = incr/3
  IF (incr == 2*(incr/2)) incr = incr + 1
  DO start = 1, incr
    finish = n

!      TEMP contains the element being compared; IT holds its current
!      location.   It is compared with the elements in locations
!      IT+INCR, IT+2.INCR, ... until a larger element is found.   All
!      smaller elements move INCR locations towards the start.   After
!      each time through the sequence, the FINISH is decreased by INCR
!      until FINISH <= INCR.

    20 i1 = start
    temp = l(i1)
    it = i1

!      I2 = location of element new to be compared with TEMP.
!      Test I2 <= FINISH.

    DO
      i2 = i1 + incr
      IF (i2 > finish) THEN
        IF (i1 > it) l(i1) = temp
        finish = finish - incr
        EXIT
      END IF
      new = l(i2)

!     If TEMP > NEW, move NEW to lower-numbered position.

      IF (temp > new) THEN
        l(i1) = new
        i1 = i2
        CYCLE
      END IF

!     TEMP <= NEW so do not swap.
!     Use NEW as the next TEMP.

      IF (i1 > it) l(i1) = temp
      i1 = i2
      temp = new
      it = i1

!     Repeat until FINISH <= INCR.
    END DO

    IF (finish > incr) GO TO 20
  END DO

!      Repeat until INCR = 1.

  IF (incr <= 1) RETURN
END DO

RETURN
END SUBROUTINE shell



FUNCTION same_vars(list1, list2, n) RESULT(same)

LOGICAL              :: same
INTEGER, INTENT(IN)  :: n, list1(:), list2(:)

same = ALL(list1(1:n) == list2(1:n))

RETURN
END FUNCTION same_vars



SUBROUTINE seq2(first, last, ier)

! Sequential replacement algorithm applied to the variables in positions
! FIRST, ..., LAST.   2 variables at a time are added or replaced.
! If FIRST > 1, variables in positions prior to this are forced in.
! If LAST < NP, variables in positions after this are left out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER  :: nv, nsize

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

nv = MIN(max_size, last-1)

!     Outer loop; SIZE = current size of subset being considered.

DO nsize = first+1, nv
  CALL replace2(first, last, nsize)
END DO

RETURN
END SUBROUTINE seq2



SUBROUTINE replace2(first, last, nsize)
! Replace 2 variables at a time from those in positions first, ..., nsize
! with 2 from positions nsize, .., last - if they reduce the RSS.

INTEGER, INTENT(IN)  :: first, last, nsize

! Local variables

INTEGER              :: ier, j1, j2, pos1, pos2, best(2), i, iwk(last)
REAL (dp)            :: smax, rssnew, rssmin, save_rss
REAL (dp), PARAMETER :: zero = 0.0_dp

10 best(1) = 0
best(2) = 0
rssmin = rss(nsize)

!     Two loops to place all pairs of variables in positions nsize-1 and nsize.
!     POS1 = destination for variable from position nsize.
!     POS2 = destination for variable from position nsize-1.

DO pos1 = first, nsize
  DO pos2 = pos1, nsize-1
    CALL add2(nsize-1, last, smax, j1, j2, ier)

    IF (j1+j2 > nsize + nsize - 1) THEN
      rssnew = MAX(rss(nsize-2) - smax, zero)
      IF (rssnew < rssmin) THEN
        best(1) = vorder(j1)
        best(2) = vorder(j2)
        iwk(1:nsize-2) = vorder(1:nsize-2)
        rssmin = rssnew
      END IF
    END IF

    CALL vmove(nsize-1, pos2, ier)
  END DO
  CALL vmove(nsize, pos1, ier)
  DO i = pos1, nsize
    CALL report(i, rss(i))
  END DO
END DO

!     If any replacement reduces the RSS, make the best one.

IF (best(1) + best(2) > 0) THEN
  iwk(nsize-1) = best(2)
  iwk(nsize) = best(1)
  save_rss = rss(nsize)
  CALL reordr(iwk, nsize, 1, ier)
  DO i = first, nsize
    CALL report(i, rss(i))
  END DO

!    The calculated value of rssmin above is only a rough approximation to
!    the real residual sum of squares, thiugh usually good enough.
!    The new value of rss(nsize) is more accurate.   It is used below
!    to avoid cycling when several subsets give the same RSS.

  IF (rss(nsize) < save_rss) GO TO 10
END IF

RETURN
END SUBROUTINE replace2



SUBROUTINE seqrep(first, last, ier)

!     Sequential replacement algorithm applied to the variables in
!     positions FIRST, ..., LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: nv, size, start, best, from, i, jmax, count, j
REAL (dp)  :: zero = 0.0_dp, ssred, ss(last), smax

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

nv = MIN(max_size, last-1)

!     Outer loop; SIZE = current size of subset being considered.

DO size = first, nv
  count = 0
  start = first
  10 ssred = zero
  best = 0
  from = 0

!     Find the best variable from those in positions SIZE+1, ..., LAST
!     to replace the one in position SIZE.   Then rotate variables in
!     positions START, ..., SIZE.

  DO i = start, size
    CALL add1(size, last, ss, smax, jmax, ier)
    IF (jmax > size) THEN
      CALL exadd1(size, smax, jmax, ss, last)
      IF (smax > ssred) THEN
        ssred = smax
        best = jmax
        IF (i < size) THEN
          from = size + start - i - 1
        ELSE
          from = size
        END IF
      END IF
    END IF
    IF (i < size) CALL vmove(size, start, ier)
    DO j = start, size-1
      CALL report(j, rss(j))
    END DO
  END DO ! i = start, size

!     If any replacement reduces the RSS, make the best one.
!     Move variable from position FROM to SIZE.
!     Move variable from position BEST to FIRST.

  IF (best > size) THEN
    IF (from < size) CALL vmove(from, size, ier)
    CALL vmove(best, first, ier)
    DO j = first, best-1
      CALL report(j, rss(j))
    END DO
    count = 0
    start = first + 1
  ELSE
    count = count + 1
  END IF

!     Repeat until COUNT = SIZE - START + 1

  IF (count <= size - start) GO TO 10
END DO

RETURN
END SUBROUTINE seqrep


SUBROUTINE xhaust(first, last, ier)

!     Exhaustive search algorithm, using leaps and bounds, applied to
!     the variables in positions FIRST, ..., LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.

INTEGER, INTENT(IN)  :: first, last
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER    :: row, i, jmax, ipt, newpos, iwk(max_size)
REAL (dp)  :: ss(last), smax, temp

!     Check call arguments

ier = 0
IF (first >= ncol) ier = 1
IF (last <= 1) ier = ier + 2
IF (first < 1) ier = ier + 4
IF (last > ncol) ier = ier + 8
IF (ier /= 0 .OR. nbest <= 0) RETURN

!     Record subsets contained in the initial ordering, including check
!     for variables which are linearly related to earlier variables.
!     This should be redundant if the user has first called SING and
!     init_subsets.

DO row = first, max_size
  IF (d(row) <= tol(row)) THEN
    ier = -999
    RETURN
  END IF
  CALL report(row, rss(row))
END DO

!     IWK(I) contains the upper limit for the I-th simulated DO-loop for
!     I = FIRST, ..., max_size-1.
!     IPT points to the current DO loop.

iwk(first:max_size) = last

!     Innermost loop.
!     Find best possible variable for position max_size from those in
!     positions max_size, .., IWK(max_size).

30 CALL add1(max_size, iwk(max_size), ss, smax, jmax, ier)
CALL exadd1(max_size, smax, jmax, ss, iwk(max_size))

!     Move to next lower numbered loop which has not been exhausted.

ipt = max_size - 1
40 IF (ipt >= iwk(ipt)) THEN
  ipt = ipt - 1
  IF (ipt >= first) GO TO 40
  RETURN
END IF

!     Lower variable from position IPT to position IWK(IPT).
!     Record any good new subsets found by the move.

newpos = iwk(ipt)
CALL vmove(ipt, newpos, ier)
DO i = ipt, MIN(max_size, newpos-1)
  CALL report(i, rss(i))
END DO

!     Reset all ends of loops for I >= IPT.

iwk(ipt:max_size) = newpos - 1

!     If residual sum of squares for all variables above position NEWPOS
!     is greater than BOUND(I), no better subsets of size I can be found
!     inside the current loop.

temp = rss(newpos-1)
DO i = ipt, max_size
  IF (temp > bound(i)) GO TO 80
END DO
IF (iwk(max_size) > max_size) GO TO 30
ipt = max_size - 1
GO TO 40

80 ipt = i - 1
IF (ipt < first) RETURN
GO TO 40

END SUBROUTINE xhaust



SUBROUTINE random_pick(first, last, npick)
! Pick npick variables at random from those in positions first, ..., last
! and move them to occupy positions starting from first.

INTEGER, INTENT(IN)  :: first, last, npick

! Local variables

INTEGER  :: first2, i, ilist(1:last), j, k, navail
REAL     :: r

navail = last + 1 - first
IF (npick >= navail .OR. npick <= 0) RETURN
DO i = first, last
  ilist(i) = vorder(i)
END DO

first2 = first
DO i = 1, npick
  CALL RANDOM_NUMBER(r)
  k = first2 + r*navail
  IF (k > first2) THEN
    j = ilist(first2)
    ilist(first2) = ilist(k)
    ilist(k) = j
  END IF
  first2 = first2 + 1
  navail = navail - 1
END DO

CALL reordr(ilist(first:), npick, first, i)

RETURN
END SUBROUTINE random_pick

END MODULE find_subsets
