MODULE d_optimal_design
IMPLICIT NONE

! Latest revision - 4 March 1998

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

CONTAINS

SUBROUTINE dopt(x, dim1, ncand, kin, n, nblock, in, blksiz, k, rstart, nrbar, &
                d, rbar, picked, lndet, xx, tol, zpz, wk, ifault)

!     Heuristic algorithm to pick N rows of X out of NCAND to maximize
!     the determinant of X'X, using the Fedorov exchange algorithm.

!     Code written by Alan Miller and Nam Nguyen from CSIRO, Melbourne,
!     Australia and published in Applied Statistics as AS 295.

!     This version differs from the official version in the following
!     ways:
!     1. When an existing design is augmented, the design points of the
!        existing design are NOT used as candidates for augmentation.
!     2. Some Fortran 90 features are used.

INTEGER, INTENT(IN)       :: dim1, ncand, kin, n, nblock, in(:),       &
                             blksiz(:), k, nrbar
INTEGER, INTENT(OUT)      :: picked(:), ifault
REAL, INTENT(IN)          :: x(:,:)
REAL (dp), INTENT(IN OUT) :: d(:), rbar(:), lndet, xx(:), tol(:), zpz(:,:), &
                             wk(:)
LOGICAL, INTENT(IN)       :: rstart

!     Local variables

INTEGER         :: i, j, nin, point, case, nb, block, l, pos, best, first,  &
                   last, cand, lastin, lstout, drop, rempos, bl, rpos,      &
                   last1, last2, first1, first2, pos1, pos2, block1,        &
                   block2, case1, case2, posi, posj, bestb1, bestb2,        &
                   bestp1, bestp2, rank, mxrank, inc, navail
REAL (dp)       :: temp, dmax, total
REAL            :: rand
LOGICAL         :: change
REAL (dp), PARAMETER :: one = 1.D0, minus1 = -1.D0, eps = EPSILON(one),   &
                        above1 = 1.0001D0, small = 1.D-04, hundrd = 100.D0
! REAL (dp)  :: det

!     Argument checks

ifault = 0
IF (dim1 < ncand) ifault = 1
IF (k > n) ifault = ifault + 2
IF (nrbar < k* (k-1)/2) ifault = ifault + 4
IF (k /= kin+nblock) ifault = ifault + 8
IF (nblock > 1) THEN
  l = 0
  DO block = 1, nblock
    l = l + blksiz(block)
  END DO
  IF (n /= l) ifault = ifault + 16
ELSE
  IF (n /= blksiz(1)) ifault = ifault + 16
END IF

!     NB = max(1, Nblock) so that we can force it to go through DO-loops
!     once.   NIN = no. of design points forced into the design.

nb = MAX(1, nblock)
nin = 0
DO i = 1, nb
  IF (in(i) < 0) GO TO 30
  IF (in(i) > 0) THEN
    IF (in(i) > blksiz(i)) GO TO 30
    nin = nin + in(i)
  END IF
END DO
navail = ncand - nin
IF (nin <= n) GO TO 40
30 ifault = ifault + 32

40 IF (ifault /= 0) RETURN
CALL clear(k, nrbar, d, rbar, ifault)

!     Set up an array of tolerances.

tol(1:k) = eps
block = 1
DO case = 1, ncand
  CALL getx(x, kin, nblock, block, xx, case)
  DO i = 1, k
    tol(i) = tol(i) + ABS(xx(i))
  END DO
END DO
temp = REAL(n) * eps / ncand
DO i = 1, k
  IF (i <= nblock) THEN
    tol(i) = eps
  ELSE
    tol(i) = tol(i) * temp
  END IF
END DO

!     Form initial Cholesky factorization.

pos = 1
DO block = 1, nb
  IF (rstart) THEN
    last1 = (in(block) + blksiz(block))/2
    inc = SQRT(REAL(ncand, KIND=dp) + small)
  END IF
  DO i = 1, blksiz(block)
    IF (rstart .AND. i > in(block)) THEN
      CALL RANDOM_NUMBER(rand)
      point = nin + 1 + navail * rand

!     If I <= LAST1, use a random point, otherwise find the candidate which
!     maximizes the rank, and then maximizes the subspace determinant for
!     that rank.

      IF (i > last1) THEN
        mxrank = 0
        lndet = -hundrd
        DO cand = nin+1, ncand, inc
          CALL getx(x, kin, nblock, block, xx, point)
          CALL modtr2(k, xx, d, rbar, tol, rank, total)
          IF (rank < mxrank) GO TO 90
          IF (rank == mxrank .AND. total < lndet) GO TO 90
          best = point
          mxrank = rank
          lndet = total * above1
          90 point = point + inc
          IF (point > ncand) point = point - navail
        END DO
        point = best
      END IF
      picked(pos) = point
    ELSE

!     Case in which a full design has been input, or points are to be
!     forced into the design.

      point = picked(pos)
    END IF

!     Augment the Cholesky factorization.

    CALL getx(x, kin, nblock, block, xx, point)
    CALL modtri(k, one, xx, d, rbar, tol)
    pos = pos + 1
  END DO
END DO

!     Adjust factorization in case of singular matrix.

CALL singm(k, nrbar, d, rbar, tol, wk, ifault)

!     If rank of input design < K, try replacing points.

IF (ifault == 0) GO TO 280

!     Find first row of Cholesky factorization with a zero multiplier.

180 DO pos = 1, k
  IF (d(pos) < tol(pos)) GO TO 200
END DO
GO TO 280

!     Find linear relationship between variable in position POS and the
!     previous variables.

200 l = pos - 1
DO i = 1, pos-1
  wk(i) = rbar(l)
  l = l + k - i - 1
END DO
CALL regcf(k, nrbar, d, rbar, wk, tol, wk, pos-1, ifault)

!     Find a candidate point which does not satisfy this linear
!     relationship.   Use a random start.

bl = 1
CALL RANDOM_NUMBER(rand)
case = nin + 1 + navail * rand
DO cand = 1, ncand
  CALL getx(x, kin, nblock, bl, xx, case)
  total = xx(pos) - DOT_PRODUCT( wk(1:pos-1), xx(1:pos-1) )
  IF (ABS(total) > hundrd * tol(pos)) GO TO 240
  case = case + 1
  IF (case > ncand) case = nin + 1
END DO

!     Failed to find any candidate point which would make the design
!     of higher rank.

ifault = -1
RETURN

!     Before adding the point, find one which it can replace without
!     lowering the rank.

240 bl = 0
temp = one - small
pos = in(1) + 1
DO block = 1, nb
  DO j = in(block)+1, blksiz(block)
    l = picked(pos)
    CALL getx(x, kin, nblock, block, xx, l)
    CALL bksub2(rbar, k, xx, wk)
    total = DOT_PRODUCT( wk(1:k), wk(1:k)/d(1:k) )
    IF (total < temp) THEN
      temp = total
      rempos = pos
      bl = block
    END IF
    pos = pos + 1
  END DO
  IF (block < nblock) pos = pos + in(block+1)
END DO

!     If BL = 0 it means that any point removed from the existing design
!     would reduce the rank.

IF (bl == 0) THEN
  ifault = -1
  RETURN
END IF

!     Add candidate CASE in block BL, then delete the design point
!     already in that position.

CALL getx(x, kin, nblock, bl, xx, case)
CALL modtri(k, one, xx, d, rbar, tol)
l = picked(rempos)
CALL getx(x, kin, nblock, bl, xx, l)
CALL modtri(k, minus1, xx, d, rbar, tol)
picked(rempos) = case
!     write(*, 8000) case, l, bl
!8000 format(' Replacing candidate ', i4, ' with candidate ', i4,
!    +    ' in block', i3)
GO TO 180
!----------------------------------------------------------------------

!     Design is now of full rank.

!     Calculate z'z for all candidate points.
!     z is the solution of R'z = x, so that z'z = x'.inv(X'X).x
!     WK holds sqrt(D) times vector z on return from BKSUB2.

280 DO block = 1, nb
  DO case = nin+1, ncand
    CALL getx(x, kin, nblock, block, xx, case)
    CALL bksub2(rbar, k, xx, wk)
    temp = DOT_PRODUCT( wk(1:k), wk(1:k)/d(1:k) )
    zpz(case, block) = temp
  END DO
END DO
!     det = one
!     do 600 i = 1, k
! 600 det = det * d(i)
!     write(*, 9020) det
!9020 format(' Det = ', g14.6)
!----------------------------------------------------------------------

!     Start of Fedorov exchange algorithm

lastin = 0
lstout = 0
320 change = .false.
last = 0
DO block = 1, nb
  first = last + 1 + in(block)
  last = last + blksiz(block)
  dmax = small
  best = 0

!     Start at a random position within the block.
!     I = no. of point being considered for deletion.

  CALL RANDOM_NUMBER(rand)
  pos = first + (blksiz(block) - in(block)) * rand
  DO case = in(block)+1, blksiz(block)
    pos = pos + 1
    IF (pos > last) pos = first
    i = picked(pos)
    IF (i == lastin) CYCLE
    CALL getx(x, kin, nblock, block, xx, i)
    CALL bksub2(rbar, k, xx, wk)
    CALL bksub1(rbar, k, wk, wk, tol, d)

!     Cycle through the candidates for exchange, using a random start.
!     J = no. of point being considered for addition.

    CALL RANDOM_NUMBER(rand)
    j = nin + 1 + navail * rand
    DO cand = nin+1, ncand
      j = j + 1
      IF (j > ncand) j = nin+1
      IF (j == i .OR. j == lstout) CYCLE
!... The Cauchy-Schwarz test.
      temp = zpz(j,block) - zpz(i,block)
      IF (temp < dmax) CYCLE
!
      CALL getx(x, kin, nblock, block, xx, j)
      total = DOT_PRODUCT( xx(1:k), wk(1:k) )
      temp = temp + total**2 - zpz(i,block)*zpz(j,block)
      IF (temp > dmax) THEN
        dmax = temp * above1
        best = j
        rempos = pos
        drop = i
      END IF
    END DO
  END DO

!     Exchange points BEST & DROP in position REMPOS, if the determinant
!     is increased.

  IF (best /= 0) THEN
    change = .true.
    IF (nb == 1) THEN
      lastin = best
      lstout = drop
    END IF

!     Add the new point, BEST, first to avoid ill-conditioning.
!     Update z'z.

    CALL getx(x, kin, nblock, block, xx, best)
    CALL bksub2(rbar, k, xx, wk)
    CALL bksub1(rbar, k, wk, wk, tol, d)
    CALL modtri(k, one, xx, d, rbar, tol)
    temp = one + zpz(best,block)
    DO bl = 1, nb
      DO case = 1, ncand
        CALL getx(x, kin, nblock, bl, xx, case)
        total = DOT_PRODUCT( xx(1:k), wk(1:k) )
        zpz(case,bl) = zpz(case,bl) - total**2 / temp
      END DO
    END DO

!     Remove the point DROP, and update z'z.

    CALL getx(x, kin, nblock, block, xx, drop)
    CALL bksub2(rbar, k, xx, wk)
    CALL bksub1(rbar, k, wk, wk, tol, d)
    CALL modtri(k, minus1, xx, d, rbar, tol)
    temp = one - zpz(drop,block)
    DO bl = 1, nb
      DO case = 1, ncand
        CALL getx(x, kin, nblock, bl, xx, case)
        total = DOT_PRODUCT( xx(1:k), wk(1:k) )
        zpz(case,bl) = zpz(case,bl) + total**2 / temp
      END DO
    END DO

!         det = one
!         do 620 i = 1, k
! 620     det = det * d(i)
!         write(*, 9010) block, best, drop, dmax, det
!9010     format(' block: ', i3, '  Add : ', i4, '  Drop: ', i4,
!    +           '  dmax: ', f9.3, '  Det: ', g14.6)
    picked(rempos) = best
  END IF
END DO

!     Repeat until there is no further improvement

IF (change) GO TO 320
!---------------------------------------------------------------------

!     If there is more than one block, try swapping treatments between
!     blocks.   This is the Cook & Nachtsheim (1989) algorithm.

IF (nblock <= 1) GO TO 500
!...  RPOS is the position of the first element in RBAR after the rows
!     for the block constants.
rpos = nblock * k - nblock * (nblock + 1)/2 + 1

430 last1 = 0
!...  POS1 & POS2 will hold the positions of the start of the means of
!     the X-variables in the two blocks being considered in RBAR.
pos1 = nblock
dmax = small
change = .false.
DO block1 = 1, nblock-1
  first1 = last1 + 1 + in(block1)
  last1 = last1 + blksiz(block1)
  last2 = last1
  pos2 = pos1 + k - 1 - block1
  DO block2 = block1+1, nblock
    first2 = last2 + 1 + in(block2)
    last2 = last2 + blksiz(block2)
    DO case1 = in(block1)+1, blksiz(block1)
      posi = first1 - 1 + case1
      i = picked(posi)
      CALL getx(x, kin, nblock, block, xx, i)
      zpz(1:kin,1) = xx(nblock+1:nblock+kin)
      DO case2 = in(block2)+1, blksiz(block2)
        posj = first2 - 1 + case2
        j = picked(posj)
        IF (i == j) CYCLE
        CALL getx(x, kin, nblock, block, xx, j)
        zpz(1:kin,2) = xx(nblock+1:nblock+kin)
!...  Pass the orthogonal factorization to DELTA with the top nblock
!     rows removed, i.e. without that part relating to the blocks.
        CALL delta(kin, zpz(:,1), zpz(:,2), rbar(pos1:), rbar(pos2:),   &
             blksiz(block1), blksiz(block2), d(nblock+1:), rbar(rpos:), temp)
        IF (temp > dmax) THEN
          dmax = temp * above1
          bestb1 = block1
          bestb2 = block2
          bestp1 = posi
          bestp2 = posj
          change = .true.
        END IF
      END DO
    END DO
    pos2 = pos2 + k - 1 - block2
  END DO
  pos1 = pos1 + k - 1 - block1
END DO

!     If CHANGE = .TRUE. then make the swap, otherwise the search ends.

IF (change) THEN
  i = picked(bestp1)
  j = picked(bestp2)
  CALL getx(x, kin, nblock, bestb2, xx, i)
  CALL modtri(k, one, xx, d, rbar, tol)
  CALL getx(x, kin, nblock, bestb1, xx, j)
  CALL modtri(k, one, xx, d, rbar, tol)
  CALL getx(x, kin, nblock, bestb1, xx, i)
  CALL modtri(k, minus1, xx, d, rbar, tol)
  CALL getx(x, kin, nblock, bestb2, xx, j)
  CALL modtri(k, minus1, xx, d, rbar, tol)
  picked(bestp1) = j
  picked(bestp2) = i
!       det = one
!       do 610 i = 1, k
! 610   det = det * d(i)
!       write(*, 4000) dmax, det
!4000   format(' DMAX = ', f9.4, 5x, 'New det. = ', g13.5)
  GO TO 430
END IF

!     Calculate log of determinant.

500 lndet = SUM( LOG(d(1:k)) )

RETURN
END SUBROUTINE dopt



SUBROUTINE modtri(np, weight, xrow, d, rbar, tol)

!     ALGORITHM AS75.1 with minor modifications

!     Modify a triangular (Cholesky) decomposition.
!     Calling this routine updates d and rbar by adding another design
!     point with weight = WEIGHT, which may be negative.

!     *** WARNING  Array XROW is overwritten.

INTEGER, INTENT(IN)       :: np
REAL (dp), INTENT(IN)     :: weight, tol(:)
REAL (dp), INTENT(IN OUT) :: xrow(:), d(:), rbar(:)

!     Local variables

INTEGER              :: i, k, nextr
REAL (dp)            :: w, xi, di, wxi, dpi, cbar, sbar, xk
REAL (dp), PARAMETER :: zero = 0.D0

w = weight
nextr = 1
DO i = 1, np

!     Skip unnecessary transformations.   Test on exact zeroes must be
!     used or stability can be destroyed.

  IF (w == zero) RETURN

  xi = xrow(i)
  IF (ABS(xi) < tol(i)) THEN
    nextr = nextr + np - i
    CYCLE
  END IF

  di = d(i)
  wxi = w*xi
  dpi = di + wxi*xi

!     Test for new singularity.

  IF (dpi < tol(i)) THEN
    dpi = zero
    cbar = zero
    sbar = zero
    w = zero
  ELSE
    cbar = di/dpi
    sbar = wxi/dpi
    w = cbar*w
  END IF

  d(i) = dpi
  DO k = i + 1, np
    xk = xrow(k)
    xrow(k) = xk - xi*rbar(nextr)
    rbar(nextr) = cbar*rbar(nextr) + sbar*xk
    nextr = nextr + 1
  END DO
END DO

RETURN
END SUBROUTINE modtri



SUBROUTINE modtr2(np, xrow, d, rbar, tol, rank, lndet)

!     ALGORITHM AS75.1 with modifications.
!     Calculate the effect of an update of a QR factorization upon the
!     rank and determinant, without changing D or RBAR.

!     *** WARNING  Array XROW is overwritten.

INTEGER, INTENT(IN)       :: np
INTEGER, INTENT(OUT)      :: rank
REAL (dp), INTENT(IN)     :: tol(:)
REAL (dp), INTENT(IN OUT) :: xrow(:), d(:), rbar(:)
REAL (dp), INTENT(OUT)    :: lndet

!     Local variables

INTEGER              :: i, k, nextr, j
REAL (dp)            :: w, xi, di, wxi, dpi, cbar, xk
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp

w = one
rank = 0
lndet = zero
nextr = 1
DO i = 1, np

!     Skip unnecessary transformations.   Test on exact zeroes must be
!     used or stability can be destroyed.

  IF (w == zero) THEN
    DO j = i, np
      IF (d(j) > tol(j)) THEN
        rank = rank + 1
        lndet = lndet + LOG(d(j))
      END IF
    END DO
    RETURN
  END IF

  xi = xrow(i)
  IF (ABS(xi) < tol(i)) THEN
    IF (d(i) > tol(i)) THEN
      rank = rank + 1
      lndet = lndet + LOG(d(i))
    END IF
    nextr = nextr + np - i
    CYCLE
  END IF

  di = d(i)
  wxi = w*xi
  dpi = di + wxi*xi

!     Test for new singularity.

  IF (dpi < tol(i)) THEN
    dpi = zero
    cbar = zero
    w = zero
  ELSE
    cbar = di/dpi
    w = cbar*w
    lndet = lndet + LOG(dpi)
    rank = rank + 1
  END IF

  DO k = i + 1, np
    xk = xrow(k)
    xrow(k) = xk - xi*rbar(nextr)
    nextr = nextr + 1
  END DO
END DO

RETURN
END SUBROUTINE modtr2



SUBROUTINE getx(x, kin, nblock, block, xx, case)

!     Copy one case from X to XX.

INTEGER, INTENT(IN)    :: kin, nblock, block, case
REAL (dp), INTENT(OUT) :: xx(:)
REAL, INTENT(IN)       :: x(:,:)

!     Local variables.

INTEGER              :: i
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp

DO i = 1, nblock
  IF (i /= block) THEN
    xx(i) = zero
  ELSE
    xx(i) = one
  END IF
END DO
xx(nblock+1:nblock+kin) = x(case, 1:kin)

RETURN
END SUBROUTINE getx



SUBROUTINE bksub1(rbar, k, rhs, soln, tol, d)

!     Solves  D R y = z for y (SOLN), where z = RHS.
!     RBAR is an upper-triangular matrix with implicit 1's on it's
!     diagonal, stored by rows.

INTEGER, INTENT(IN)    :: k
REAL (dp), INTENT(IN)  :: rbar(:), rhs(:), tol(:), d(:)
REAL (dp), INTENT(OUT) :: soln(:)

!     Local variables

INTEGER              :: row, col, pos
REAL (dp)            :: temp
REAL (dp), PARAMETER :: zero = 0.0_dp

pos = k * (k - 1) / 2
DO row = k, 1, -1
  IF (d(row) > tol(row)) THEN
    temp = rhs(row) / d(row)
    DO col = k, row+1, -1
      temp = temp - rbar(pos) * soln(col)
      pos = pos - 1
    END DO
    soln(row) = temp
  ELSE
    pos = pos - k + row
    soln(row) = zero
  END IF
END DO

RETURN
END SUBROUTINE bksub1



SUBROUTINE singm(np, nrbar, d, rbar, tol, work, ifault)

!     Modified from:
!     ALGORITHM AS274.5  APPL. STATIST. (1992) VOL.41, NO.2

!     Checks for singularities, and adjusts orthogonal reductions produced by
!     AS75.1.

!     Auxiliary routine called: MODTRI

INTEGER, INTENT(IN)       :: np, nrbar
INTEGER, INTENT(OUT)      :: ifault
REAL (dp), INTENT(IN OUT) :: d(:), rbar(:), work(:)
REAL (dp), INTENT(IN)     :: tol(:)

!     Local variables

REAL (dp)            :: temp
REAL (dp), PARAMETER :: zero = 0.0_dp
INTEGER              :: col, pos, row, np2, pos2

!     Check input parameters

ifault = 0
IF (np <= 0) ifault = 1
IF (nrbar < np*(np-1)/2) ifault = ifault + 2
IF (ifault /= 0) RETURN

work(1:np) = SQRT(d(1:np))

DO col = 1, np

!     Set elements within RBAR to zero if they are less than TOL(COL) in
!     absolute value after being scaled by the square root of their row
!     multiplier.

  temp = tol(col)
  pos = col - 1
  DO row = 1, col-1
    IF (ABS(rbar(pos)) * work(row) < temp) rbar(pos) = zero
    pos = pos + np - row - 1
  END DO

!     If diagonal element is near zero, set it to zero, and use MODTRI
!     to augment the projections in the lower rows of the factorization.

  IF (work(col) <= temp) THEN
    ifault = ifault - 1
    IF (col < np) THEN
      np2 = np - col
      pos2 = pos + np - col + 1
      IF (np2 > 1) THEN
        CALL modtri(np2, d(col), rbar(pos+1:), d(col+1:), rbar(pos2:), tol)
      ELSE
        CALL modtri(1, d(col), rbar(pos+1:), d(col+1:), rbar, tol)
      END IF
      rbar(pos+1:pos2-1) = zero
    END IF
    d(col) = zero
  END IF
END DO

RETURN
END SUBROUTINE singm



SUBROUTINE xxtr(np, d, rbar, nreq, trace, rinv)

!     Calculate the trace of the inverse of X'X (= R'R).

INTEGER, INTENT(IN)    :: np, nreq
REAL (dp), INTENT(IN)  :: d(:), rbar(:)
REAL (dp), INTENT(OUT) :: trace, rinv(:)

!     Local variables

INTEGER              :: pos, row, col
REAL (dp), PARAMETER :: one = 1.D0, zero = 0.D0

!     Get the inverse of R

CALL inv(np, rbar, nreq, rinv)

!     Trace = the sum of the diagonal elements of RINV * (1/D) * (RINV)'

trace = zero
pos = 1
DO row = 1, nreq
  trace = trace + one / d(row)
  DO col = row+1, nreq
    trace = trace + rinv(pos)**2 / d(col)
    pos = pos + 1
  END DO
END DO

RETURN
END SUBROUTINE xxtr


SUBROUTINE inv(np, rbar, nreq, rinv)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Invert first NREQ rows and columns of Cholesky factorization
!     produced by AS 75.1.

INTEGER, INTENT(IN)    :: np, nreq
REAL (dp), INTENT(IN)  :: rbar(:)
REAL (dp), INTENT(OUT) :: rinv(:)

!     Local variables.

INTEGER              :: pos, row, col, start, k, pos1, pos2
REAL (dp)            :: sum
REAL (dp), PARAMETER :: zero = 0.0_dp

!     Invert RBAR ignoring row multipliers, from the bottom up.

pos = nreq * (nreq-1)/2
DO row = nreq-1, 1, -1
  start = (row-1) * (np+np-row)/2 + 1
  DO col = nreq, row+1, -1
    pos1 = start
    pos2 = pos
    sum = zero
    DO k = row+1, col-1
      pos2 = pos2 + nreq - k
      sum = sum - rbar(pos1) * rinv(pos2)
      pos1 = pos1 + 1
    END DO
    rinv(pos) = sum - rbar(pos1)
    pos = pos - 1
  END DO
END DO

RETURN
END SUBROUTINE inv



SUBROUTINE clear(np, nrbar, d, rbar, ier)

!     Modified version of:
!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Sets arrays to zero prior to calling INCLUD

INTEGER, INTENT(IN)    :: np, nrbar
INTEGER, INTENT(OUT)   :: ier
REAL (dp), INTENT(OUT) :: d(:), rbar(:)

!     Local variable

REAL (dp), PARAMETER :: zero = 0.0_dp

!     Some checks.

ier = 0
IF (np < 1) ier = 1
IF (nrbar < np*(np-1)/2) ier = ier + 2
IF (ier /= 0) RETURN

d(1:np) = zero
rbar(1:nrbar) = zero

RETURN
END SUBROUTINE clear



SUBROUTINE regcf(np, nrbar, d, rbar, thetab, tol, beta, nreq, ier)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Modified version of AS75.4 to calculate regression coefficients for the
!     first NREQ variables, given an orthogonal reduction from AS75.1.

INTEGER, INTENT(IN)       :: np, nrbar, nreq
INTEGER, INTENT(OUT)      :: ier
REAL (dp), INTENT(IN)     :: rbar(:), thetab(:), tol(:)
REAL (dp), INTENT(IN OUT) :: d(:)
REAL (dp), INTENT(OUT)    :: beta(:)

!     Local variables

INTEGER              :: i, j, nextr
REAL (dp), PARAMETER :: zero = 0.0_dp

!     Some checks.

ier = 0
IF (np < 1) ier = 1
IF (nrbar < np*(np-1)/2) ier = ier + 2
IF (nreq < 1 .OR. nreq > np) ier = ier + 4
IF (ier /= 0) RETURN

DO i = nreq, 1, -1
  IF (SQRT(d(i)) < tol(i)) THEN
    beta(i) = zero
    d(i) = zero
    CYCLE
  END IF
  beta(i) = thetab(i)
  nextr = (i-1) * (np+np-i)/2 + 1
  DO j = i+1, nreq
    beta(i) = beta(i) - rbar(nextr) * beta(j)
    nextr = nextr + 1
  END DO
END DO

RETURN
END SUBROUTINE regcf


SUBROUTINE bksub2(rbar, k, rhs, soln)

!     Solves  R'(sqrt(D).z) = x where (sqrt(D).z) = SOLN and x = RHS.
!     RBAR is an upper-triangular matrix with implicit 1's on it's diagonal,
!     stored by rows.

INTEGER, INTENT(IN)    :: k
REAL (dp), INTENT(IN)  :: rbar(:), rhs(:)
REAL (dp), INTENT(OUT) :: soln(:)

!     Local variables

INTEGER   :: row, col, pos
REAL (dp) :: temp

soln(1) = rhs(1)
DO row = 2, k
  temp = rhs(row)
  pos = row - 1
  DO col = 1, row-1
    temp = temp - rbar(pos) * soln(col)
    pos = pos + k - col - 1
  END DO
  soln(row) = temp
END DO

RETURN
END SUBROUTINE bksub2



SUBROUTINE delta(k, xj, xl, xbari, xbark, ni, nk, d, rbar, fn_val)

!     Calculate the delta function for the swap of case J in block I
!     with case L in block K.   Uses the method of Cook & Nachtsheim.

INTEGER, INTENT(IN)    :: k, ni, nk
REAL (dp), INTENT(IN ) :: xj(:), xl(:), xbari(:), xbark(:), d(:), rbar(:)
REAL (dp), INTENT(OUT) :: fn_val

!     Local variables

INTEGER              :: i
REAL (dp)            :: const, temp, e11, e12, e21, e22
REAL (dp)            :: z(k,3), a(k), b(k), diff(k)
REAL (dp), PARAMETER :: one = 1.0_dp, two = 2.0_dp

!     Calculate vectors DIFF, A & B.

const = two - one/ni - one/nk
DO i = 1, k
  temp = xj(i) - xl(i)
  diff(i) = -temp
  a(i) = temp - xbari(i) + xbark(i)
  b(i) = a(i) - const * temp
END DO

!     Calculate the z-vectors by back-substitution.
!     Z1 for A, Z2 for B and Z3 for DIFF.
!     N.B. The solutions returned from BKSUB2 have the I-th element
!     multiplied by sqrt(D(I)).

CALL bksub2(rbar, k, a, z(:,1))
CALL bksub2(rbar, k, b, z(:,2))
CALL bksub2(rbar, k, diff, z(:,3))

!     Calculate the elements E11, E12, E21 & E22 as dot-products of the
!     appropriate z-vectors.

e11 = DOT_PRODUCT( z(:,3), z(:,1)/d )
e12 = DOT_PRODUCT( z(:,3), z(:,3)/d )
e21 = DOT_PRODUCT( z(:,2), z(:,1)/d )
e22 = DOT_PRODUCT( z(:,2), z(:,3)/d )

!     Return the determinant of the matrix:   E11+1    E12
!     (minus 1)                                E21    E22+1

fn_val = (e11+one) * (e22+one) - e12 * e21 - one

RETURN
END SUBROUTINE delta


END MODULE d_optimal_design
