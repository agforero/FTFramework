MODULE kmeans
IMPLICIT NONE

!     Define BIG to be a very large positive number

REAL, PARAMETER, PRIVATE  :: zero = 0.0, one = 1.0, big = HUGE(zero)

PUBLIC   :: kmns
PRIVATE  :: optra, qtran


CONTAINS


SUBROUTINE kmns(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d,  &
                itran, live, iter, wss, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-11-01  Time: 10:00:51

!     ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1

!     Divide M points in N-dimensional space into K clusters so that
!     the within cluster sum of squares is minimized.

REAL, INTENT(IN)         :: a(:,:)   ! a(m,n)
INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN OUT)     :: c(:,:)   ! c(k,n)
INTEGER, INTENT(IN)      :: k
INTEGER, INTENT(OUT)     :: ic1(:)
INTEGER, INTENT(OUT)     :: ic2(:)
INTEGER, INTENT(OUT)     :: nc(:)
REAL, INTENT(OUT)        :: an1(:)
REAL, INTENT(OUT)        :: an2(:)
INTEGER, INTENT(OUT)     :: ncp(:)
REAL, INTENT(IN OUT)     :: d(:)
INTEGER, INTENT(OUT)     :: itran(:)
INTEGER, INTENT(IN OUT)  :: live(:)
INTEGER, INTENT(IN)      :: iter
REAL, INTENT(OUT)        :: wss(:)
INTEGER, INTENT(OUT)     :: ifault

! Local variables

INTEGER  :: i, ii, ij, il, indx, j, l
REAL     :: aa, da, db, dc, dt(2), temp


ifault = 3
IF (k <= 1 .OR. k >= m) RETURN

!     For each point I, find its two closest centres, IC1(I) and
!     IC2(I).     Assign it to IC1(I).

DO  i = 1, m
  ic1(i) = 1
  ic2(i) = 2
  DO  il = 1, 2
    dt(il) = zero
    DO  j = 1, n
      da = a(i,j) - c(il,j)
      dt(il) = dt(il) + da*da
    END DO
  END DO
  IF (dt(1) > dt(2)) THEN
    ic1(i) = 2
    ic2(i) = 1
    temp = dt(1)
    dt(1) = dt(2)
    dt(2) = temp
  END IF

  loop50:  DO  l = 3, k
    db = zero
    DO  j = 1, n
      dc = a(i,j) - c(l,j)
      db = db + dc*dc
      IF (db >= dt(2)) CYCLE loop50
    END DO
    IF (db >= dt(1)) THEN
      dt(2) = db
      ic2(i) = l
    ELSE
      dt(2) = dt(1)
      ic2(i) = ic1(i)
      dt(1) = db
      ic1(i) = l
    END IF
  END DO loop50
END DO

!     Update cluster centres to be the average of points contained
!     within them.

DO  l = 1, k
  nc(l) = 0
  c(l,1:n) = zero
END DO
DO  i = 1, m
  l = ic1(i)
  nc(l) = nc(l) + 1
  c(l,1:n) = c(l,1:n) + a(i,1:n)
END DO

!     Check to see if there is any empty cluster at this stage

DO  l = 1, k
  IF (nc(l) == 0) THEN
    ifault = 1
    RETURN
  END IF
  aa = nc(l)
  c(l,1:n) = c(l,1:n) / aa
  
!     Initialize AN1, AN2, ITRAN & NCP
!     AN1(L) = NC(L) / (NC(L) - 1)
!     AN2(L) = NC(L) / (NC(L) + 1)
!     ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
!              = 0 otherwise
!     In the optimal-transfer stage, NCP(L) stores the step at which
!     cluster L is last updated.
!     In the quick-transfer stage, NCP(L) stores the step at which
!     cluster L is last updated plus M.
  
  an2(l) = aa / (aa + one)
  an1(l) = big
  IF (aa > one) an1(l) = aa / (aa - one)
  itran(l) = 1
  ncp(l) = -1
END DO
indx = 0
DO  ij = 1, iter
  
!     In this stage, there is only one pass through the data.   Each
!     point is re-allocated, if necessary, to the cluster that will
!     induce the maximum reduction in within-cluster sum of squares.
  
  CALL optra(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, indx)
  
!     Stop if no transfer took place in the last M optimal transfer steps.
  
  IF (indx == m) GO TO 150
  
!     Each point is tested in turn to see if it should be re-allocated
!     to the cluster to which it is most likely to be transferred,
!     IC2(I), from its present cluster, IC1(I).   Loop through the
!     data until no further change is to take place.
  
  CALL qtran(a, m, n, c, ic1, ic2, nc, an1, an2, ncp, d, itran, indx)
  
!     If there are only two clusters, there is no need to re-enter the
!     optimal transfer stage.
  
  IF (k == 2) GO TO 150
  
!     NCP has to be set to 0 before entering OPTRA.
  
  ncp(1:k) = 0
END DO

!     Since the specified number of iterations has been exceeded, set
!     IFAULT = 2.   This may indicate unforeseen looping.

ifault = 2

!     Compute within-cluster sum of squares for each cluster.

150 DO  l = 1, k
  wss(l) = zero
  c(l,1:n) = zero
END DO
DO  i = 1, m
  ii = ic1(i)
  c(ii,1:n) = c(ii,1:n) + a(i,1:n)
END DO
DO  j = 1, n
  c(1:k,j) = c(1:k,j) / REAL(nc(l))
  DO  i = 1, m
    ii = ic1(i)
    da = a(i,j) - c(ii,j)
    wss(ii) = wss(ii) + da*da
  END DO
END DO

RETURN
END SUBROUTINE kmns



SUBROUTINE optra(a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d,  &
                 itran, live, indx)

!     ALGORITHM AS 136.1  APPL. STATIST. (1979) VOL.28, NO.1

!     This is the optimal transfer stage.

!     Each point is re-allocated, if necessary, to the cluster that will
!     induce a maximum reduction in the within-cluster sum of squares.

REAL, INTENT(IN)         :: a(:,:)   ! a(m,n)
INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN OUT)     :: c(:,:)   ! c(k,n)
INTEGER, INTENT(IN)      :: k
INTEGER, INTENT(IN OUT)  :: ic1(:)
INTEGER, INTENT(IN OUT)  :: ic2(:)
INTEGER, INTENT(IN OUT)  :: nc(:)
REAL, INTENT(IN OUT)     :: an1(:)
REAL, INTENT(IN OUT)     :: an2(:)
INTEGER, INTENT(IN OUT)  :: ncp(:)
REAL, INTENT(OUT)        :: d(:)
INTEGER, INTENT(IN OUT)  :: itran(:)
INTEGER, INTENT(OUT)     :: live(:)
INTEGER, INTENT(OUT)     :: indx

! Local variables

INTEGER  :: i, j, l, l1, l2, ll
REAL     :: al1, al2, alt, alw, da, db, dc, dd, de, df, rr, r2

!     If cluster L is updated in the last quick-transfer stage, it
!     belongs to the live set throughout this stage.   Otherwise, at
!     each step, it is not in the live set if it has not been updated
!     in the last M optimal transfer steps.

DO  l = 1, k
  IF (itran(l) == 1) live(l) = m + 1
END DO
DO  i = 1, m
  indx = indx + 1
  l1 = ic1(i)
  l2 = ic2(i)
  ll = l2
  
!     If point I is the only member of cluster L1, no transfer.
  
  IF (nc(l1) == 1) THEN
    IF (indx == m) RETURN
    CYCLE
  END IF
  
!     If L1 has not yet been updated in this stage, no need to re-compute D(I).
  
  IF (ncp(l1) == 0) GO TO 30
  de = zero
  DO  j = 1, n
    df = a(i,j) - c(l1,j)
    de = de + df*df
  END DO
  d(i) = de * an1(l1)
  
!     Find the cluster with minimum R2.
  
  30 da = zero
  DO  j = 1, n
    db = a(i,j) - c(l2,j)
    da = da + db*db
  END DO
  r2 = da * an2(l2)
  loop60:  DO  l = 1, k
    
!    If I >= LIVE(L1), then L1 is not in the live set.   If this is true, we
!    only need to consider clusters that are in the live set for possible
!    transfer of point I.  Otherwise, we need to consider all possible clusters.
    
    IF (i >= live(l1) .AND. i >= live(l) .OR. l == l1 .OR.  &
        l == ll) CYCLE loop60
    rr = r2 / an2(l)
    dc = zero
    DO  j = 1, n
      dd = a(i,j) - c(l,j)
      dc = dc + dd*dd
      IF (dc >= rr) CYCLE loop60
    END DO
    r2 = dc * an2(l)
    l2 = l
  END DO loop60
  IF (r2 >= d(i)) THEN
  
!     If no transfer is necessary, L2 is the new IC2(I).
  
    ic2(i) = l2
    IF (indx == m) RETURN
  
!     Update cluster centres, LIVE, NCP, AN1 & AN2 for clusters L1 and
!     L2, and update IC1(I) & IC2(I).

  ELSE
    indx = 0
    live(l1) = m + i
    live(l2) = m + i
    ncp(l1) = i
    ncp(l2) = i
    al1 = nc(l1)
    alw = al1 - one
    al2 = nc(l2)
    alt = al2 + one
    DO  j = 1, n
      c(l1,j) = (c(l1,j) * al1 - a(i,j)) / alw
      c(l2,j) = (c(l2,j) * al2 + a(i,j)) / alt
    END DO
    nc(l1) = nc(l1) - 1
    nc(l2) = nc(l2) + 1
    an2(l1) = alw / al1
    an1(l1) = big
    IF (alw > one) an1(l1) = alw / (alw - one)
    an1(l2) = alt / al2
    an2(l2) = alt / (alt + one)
    ic1(i) = l2
    ic2(i) = l1
  END IF
END DO
DO  l = 1, k
  
!     ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
!     decreased by M before re-entering OPTRA.
  
  itran(l) = 0
  live(l) = live(l) - m
END DO

RETURN
END SUBROUTINE optra



SUBROUTINE qtran(a, m, n, c, ic1, ic2, nc, an1, an2, ncp, d, itran, indx)

! N.B. Argument K has been removed.

!     ALGORITHM AS 136.2  APPL. STATIST. (1979) VOL.28, NO.1

!     This is the quick transfer stage.
!     IC1(I) is the cluster which point I belongs to.
!     IC2(I) is the cluster which point I is most likely to be
!         transferred to.
!     For each point I, IC1(I) & IC2(I) are switched, if necessary, to
!     reduce within-cluster sum of squares.  The cluster centres are
!     updated after each step.


REAL, INTENT(IN)         :: a(:,:)   ! a(m,n)
INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN OUT)     :: c(:,:)   ! c(k,n)
INTEGER, INTENT(IN OUT)  :: ic1(:)
INTEGER, INTENT(IN OUT)  :: ic2(:)
INTEGER, INTENT(IN OUT)  :: nc(:)
REAL, INTENT(IN OUT)     :: an1(:)
REAL, INTENT(IN OUT)     :: an2(:)
INTEGER, INTENT(IN OUT)  :: ncp(:)
REAL, INTENT(OUT)        :: d(:)
INTEGER, INTENT(OUT)     :: itran(:)
INTEGER, INTENT(OUT)     :: indx

! Local variables

INTEGER  :: i, icoun, istep, j, l1, l2
REAL     :: al1, al2, alt, alw, da, db, dd, de, r2

!     In the optimal transfer stage, NCP(L) indicates the step at which
!     cluster L is last updated.   In the quick transfer stage, NCP(L)
!     is equal to the step at which cluster L is last updated plus M.

icoun = 0
istep = 0
10 DO  i = 1, m
  icoun = icoun + 1
  istep = istep + 1
  l1 = ic1(i)
  l2 = ic2(i)
  
!     If point I is the only member of cluster L1, no transfer.
  
  IF (nc(l1) == 1) GO TO 60
  
!     If ISTEP > NCP(L1), no need to re-compute distance from point I to
!     cluster L1.   Note that if cluster L1 is last updated exactly M steps
!     ago, we still need to compute the distance from point I to cluster L1.
  
  IF (istep > ncp(l1)) GO TO 30
  da = zero
  DO  j = 1, n
    db = a(i,j) - c(l1,j)
    da = da + db*db
  END DO
  d(i) = da * an1(l1)
  
!     If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer of
!     point I at this step.
  
  30 IF (istep >= ncp(l1) .AND. istep >= ncp(l2)) GO TO 60
  r2 = d(i) / an2(l2)
  dd = zero
  DO  j = 1, n
    de = a(i,j) - c(l2,j)
    dd = dd + de*de
    IF (dd >= r2) GO TO 60
  END DO
  
!     Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters
!     L1 & L2.   Also update IC1(I) & IC2(I).   Note that if any
!     updating occurs in this stage, INDX is set back to 0.
  
  icoun = 0
  indx = 0
  itran(l1) = 1
  itran(l2) = 1
  ncp(l1) = istep + m
  ncp(l2) = istep + m
  al1 = nc(l1)
  alw = al1 - one
  al2 = nc(l2)
  alt = al2 + one
  DO  j = 1, n
    c(l1,j) = (c(l1,j) * al1 - a(i,j)) / alw
    c(l2,j) = (c(l2,j) * al2 + a(i,j)) / alt
  END DO
  nc(l1) = nc(l1) - 1
  nc(l2) = nc(l2) + 1
  an2(l1) = alw / al1
  an1(l1) = big
  IF (alw > one) an1(l1) = alw / (alw - one)
  an1(l2) = alt / al2
  an2(l2) = alt / (alt + one)
  ic1(i) = l2
  ic2(i) = l1
  
!     If no re-allocation took place in the last M steps, return.
  
  60 IF (icoun == m) RETURN
END DO
GO TO 10
END SUBROUTINE qtran

END MODULE kmeans
