SUBROUTINE gllm(ni, nid, nj, nk, nkp, ji, y, c, conv, e, f, cspr, cslr, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-04-26  Time: 11:05:17

! N.B. Arguments W & V have been removed - they were workspaces.

!    Algorithm AS 207 Appl. Statist. (1984) vol.33, no.3, pp. 358-362.
!    by Michael Haber

!    Fitting a generalized log-linear model
!    to fully or partially classified frequencies.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)        :: ni, nid, nj, nk, nkp, ji(ni)
REAL (dp), INTENT(IN)      :: y(nj)
REAL (dp), INTENT(IN OUT)  :: c(nid,nkp)
REAL (dp), INTENT(IN)      :: conv
REAL (dp), INTENT(OUT)     :: e(ni), f(nj), cspr, cslr
INTEGER, INTENT(OUT)       :: ifault

REAL (dp)  :: w(4,ni), v(2,nkp)
REAL (dp)  :: cmin, csum, ctmax, ff, yy
INTEGER    :: i, j, k, nkk
LOGICAL    :: lconv
REAL (dp), PARAMETER  :: eps = 0.00001_dp, zero = 0.0_dp, one = 1.0_dp

ifault = 1

!       Check the ji array

DO  i = 1, ni
  IF (ji(i) < 1 .OR. ji(i) > nj) RETURN
END DO
ifault = 0

!       initialize

e(1:ni) = one
w(3,1:nj) = zero

!       standardize the C matrix

cmin = zero
DO  i = 1, ni
  DO  k = 1, nk
    cmin = MIN(cmin,c(i,k))
  END DO
END DO
IF (cmin /= zero) THEN
    c(1:ni,1:nk) = c(1:ni,1:nk) - cmin
END IF

80 ctmax = zero
DO  i = 1, ni
  csum = zero
  DO  k = 1, nk
    csum = csum + c(i,k)
  END DO
  ctmax = MAX(ctmax,csum)
  w(4,i) = csum
END DO

IF (ctmax <= eps) THEN
  ifault = 2
  RETURN
END IF
IF (ABS(ctmax-one) > eps) THEN
  DO  i = 1, ni
    DO  k = 1, nk
      c(i,k) = c(i,k) / ctmax
    END DO
  END DO
  GO TO 80
END IF
DO  i = 1, ni
  IF (ABS(w(4,i)-one) > eps) GO TO 140
END DO
nkk = nk
GO TO 160
140 nkk = nkp
DO  i = 1, ni
  c(i,nkk) = one - w(4,i)
END DO

!       Enter the EM algorithm

160 f(1:nj) = zero
DO  i = 1, ni
  j = ji(i)
  f(j) = f(j) + e(i)
END DO

!       Check for convergence

lconv = .true.
DO  j = 1, nj
  IF (ABS(f(j)-w(3,j)) > conv) lconv = .false.
  w(3,j) = f(j)
END DO
IF (.NOT.lconv) THEN
  
  DO  i = 1, ni
    j = ji(i)
    w(1,i) = y(j)
    IF (f(j) > eps) w(1,i) = e(i) * y(j) / f(j)
  END DO
  DO  k = 1, nkk
    v(1,k) = SUM( c(1:ni,k) * w(1,1:ni) )
  END DO
  
!       Enter the IPF algorithm
  
  230 w(2,1:ni) = e(1:ni)
  DO  k = 1, nkk
    v(2,k) = SUM( c(1:ni,k) * e(1:ni) )
    DO  i = 1, ni
      IF (c(i,k) > eps .AND. v(2,k) > eps) e(i) = e(i) * (v(1,k)/  &
          v(2,k)) ** c(i,k)
    END DO
  END DO
  DO  i = 1, ni
    IF (ABS(e(i)-w(2,i)) > conv) GO TO 230
  END DO
  GO TO 160
END IF

!       calculate the goodness-of-fit statistics

cspr = zero
cslr = zero
DO  j = 1, nj
  yy = y(j)
  ff = f(j)
  IF (ff > eps) THEN
    cspr = cspr + (yy-ff) ** 2 / ff
    IF (yy > eps) THEN
      cslr = cslr + yy * LOG(yy/ff)
    END IF
  END IF
END DO
cslr = 2.0D0 * cslr
RETURN
END SUBROUTINE gllm
