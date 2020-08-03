! Code converted using TO_F90 by Alan Miller
! Date: 2002-03-06  Time: 08:36:31

! Converted with permission, from the F77 code in the CERN MATHLIB library.
! $Id: assndx.F,v 1.1.1.1 1996/04/01 15:02:49 mclareni Exp $

! $Log: assndx.F,v $
! Revision 1.1.1.1  1996/04/01 15:02:49  mclareni
! Mathlib gen/H (H301)
! Author: F. Bourgeois, 15 February 1994

SUBROUTINE assndx(mode, a, n, m, k, sum)
! N.B. Arguments IDA, IW & IDW have been removed.

! If MODE = 1, then it finds k(1), k(2), ..., k(n) to minimize
!        S = Sum(i=1, .., n) a(i, k(i))
! If MODE = 2,  then it finds k(1), k(2), ..., k(m) to minimize
!        S = Sum(j=1, .., m) a(k(j), j)
! given the array a(n,m).

! References:
! Munkres, J. (1957) `Algorithms for the assignment and transportation problems',
!                    J. SIAM, vol.5, 32-38.
! Silver, R. (1960) `An algorithm for the assignment problem', Comm. ACM, vol.3,
!                   605-606.   The algorithm (CACM 27) is in Algol.

IMPLICIT NONE

INTEGER, INTENT(IN)   :: mode
REAL, INTENT(IN OUT)  :: a(:,:)
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: m
INTEGER, INTENT(OUT)  :: k(:)
REAL, INTENT(OUT)     :: sum

LOGICAL  :: lsw
INTEGER  :: i, icbl, icl, icl0, iflag, imax, imin, ipp, irl, irs, &
            j, j1, jsv, new
REAL     :: rmin
INTEGER, ALLOCATABLE  :: iw(:,:)

IF (n < 1 .OR. m < 1) THEN
  WRITE(*, '(a, 2i8)') ' ** Error in call to ASSNDX; m, n = ', m, n
  RETURN
END IF

imax = MAX(n,m)
imin = MIN(n,m)
ALLOCATE( iw(imax,6) )
sum = 0.0
IF (n <= m) THEN
  DO  i = 1, n
    rmin = a(i,1)
    DO  j = 1, m
      rmin = MIN(rmin, a(i,j))
    END DO
    sum = sum + rmin
    a(i,1:m) = a(i,1:m) - rmin
  END DO
END IF
IF (n >= m) THEN
  DO  j = 1, m
    rmin = a(1,j)
    DO  i = 1, n
      rmin = MIN(rmin,a(i,j))
    END DO
    sum = sum + rmin
    a(1:n,j) = a(1:n,j) - rmin
  END DO
END IF

DO  i = 1, imax
  k(i) = 0
  iw(i,1) = 0
END DO

loop90:  DO  i = 1, n
  DO  j = 1, m
    IF (a(i,j)+iw(j,1) == 0) THEN
      k(i) = j
      iw(j,1) = i
      CYCLE loop90
    END IF
  END DO
END DO loop90

100 iflag = n
irl = 0
icl = 0
irs = 1

DO  i = 1, n
  iw(i,5) = 0
  IF (k(i) == 0) THEN
    irl = irl + 1
    iw(irl,6) = i
    iw(i,5) = -1
    iflag = iflag - 1
  END IF
END DO
IF (iflag == imin) THEN
  IF (mode == 2) k(1:imax) = iw(1:imax,1)
  RETURN
END IF

iw(1:m,4) = 0

140 i = iw(irs,6)
irs = irs + 1
DO  j = 1, m
  IF (a(i,j)+iw(j,4) == 0) THEN
    iw(j,4) = i
    icl = icl + 1
    iw(icl,2) = j
    NEW = iw(j,1)
    IF (NEW == 0) THEN
      j1 = j
      DO
        iw(j1,1) = iw(j1,4)
        i = iw(j1,4)
        IF (k(i) == 0) THEN
          k(i) = j1
          GO TO 100
        END IF
        jsv = j1
        j1 = k(i)
        k(i) = jsv
      END DO
    END IF
    irl = irl + 1
    iw(irl,6) = NEW
    iw(NEW,5) = i
  END IF
END DO
IF (irs <= irl) GO TO 140

lsw = .true.
icl0 = icl
icbl = 0
DO  j = 1, m
  IF (iw(j,4) == 0) THEN
    icbl = icbl + 1
    iw(icbl,3) = j
  END IF
END DO
rmin = a(iw(1,6),iw(1,3))
DO  i = 1, irl
  DO  j = 1, icbl
    rmin = MIN(rmin, a(iw(i,6), iw(j,3)))
  END DO
END DO
sum = sum + rmin * (irl+icbl-imax)

DO  i = 1, n
  IF (iw(i,5) == 0) THEN
    DO  ipp = 1, icl0
      a(i,iw(ipp,2)) = a(i,iw(ipp,2)) + rmin
    END DO
    CYCLE
  END IF
  DO  ipp = 1, icbl
    NEW = iw(ipp,3)
    a(i,NEW) = a(i,NEW) - rmin
    IF (lsw.AND.a(i,NEW)+iw(NEW,4) == 0) THEN
      iw(NEW,4) = i
      IF (iw(NEW,1) == 0) THEN
        j1 = NEW
        lsw = .false.
      ELSE
        icl = icl + 1
        iw(icl,2) = NEW
        irl = irl + 1
        iw(irl,6) = iw(NEW,1)
      END IF
    END IF
  END DO
END DO

IF (lsw) THEN
  DO  i = icl0 + 1, icl
    iw(iw(iw(i,2),1),5) = iw(i,2)
  END DO
  GO TO 140
ELSE
  DO
    iw(j1,1) = iw(j1,4)
    i = iw(j1,4)
    IF (k(i) == 0) THEN
      k(i) = j1
      GO TO 100
    END IF
    jsv = j1
    j1 = k(i)
    k(i) = jsv
  END DO
END IF

RETURN
END SUBROUTINE assndx
