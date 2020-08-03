PROGRAM ibquad
!
!     Use DOPT to find good incomplete-block designs for fitting
!     quadratic surfaces.
!
USE d_optimal_design

IMPLICIT NONE

INTEGER, PARAMETER :: maxf=6, maxb=10, maxcol=maxf*(maxf+3)/2, &
                      kmax=maxcol+maxb, nrmax=kmax*(kmax-1)/2, &
                      mxcand=3**maxf, nmax=100
INTEGER            :: nfact, nblock, blksiz(maxb), n, kin, kfull, nrbar, &
                      ncand, picked(nmax), ifault, point, i, j, ind(maxcol), &
                      pos, nrep, design(nmax), i1, i2
REAL (dp)          :: d(kmax), rbar(nrmax), lndet, xx(kmax), &
                      tol(kmax), zpz(mxcand,maxb), wk(kmax), dbest, stdet
REAL               :: x(mxcand,maxcol)
LOGICAL, PARAMETER :: rstart = .TRUE.
INTEGER, PARAMETER :: in(maxb) = (/ (0,i=1,maxb) /)
INTEGER, ALLOCATABLE :: seed(:)
CHARACTER (LEN=1)    :: bel
REAL (dp), PARAMETER :: zero = 0.0_dp


WRITE(*, *)'   +++   Incomplete blocks for quad. surfaces   +++'
bel = CHAR(7)
10 WRITE(*, *) 'Enter no. of factors: '
READ(*, *) nfact
IF (nfact > maxf .OR. nfact < 2) THEN
  WRITE(*, 900) bel
  900   FORMAT(' ', a, '** Illegal value entered **')
  GO TO 10
END IF

WRITE(*, *) 'Enter no. of blocks: '
READ(*, *) nblock
IF (nblock > maxb .OR. nblock < 0) THEN
  WRITE(*, 900) bel
  GO TO 10
END IF

WRITE(*, *)'Enter size of each block: '
READ(*, *) blksiz(1:nblock)
n = SUM( blksiz(1:nblock) )
IF (n > nmax) THEN
  WRITE(*, *) bel, '** Sorry, design is too large **'
  GO TO 10
END IF

WRITE(*, *)'How many tries?: '
READ(*, *) nrep

CALL RANDOM_SEED(size=i)
ALLOCATE( seed(i) )
WRITE(*, *)'Enter', i, ' integers for random no. seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

kin = nfact * (nfact + 3) / 2
kfull = kin + nblock
IF (n < kfull) THEN
  WRITE(*, *) bel, '** Design too small to fit model **'
  GO TO 10
END IF
nrbar = kfull * (kfull - 1) / 2
ncand = 3**nfact


!     Open a file for the output
!     Default file name:  DESIGN.OPT

OPEN(10, FILE='DESIGN.OPT')
WRITE(10, 1000) nfact
1000 FORMAT('Output from DOPT for fitting quadratic surfaces'/   &
            i2, ' factors at 2 levels')
IF (nblock .GT. 1) THEN
  WRITE(10, '(a, 10i4)') 'Block sizes: ', blksiz
ELSE
  WRITE(10, '(a, i5)') 'No. of experimental runs = ', n
END IF

!     Generate the candidate design points.
!     ind(i) stores the level (-1,0,1) of the i-th factor.

point = 1
ind(1:nfact) = -1
GO TO 55
30 i = 1
40 ind(i) = ind(i) + 1
IF (ind(i) > 1) THEN
  IF (i == nfact) GO TO 90
  ind(1:i) = -1
  i = i + 1
  GO TO 40
END IF
!
!     Calculate row of X.
!
55 x(point,1:nfact) = ind(1:nfact)
pos = nfact + 1
DO i = 1, nfact
  DO j = 1, i
    x(point,pos) = ind(i) * ind(j)
    pos = pos + 1
  END DO
END DO
point = point + 1
GO TO 30

90 dbest = zero
DO i = 1, nrep
  CALL dopt(x, mxcand, ncand, kin, n, nblock, in, blksiz, kfull, &
            rstart, nrbar, d, rbar, picked, lndet, xx, tol, zpz, wk, &
  ifault)
  IF (ifault /= 0) THEN
    WRITE(*, *) bel, 'IFAULT = ', ifault
    ELSE
    WRITE(*, *) i, ' Log Det. = ', lndet
  END IF
  IF (lndet > dbest) THEN
    design(1:n) = picked(1:n)
    dbest = lndet
  END IF
END DO

WRITE(*, *)
stdet = EXP(dbest - kfull*LOG(REAL(n, KIND=dp)))
WRITE(*, *) 'Max. log det. = ', dbest, '   Std. det = ', stdet
WRITE(*, *) 'Design:'
WRITE(10, *) 'Max. log det. = ', dbest, '   Std. det = ', stdet
WRITE(10, *) 'Design:'
i2 = 0
DO i = 1, nblock
  i1 = i2 + 1
  i2 = i2 + blksiz(i)
  WRITE(*, 960) i, (j, j=1,nfact)
  WRITE(10, 960) i, (j, j=1,nfact)
  960 FORMAT(' BLOCK', i5 / t15, 'F A C T O R' / '   ', 15I5)
  DO j = i1, i2
    WRITE(*, 970) x(design(j), 1:nfact)
    WRITE(10, 970) x(design(j), 1:nfact)
    970 FORMAT('    ', 15F5.0)
  END DO
END DO
WRITE(*, *)
GO TO 10

STOP
END PROGRAM ibquad
