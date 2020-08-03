PROGRAM driver

!     Use DOPT to find good incomplete-block designs for fitting
!      p   q
!     2 x 3  fractional factorials in blocks which which may
!     be of unequal sizes.
!
!     Either p or q must be > 0, and there may be only one block.
!
!     The design is for fitting all main effects, all 2-factor
!     interactions amongst the factors at 2 levels, a quadratic surface
!     in the factors at 3 levels, and all interactions between the
!     linear terms of the quadratic surface and the factors at 2 levels.
!
!     For large designs it may be necessary to increase one or more of
!     the values in the parameter statement below.   Notice that there
!     are some very large 2-D arrays.
!     maxf = maximum number of factors (p + q)
!     maxb = maximum number of blocks
!     maxcol = maximum number of columns in the design matrix
!     mxcand = maximum number of candidate points

USE d_optimal_design

IMPLICIT NONE

INTEGER, ALLOCATABLE   :: picked(:), blksiz(:), in(:), design(:), ind(:), &
                          seed(:)
INTEGER                :: p, q, nfact, nblock, n, kin, kfull, nrbar,  &
                          ncand, ifault, point, i, j, pos, nrep, i1, i2
REAL, ALLOCATABLE      :: x(:,:)
REAL (dp), ALLOCATABLE :: d(:), rbar(:), xx(:), tol(:), wk(:), zpz(:,:)
REAL (dp)              :: dbest, zero = 0.D0, stdet, lndet
LOGICAL                :: rstart = .true.
CHARACTER (LEN = 1)    :: bel

WRITE(*, *)'   +++  Driver program for DOPT program  +++'
WRITE(*, *)
bel = CHAR(7)

10 WRITE(*, *)'Enter no. of factors with 2 levels: '
READ(*, *) p
IF (p < 0) THEN
  WRITE(*, 900) bel
  900 FORMAT(' ', a, '  **  Number must be positive or zero **')
  GO TO 10
END IF
WRITE(*, *)'Enter no. of factors with 3 levels: '
READ(*, *) q
IF (q < 0) THEN
  WRITE(*, 900) bel
  GO TO 10
END IF
nfact = p + q
IF (nfact <= 0) THEN
  WRITE(*, *) bel, '  **  No. of factors must be > 0 **'
  GO TO 10
END IF

!     Ask for the block specification

WRITE(*, *) 'Enter no. of blocks: '
READ(*, *) nblock
IF (nblock < 0) THEN
  WRITE(*, *) bel, '  ** Negative no. of blocks **'
  GO TO 10
END IF
nblock = MAX(1, nblock)
ALLOCATE( blksiz(nblock), in(nblock) )
in = 0

IF (nblock == 1) THEN
  WRITE(*, *)'Enter size of the experiment: '
ELSE
  WRITE(*, *)'Enter size of each block: '
END IF
READ(*, *) blksiz(1:nblock)
n = SUM( blksiz(1:nblock) )

kin = nfact + q + nfact*(nfact-1)/2
kfull = kin + nblock
WRITE(*, 980) kfull
980 FORMAT(' No. of parameters in model = ', i5)

!     Check that the number of design points is at least equal to the
!     number of parameters

IF (n < kfull) THEN
  WRITE(*, *) bel, '** Design too small to fit model **'
  DEALLOCATE( blksiz, in )
  GO TO 10
END IF
nrbar = kfull * (kfull - 1) / 2
ncand = (2**p) * (3**q)
WRITE(*, *) 'No. of candidate points = ', ncand

WRITE(*, *) 'How many tries?: '
READ(*, *) nrep

CALL RANDOM_SEED(size=i)
ALLOCATE( seed(i) )
WRITE(*, *)'Enter', i, ' integers for random no. seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

!     Open a file for the output
!     Default file name:  DESIGN.OPT

OPEN(10, FILE='DESIGN.OPT')
WRITE(10, 1000) p, q
1000 FORMAT('Output from DOPT for fractional factorial'/   &
            i2, ' factors at 2 levels, and', i3, ' factors at 3 levels')
IF (nblock .GT. 1) THEN
  WRITE(10, '(a, 10i4)') 'Block sizes: ', blksiz
ELSE
  WRITE(10, '(a, i5)') 'No. of experimental runs = ', n
END IF
WRITE(10, 1020) kfull, ncand, nrep, seed
1020 FORMAT('No. of parameters in full model incl. blocks = ', i4/ &
            'No. of candidate points = ', i6/                      &
            'No. of tries to find optimum design = ', i5/          &
            ('Starting random number seeds = ', 7I10) )
DEALLOCATE( seed )

!     End of user input.

!---------------------------------------------------------------------
!     Generate the candidate design points.
!     ind(i) stores the level of the i-th factor.
!     For the first p factors, ind(i) = -1 or +1
!     For the next  q factors, ind(i) = -1, 0 or +1

point = 1
ALLOCATE( ind(nfact), x(ncand,kfull) )
ind = -1

!     Calculate row of X for the current candidate point
!     1. Main effects

60 x(point, 1:nfact) = ind(1:nfact)
pos = nfact

!     2. Quadratic terms

DO i = 1, q
  x(point, pos+i) = ind(i+p)**2
END DO

!     3. Interaction terms

pos = pos + q
DO i = 1, nfact-1
  DO j = i+1, nfact
    pos = pos + 1
    x(point, pos) = x(point, i) * x(point, j)
  END DO
END DO

!     Find next candidate point by incrementing indices

i = 1
120 IF (i <= p) THEN
  ind(i) = - ind(i)                        ! 2-level factor
  IF (ind(i) == -1) THEN
    i = i + 1
    GO TO 120
  END IF
ELSE
  ind(i) = ind(i) + 1                      ! 3-level factor
  IF (ind(i) > 1) THEN
    ind(i) = -1
    i = i + 1
    IF (i > nfact) GO TO 150
    GO TO 120
  END IF
END IF

point = point + 1
GO TO 60
!---------------------------------------------------------------------
!
!     Call DOPT to try to find suitable designs

150 DEALLOCATE( ind )
ALLOCATE( d(kfull), rbar(nrbar), picked(n), xx(kfull), tol(kfull),         &
          zpz(ncand,nblock), wk(kfull), design(n) )
dbest = zero
DO i = 1, nrep
  CALL dopt(x, ncand, ncand, kin, n, nblock, in, blksiz, kfull, rstart,    &
            nrbar, d, rbar, picked, lndet, xx, tol, zpz, wk, ifault)
  IF (ifault /= 0) THEN
    WRITE(*, *) bel, 'IFAULT = ', ifault
  ELSE
    WRITE(*, 950) i, lndet
    WRITE(10, 950) i, lndet
    950 FORMAT(' Try no. ', i4, '     Log of determinant = ', f9.3)
  END IF
  IF (lndet > dbest) THEN
    design = picked
    dbest = lndet
  END IF
END DO
!---------------------------------------------------------------------

!     Output the design with the highest determinant

WRITE(*, *)
WRITE(10, *)
stdet = EXP(dbest - kfull*LOG(REAL(n)))
WRITE(*, 990) dbest, stdet
WRITE(10, 990) dbest, stdet
990 FORMAT(' Max. log det. = ', g13.5, '     Std. det = ', g13.5/    &
           ' Design:')
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

STOP
END PROGRAM driver
