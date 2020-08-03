PROGRAM freq2d
! Generate pairs of random numbers and put into 256 x 256 bins.

IMPLICIT NONE
INTEGER              :: i, j, k, l, lower, ndf, nrows, ncols, upper
REAL                 :: x(2), aver, stdev, chisq
INTEGER, ALLOCATABLE :: freq(:,:), seed(:)
CHARACTER (LEN=200)  :: text

WRITE(*, *) 'Enter number of rows & columns (e.g. 24  64): '
READ(*, *) nrows, ncols
ALLOCATE( freq(0:nrows-1,0:ncols-1) )

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
CALL RANDOM_SEED(get=seed)
WRITE(*, *)'Old random number seeds: ', seed

WRITE(*, '(a, i4, a)') ' Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

! Generate 4096 x 4096 pairs of random numbers

freq = 0
DO i = 1, 4096
  DO j = 1, 4096
    CALL RANDOM_NUMBER(x)
    k = nrows * x(1)
    l = ncols * x(2)
    freq(k,l) = freq(k,l) + 1
  END DO
END DO

! Average & variance of number in each cell = 4096 x 4096 / (nrows x ncols)
aver = 4096. * 4096. / (nrows * ncols)
stdev = SQRT(aver)
upper = aver + 3.*stdev + 0.5
lower = aver - 3.*stdev + 0.5

chisq = 0.0
DO k = 0, nrows-1
  text = ' '
  DO l = 0, ncols-1
    IF (freq(k,l) >= lower .AND. freq(k,l) <= upper) THEN
      text(l+3:l+3) = '+'
    ELSE IF (freq(k,l) > upper) THEN
      text(l+3:l+3) = CHAR(219)
    END IF
    chisq = chisq + (freq(k,l) - aver)**2
  END DO
  WRITE(*, '(a)') TRIM(text)
END DO

WRITE(*, *) '  If most of the cells are NOT pluses, its a bad generator'

chisq = chisq / aver
ndf = nrows * ncols - 1
WRITE(*, '(a, f10.1, a, i6, a)')  &
         ' Chi-squared =', chisq, ' with', ndf, ' deg. of freedom'
chisq = ndf + 2.0*SQRT(2.0*ndf)
WRITE(*, '(a, f10.1)') ' Chi-squared should be less than', chisq

STOP
END PROGRAM freq2d
