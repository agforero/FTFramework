PROGRAM t_random

!     Test of module RANDOM
!     Latest revision - 4 January 1999
!     Programmer - Alan.Miller @ vic.cmis.csiro.au

USE random

IMPLICIT NONE

CHARACTER (LEN = 2)   :: option
INTEGER, ALLOCATABLE  :: ix(:), seed(:)
REAL, ALLOCATABLE     :: x(:), xmean(:), cov(:), chol_f(:)
INTEGER               :: i, n, ndf, n_binom, k, ier, pos1, pos2, which
REAL                  :: average, stdvn, shape, a, b, mu, p, pop_mean, pop_sd,   &
                         one = 1.0, zero = 0.0, sk, two = 2.0, pcntile(9), &
                         middle, xmax, xmin
LOGICAL               :: first

INTERFACE
  SUBROUTINE statistics(x, n, mean, median, xmax, xmin, stdev, percentile)
    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: n
    REAL, INTENT(IN)            :: x(n)
    REAL, OPTIONAL, INTENT(OUT) :: mean, median, xmax, xmin, stdev, percentile(9)
  END SUBROUTINE statistics
END INTERFACE

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)

DO
  WRITE(*, *)' 1.  Normal (Gaussian)                  2.  Gamma'
  WRITE(*, *)' 3.  Exponential                        4.  Beta'
  WRITE(*, *)' 5.  t                                  6.  Multivariate normal'
  WRITE(*, *)' 7.  Generalized inverse Gaussian       8.  Poisson'
  WRITE(*, *)' 9.  Binomial                          10.  Negative binomial'
  WRITE(*, *)'11.  von Mises                         12.  Cauchy'
  WRITE(*, *)'13.  Weibull                           14.  Chi-squared'
  WRITE(*, *)'15.  Random order of integers'
  WRITE(*, *)
  WRITE(*, *)' Q   Quit'
  WRITE(*, *)
  WRITE(*, *)'            Enter your choice: '
  READ(*, '(a)') option

  option = ADJUSTL(option)
  IF (option .EQ. 'q ' .OR. option .EQ. 'Q ') STOP

  SELECT CASE (option)
    CASE ('1 ')
      WRITE(*, *) 'Distribution:  Normal (Gaussian)'
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      DO i = 1, n
        x(i) = random_normal()
      END DO
      CALL statistics(x, n, average, middle, xmax, xmin, stdvn, pcntile)
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4, 2x, a, f9.4)')                    &
               'Mean = ', average, 'Median = ', middle, 'Sample st.dev. = ', stdvn
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4)')                                 &
               'Smallest x = ', xmin, 'Largest x = ', xmax
      WRITE(*, 900) pcntile
      900 FORMAT("     10%     20%     30%     40%     50%     60%     70%     80%     90% "/ &
                 9f8.2)
      WRITE(*, *)
      DEALLOCATE (x)

    CASE ('2 ')
      WRITE(*, *) 'Distribution:  Gamma'
      WRITE(*, *) 'Enter +ve value for shape parameter: '
      READ(*, *) shape
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      first = .true.
      DO i = 1, n
        x(i) = random_gamma(shape, first)
        first = .false.
      END DO
      CALL statistics(x, n, average, middle, stdev=stdvn)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Median = ', middle, '  Std.devn. = ', stdvn
      pop_mean = shape
      pop_sd = SQRT(shape)
      WRITE(*, '(1x, a, 2g14.6)') 'Population mean & st.devn. = ', pop_mean,  &
                                   pop_sd
      WRITE(*, *)

    CASE ('3 ')
      WRITE(*, *) 'Distribution:  Exponential'
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      DO i = 1, n
        x(i) = random_exponential()
      END DO
      CALL statistics(x, n, average, stdev=stdvn, percentile=pcntile)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Std.devn. = ', stdvn
      pop_mean = one
      pop_sd = one
      WRITE(*, '(1x, a, 2g14.6)') 'Population mean & st.devn. = ', pop_mean,  &
                                   pop_sd
      WRITE(*, 900) pcntile
      WRITE(*, *)

    CASE ('4 ')
      WRITE(*, *) 'Distribution:  Beta'
      WRITE(*, *) 'Enter values of a, b: '
      READ(*, *) a, b
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      first = .true.
      DO i = 1, n
        x(i) = random_beta(a, b, first)
        first = .false.
      END DO
      CALL statistics(x, n, average, stdev=stdvn, percentile=pcntile)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Std.devn. = ', stdvn
      pop_mean = a/(a+b)
      pop_sd = SQRT(a*b/(a+b+1)) / (a+b)
      WRITE(*, '(1x, a, 2g14.6)') 'Population mean & st.devn. = ', pop_mean,  &
                                   pop_sd
      WRITE(*, 900) pcntile
      WRITE(*, *)

    CASE ('5 ')
      WRITE(*, *) 'Distribution:  t'
      WRITE(*, *) 'How many degrees of freedom?: '
      READ(*, *) ndf
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      DO i = 1, n
        x(i) = random_t(ndf)
      END DO
      CALL statistics(x, n, average, stdev=stdvn, percentile=pcntile)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Std.devn. = ', stdvn
      IF (ndf .GT. 2) THEN
        pop_mean = zero
        pop_sd = SQRT(ndf / REAL(ndf-2))
        WRITE(*, '(a, 2g14.6)') ' Population mean & st.devn. = ', pop_mean, &
                                pop_sd
      END IF
      WRITE(*, 900) pcntile
      WRITE(*, *)

    CASE ('6 ')
      WRITE(*, *) 'Distribution:  Multivariate normal'
      WRITE(*, *) 'Enter no. of variables: '
      READ(*, *) k
      n = k*(k+1)/2
      ALLOCATE ( x(k), xmean(k), cov(n), chol_f(n) )
      WRITE(*, *) 'Enter means of the variables: '
      READ(*, *) xmean
      WRITE(*, *) 'Enter the covariance matrix, row-by-row'
      pos2 = 0
      DO i = 1, k
        pos1 = pos2 + 1
        pos2 = pos2 + i
        WRITE(*, '(1x, a, i2, a, i2, a, i2, a)') 'Row ', i, ' Columns ', 1,   &
                                                 ' - ', i, ' : '
        READ(*, *) cov(pos1:pos2)
      END DO
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      first = .true.
      DO i = 1, n
        CALL random_mvnorm(k, xmean, cov, chol_f, first, x, ier)
        IF (ier .NE. 0) THEN
          WRITE(*, *) '** Covariance matrix is not +ve definite **'
          EXIT
        END IF
        first = .false.
        WRITE(*, '(1x, i5, 10f7.3/ (6x, 10f7.3))') i, x
      END DO
      DEALLOCATE (x, xmean, cov, chol_f)
      WRITE(*, *)

    CASE ('7 ')
      WRITE(*, *) 'Distribution:  Generalized inverse Gaussian'
      WRITE(*, *) 'Enter values of h, b: '
      READ(*, *) a, b
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      first = .true.
      DO i = 1, n
        x(i) = random_inv_gauss(a, b, first)
        first = .false.
      END DO
      CALL statistics(x, n, average, stdev=stdvn, percentile=pcntile)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Std.devn. = ', stdvn
      WRITE(*, 900) pcntile
      WRITE(*, *)

    CASE ('8 ')
      WRITE(*, *) 'Distribution:  Poisson'
      WRITE(*, *) 'Enter +ve value for mean: '
      READ(*, *) mu
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      first = .true.
      DO i = 1, n
        x(i) = random_Poisson(mu, first)
        first = .false.
      END DO
      CALL statistics(x, n, average, middle, xmax, xmin, stdvn, pcntile)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Std.devn. = ', stdvn
      pop_mean = mu
      pop_sd = SQRT(mu)
      WRITE(*, '(1x, a, 2g14.6)') 'Population mean & st.devn. = ', pop_mean,  &
                                   pop_sd
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4)')                                 &
               'Smallest x = ', xmin, 'Largest x = ', xmax
      WRITE(*, 900) pcntile
      WRITE(*, *)

    CASE ('9 ')
      WRITE(*, *) 'Distribution:  Binomial'
      DO
        WRITE(*, *) 'Enter 1 or 2: '
        READ(*, *) which
        IF (which == 1 .OR. which == 2) EXIT
      END DO
      WRITE(*, *) 'Enter parameters N & p: '
      READ(*, *) n_binom, p
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      first = .true.
      IF (which == 1) THEN
        DO i = 1, n
          x(i) = random_binomial1(n_binom, p, first)
          first = .false.
        END DO
      ELSE
        DO i = 1, n
          x(i) = random_binomial2(n_binom, p, first)
          first = .false.
        END DO
      END IF
      CALL statistics(x, n, average, middle, xmax, xmin, stdvn, pcntile)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Std.devn. = ', stdvn
      pop_mean = n_binom * p
      pop_sd = SQRT(pop_mean * (one-p))
      WRITE(*, '(1x, a, 2g14.6)') 'Population mean & st.devn. = ', pop_mean,  &
                                   pop_sd
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4)')                                 &
               'Smallest x = ', xmin, 'Largest x = ', xmax
      WRITE(*, 900) pcntile
      WRITE(*, *)

    CASE ('10')
      WRITE(*, *) 'Distribution:  Negative binomial'
      WRITE(*, *) 'Enter parameters k & p: '
      READ(*, *) sk, p
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      first = .true.
      DO i = 1, n
        x(i) = random_neg_binomial(sk, p)
        first = .false.
      END DO
      CALL statistics(x, n, average, middle, xmax, xmin, stdvn, pcntile)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Std.devn. = ', stdvn
      pop_mean = sk * p / (one-p)
      pop_sd = SQRT(pop_mean / (one-p))
      WRITE(*, '(1x, a, 2g14.6)') 'Population mean & st.devn. = ', pop_mean,  &
                                   pop_sd
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4)')                                 &
               'Smallest x = ', xmin, 'Largest x = ', xmax
      WRITE(*, 900) pcntile
      WRITE(*, *)

    CASE ('11')
      WRITE(*, *) 'Distribution:  von Mises'
      WRITE(*, *) 'Enter +ve value for scale parameter: '
      READ(*, *) sk
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      first = .true.
      DO i = 1, n
        x(i) = random_von_Mises(sk, first)
        first = .false.
      END DO
      CALL statistics(x, n, average, middle, xmax, xmin, stdvn, pcntile)
      DEALLOCATE (x)
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4, 2x, a, f9.4)')                    &
               'Mean = ', average, 'Median = ', middle, 'Sample st.dev. = ', stdvn
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4)')                                 &
               'Smallest x = ', xmin, 'Largest x = ', xmax
      WRITE(*, 900) pcntile
      WRITE(*, *)

    CASE ('12')
      WRITE(*, *) 'Distribution: Cauchy'
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE ( x(n) )
      DO i = 1, n
        x(i) = random_Cauchy()
      END DO
      CALL statistics(x, n, average, middle, xmax, xmin, percentile=pcntile)
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4)')                                 &
               'Mean = ', average, 'Median = ', middle
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4)')                                 &
               'Smallest x = ', xmin, 'Largest x = ', xmax
      WRITE(*, 900) pcntile
      WRITE(*, *)
      DEALLOCATE ( x )

    CASE ('13')
      WRITE(*, *) 'Distribution:  Weibull'
      WRITE(*, *) 'Enter +ve value for scale parameter: '
      READ(*, *) sk
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE (x(n))
      DO i = 1, n
        x(i) = random_Weibull(sk)
      END DO
      CALL statistics(x, n, average, middle, stdev=stdvn)
      DEALLOCATE (x)
      WRITE(*, *)'Mean = ', average, '  Median = ', middle, '  Std.devn. = ', stdvn
      pop_mean = EXP( lngamma( DBLE(one/sk + one) ) )
      pop_sd = SQRT( lngamma( DBLE(two/sk + one) - pop_mean**2) )
      WRITE(*, '(1x, a, 2g14.6)') 'Population mean & st.devn. = ', pop_mean,  &
                                   pop_sd
      WRITE(*, *)

    CASE ('14')
      WRITE(*, *) 'Distribution: Chi-squared'
      WRITE(*, *) 'Enter number of degrees of freedom (integer): '
      READ(*, *) ndf
      WRITE(*, *) 'How many random deviates?: '
      READ(*, *) n
      ALLOCATE ( x(n) )
      first = .true.
      DO i = 1, n
        x(i) = random_chisq(ndf, first)
        first = .false.
      END DO
      CALL statistics(x, n, average, middle, xmax, xmin, stdvn, pcntile)
      WRITE(*, *)'Mean = ', average, '  Median = ', middle, '  Std.devn. = ', stdvn
      pop_mean = ndf
      pop_sd = SQRT(two * ndf)
      WRITE(*, '(1x, a, 2g14.6)') 'Population mean & st.devn. = ', pop_mean,  &
                                   pop_sd
      WRITE(*, '(1x, a, f9.4, 2x, a, f9.4)')                                 &
               'Smallest x = ', xmin, 'Largest x = ', xmax
      WRITE(*, 900) pcntile
      WRITE(*, *)
      DEALLOCATE ( x )

    CASE ('15')
      WRITE(*, *) 'Random order of integers 1 .. N'
      WRITE(*, *) 'Enter N: '
      READ(*, *) n
      ALLOCATE ( ix(n) )
      CALL random_order(ix, n)
      WRITE(*, '(1x, 15I5)') ix
      WRITE(*, *)
      DEALLOCATE ( ix )

    CASE DEFAULT
      WRITE(*, *) 'Option ', option, ' not allowed'
  END SELECT

END DO

STOP
END PROGRAM t_random



SUBROUTINE statistics(x, n, mean, median, xmax, xmin, stdev, percentile)
!     Calculate the requested statistics for the sample x(1), ..., x(n).
!     The percentiles estimated are 10%  , 20%  , ..., 90%.

IMPLICIT NONE

INTEGER, INTENT(IN)         :: n
REAL, INTENT(IN)            :: x(n)
REAL, OPTIONAL, INTENT(OUT) :: mean, median, xmax, xmin, stdev, percentile(9)

INTERFACE
  SUBROUTINE qsort(a, n, t)
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: n
    REAL, INTENT(IN OUT)    :: a(:)
    INTEGER, INTENT(IN OUT) :: t(:)
  END SUBROUTINE qsort
END INTERFACE

!     Local variables

INTEGER, ALLOCATABLE        :: original_pos(:)
REAL, ALLOCATABLE           :: xx(:)
INTEGER                     :: i, j
REAL                        :: pos, remainder, half = 0.5, one = 1.0

!     n must be > 0

IF (n <= 0) THEN
  WRITE(*, *)'N must be > 0, value of n input to subroutine statistics = ', n
  RETURN
END IF

!     If requested, calculate the mean & standard deviation

IF (PRESENT(mean) .OR. PRESENT(stdev)) THEN
  mean = SUM(x)/n
  IF(PRESENT(stdev)) THEN
    stdev = SUM( (x-mean)**2 )
    IF (n > 1) stdev = SQRT(stdev/(n-1))
  END IF
END IF

!     If median or percentiles have been requested, copy x to xx and sort.

IF (PRESENT(median) .OR. PRESENT(percentile)) THEN
  ALLOCATE (xx(n), original_pos(n))
  xx = x
  DO i = 1, n
    original_pos(i) = i
  END DO
  CALL qsort(xx, n, original_pos)

  IF (PRESENT(median)) THEN
    IF (n == 2*(n/2)) THEN
      i = n/2
      median = half * (xx(i) + xx(i+1))
    ELSE
      i = (n+1)/2
      median = xx(i)
    END IF
  END IF

  IF (PRESENT(percentile)) THEN
    DO i = 1, 9
      pos = i * REAL(n+1) / 10.
      j = pos
      remainder = pos - j
      percentile(i) = (one - remainder)*xx(j) + remainder*xx(j+1)
    END DO
  END IF
END IF

IF (PRESENT(xmin) .OR. PRESENT(xmax)) THEN
  xmin = MINVAL(x)
  xmax = MAXVAL(x)
END IF

RETURN
END SUBROUTINE statistics



SUBROUTINE qsort(a, n, t)

!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N. WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!     ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
!     30 stacks are enough for about 150,000 < n < 200,000

! Corrected by Jaoa Martins, Oeiras, Portugal
! Latest revision - 4 January 1999

IMPLICIT NONE
INTEGER, INTENT(IN)     :: n
REAL, INTENT(IN OUT)    :: a(:)
INTEGER, INTENT(IN OUT) :: t(:)

!     Local Variables

INTEGER   :: k, low, high, low0, high0, s, stackl(30), stackr(30), ww
REAL      :: w, x

    stackl(1) = 1
    stackr(1) = n
    s = 1
    low0 = 1
    high0 = n

10  low = low0
    high = high0
    k = (low + high)/2
    x = a(k)

    DO WHILE (low < high)
      DO WHILE (a(low) < x)
        low = low + 1
      END DO
      DO WHILE (a(high) > x)
        high = high - 1
      END DO

      IF (low <= high) THEN
        w = a(low)
        ww = t(low)
        a(low) = a(high)
        t(low) = t(high)
        a(high) = w
        t(high) = ww
      END IF
    END DO

    IF (low0 < high) THEN
      stackl(s) = high
      stackr(s) = high0
      high0 = stackl(s) - 1
      low0 = low0 - 1
      IF (low0 <= 0) low0 = 1
      s = s + 1
      GO TO 10
    END IF

    IF (low < high0) THEN
      low0 = low + 1
      GO TO 10
    END IF

    IF (s > 1) THEN
      s = s - 1
      low0 = stackl(s)
      high0 = stackr(s)
      GO TO 10
    END IF

RETURN
END SUBROUTINE qsort
