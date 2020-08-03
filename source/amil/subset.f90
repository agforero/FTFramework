PROGRAM subset

!  Interactive program for finding best-fitting subsets of variables.
!  Author:  Alan Miller
!  Retired from CSIRO Mathematical & Information Sciences, Melbourne, Australia
!  amiller @ bigpond.net.au
!  http://www.ozemail.com.au/~milleraj
!  http://users.bigpond.net.au/amiller/

!  Latest revision - 7 August 2002
!  Thanks to Paul Mather, Nottingham, UK for:
!  1. i0 undefined when form_red_file is called.
!  2. Suggesting that END FILE be replaced with CLOSE.
!  3. Several changes to form_red_file to correctly handle the case in which
!     a constant is not being fitted.
!  Thanks to Wei-Yin Loh, Univ. of Wisconsin, USA for pointing out an error
!  in a call to REORDR in the code for cross-validation.
!  Other changes:
!  1. In SUBROUTINE update, the INTENT of xmin, xmax, xmean & sxx has been
!     changed to IN OUT (previously OUT).

USE lsq
USE find_subsets

IMPLICIT NONE

CHARACTER (LEN = 40) :: fname_dat, fname_red, fname_rpt
CHARACTER (LEN = 1)  :: ans, bel = CHAR(7), yesno
CHARACTER (LEN = 10) :: string
INTEGER              :: nvar, nvar_max, first, last, ier, i, nsize, pos, j,  &
                        ypos, in, dimc, nv, i0, best_size, line1, nrepl,    &
                        search_method, criterion, iostatus, rank_deficit, ndf
INTEGER, ALLOCATABLE :: list(:), seed(:), order_copy(:)
LOGICAL, ALLOCATABLE :: lindep(:)

REAL (dp)            :: one = 1.0, fin, fout, total_sumsq, r2, msep, press, &
                        y, e, h
REAL (dp), ALLOCATABLE  :: cormat(:), ycorr(:), beta(:), x(:), xcopy(:)
CHARACTER (LEN = 8), ALLOCATABLE :: vname(:)
CHARACTER (LEN = 8)  :: yname
CHARACTER (LEN = 42) :: version = 'Subset version 1.11, date 9 May 2001'
CHARACTER (LEN = 32) :: method(4) = (/ 'Stepwise (Efroymson)            ',  &
                                       'Sequential replacement          ',  &
                                       '2-at-a-time replacement         ',  &
                                       'Best subsets (exhaustive search)' /)
CHARACTER (LEN = 14) :: crit_name(5) = (/ 'AIC (Akaike)  ', 'BIC (Bayesian)',&
                                          'Mallows Cp    ', 'Hannan-Quinn  ',&
                                          'F-ratio = 4.0 ' /)
LOGICAL              :: fit_const, lsel = .false., OK
REAL                 :: var, Cp, Cp_last, fmax, f1, f5, f10, zero = 0.0

CALL start()                 ! Read in data and form QR reduction, if necessary

WRITE(*, *)'Enter maximum size of subset to be considered: '
READ(*, *) nvar_max
IF (nvar_max > nvar) THEN
  WRITE(*, *)'Sorry, too many.   Reset to ', nvar
  nvar_max = nvar
END IF
WRITE(11, *)'Max. subset size = ', nvar_max
WRITE(*, *)'How many subsets of each size do you want recorded?: '
READ(*, *) nbest
WRITE(11, *)'No. of best subsets = ', nbest
WRITE(11, 950) (i, vname(i), i=1, nvar)

CALL init_subsets(nvar_max, fit_const)
IF (fit_const) THEN
  i0 = 1
  total_sumsq = rss(1)
ELSE
  i0 = 0
  total_sumsq = rss(1) + d(1)*rhs(1)**2
END IF

!     Display subset selection menu and ask for choice.

DO
  WRITE(*, *)
  WRITE(*, '(30x, a)') 'Subset selection menu'
  WRITE(*, *) 'C Correlations & partial correlations    F Forward selection'
  WRITE(*, *) 'E Efroymson stepwise regression          B Backward elimination'
  WRITE(*, *) 'R Replacement sequentially               X Exhaustive search'
  WRITE(*, *) '2 Two-at-a-time replacement              D Display best found'
  WRITE(*, *) 'I Force variables IN                     O Force variables OUT'
  WRITE(*, *) 'L Least-squares regression coeffs.       V Show variable names'
  WRITE(*, *) 'M Mallows Cp for best subsets            P PRESS statistic'
  WRITE(*, *) 'S Stochastic + replacement of pairs      Q Quit this menu'
  WRITE(*, *)
  WRITE(*, *) 'N.B. Exhaustive search can be very slow.'
  WRITE(*, *) 'Use E, F, B, R and/or 2 first to establish good bounds'
  WRITE(*, *)
  WRITE(*, *) 'Enter your choice (upper or lower case OK): '
  READ(*, *) ans

  SELECT CASE (ans)
    CASE ('c', 'C')          ! Correlations & partial correlations
      IF (fit_const) THEN
        in = 1
      ELSE
        in = 0
      END IF
      nv = 0
      WRITE(*, *)'Do you want partial correlations? (Y/N): '
      READ(*, *) yesno
      IF (yesno .EQ. 'y' .OR. yesno .EQ. 'Y') THEN
        WRITE(*, *)'How many variables are to be forced in?: '
        READ(*, *) nv
        IF (nv > 0) THEN
          WRITE(*, 950) (i, vname(i), i=1, nvar)
          ALLOCATE(list(nv))
          WRITE(*, *)'Enter numbers of variables to be forced in: '
          READ(*, *) list
          CALL reordr(list, nv, in+1, ier)
          in = in + nv
        END IF
      END IF
      dimc = (last-in)*(last-in-1)/2
      ALLOCATE(cormat(dimc), ycorr(last))
      CALL partial_corr(in, cormat, dimc, ycorr, ier)
      CALL print_correlations()

      IF (nv > 0) DEALLOCATE(list)
      DEALLOCATE(cormat, ycorr)

    CASE ('f', 'F')          ! Forward selection
      CALL forwrd(first, last, ier)
      WRITE(*, *) 'Forward selection:'
      WRITE(11, *) 'Forward selection:'
      IF (fit_const) THEN
        WRITE(*, 900) (i, vname(vorder(i+1)), rss(i+1), i=0,nvar_max)
        WRITE(11, 900) (i, vname(vorder(i+1)), rss(i+1), i=0,nvar_max)
        900 FORMAT('     Order  Variable  Resid.sum of sq.'/      &
                   ('       ', i3, '  ', a8, ' ', g15.7))
      ELSE
        WRITE(*, 900) (i, vname(vorder(i)), rss(i), i=1,nvar_max)
        WRITE(11, 900) (i, vname(vorder(i)), rss(i), i=1,nvar_max)
      END IF

    CASE ('e', 'E')          ! Efroymson stepwise regression
      WRITE(*, *)'Enter F-to-add value (default = 4.0): '
      READ(*, '(a)') string
      IF (LEN_TRIM(string) .EQ. 0) THEN
        fin = 4.0
      ELSE
        READ(string, *) fin
      END IF

      WRITE(*, *)'Enter F-to-remove value (default = 4.0): '
      READ(*, '(a)') string
      IF (LEN_TRIM(string) .EQ. 0) THEN
        fout = 4.0
      ELSE
        READ(string, *) fout
      END IF
      WRITE(*, 910) fin, fout
      WRITE(11, 910) fin, fout
      910 FORMAT(' Efroymson stepwise regression algorithm'/        &
                 '     F-to-add = ', f8.2, '   F-to-remove = ', f8.2)

      CALL efroym(first, last, fin, fout, nsize, ier, 11)
      WRITE(*, *) 'Size =', nsize
      WRITE(*, 920) vname(vorder(1:nsize))
      WRITE(11, 920) vname(vorder(1:nsize))
      920 FORMAT(' Selected variables'/ (' ', 8a9))
      WRITE(*, 930) rss(nsize)
      WRITE(11, 930) rss(nsize)
      930 FORMAT(' Residual sum of squares for this model = ', g15.7/)

    CASE ('b', 'B')          ! Backward elimination
      CALL bakwrd(first, last, ier)
      WRITE(*, *) 'Backward elimination:'
      WRITE(11, *) 'Backward elimination:'
      IF (fit_const) THEN
        WRITE(*, 900) (i, vname(vorder(i+1)), rss(i+1), i=0,nvar)
        WRITE(11, 900) (i, vname(vorder(i+1)), rss(i+1), i=0,nvar)
      ELSE
        WRITE(*, 900) (i, vname(vorder(i)), rss(i), i=1,nvar)
        WRITE(11, 900) (i, vname(vorder(i)), rss(i), i=1,nvar)
      END IF

    CASE ('r', 'R')          ! Replacement sequentially
      CALL seqrep(first, last, ier)
      WRITE(*, *) 'Sequential replacement:'
      WRITE(11, *) 'Sequential replacement:'
      WRITE(*, 920) vname(vorder(1:max_size))
      WRITE(11, 920) vname(vorder(1:max_size))
      WRITE(*, 930) rss(max_size)
      WRITE(11, 930) rss(max_size)

    CASE ('x', 'X')          ! Exhaustive search
      CALL xhaust(first, last, ier)
      WRITE(*, *) 'Exhaustive search:'
      WRITE(11, *) 'Exhaustive search:'
      j = max_size*(max_size-1)/2 + 1
      WRITE(*, 920) vname(lopt(j:j+max_size-1, 1))
      WRITE(11, 920) vname(lopt(j:j+max_size-1, 1))
      WRITE(*, 930) ress(max_size,1)
      WRITE(11, 930) ress(max_size,1)

    CASE ('2')               ! Two-at-a-time replacement
      CALL seq2(first, last, ier)
      WRITE(*, *) 'Two-at-a-time sequential replacement:'
      WRITE(11, *) 'Two-at-a-time sequential replacement:'
      WRITE(*, 920) vname(vorder(1:max_size))
      WRITE(11, 920) vname(vorder(1:max_size))
      WRITE(*, 930) rss(max_size)
      WRITE(11, 930) rss(max_size)

    CASE ('s', 'S')          ! Stochastic + replacement of pairs
      WRITE(*, *) 'Random start, two-at-a-time replacement:'
      WRITE(11, *) 'Random start, two-at-a-time replacement:'
      CALL set_seed()
      WRITE(*, '(a)', ADVANCE='NO') ' How many variables in subset? '
      READ(*, *) nv
      WRITE(11, *) 'No. of variables (excl. constant) = ', nv
      WRITE(*, '(a)', ADVANCE='NO') ' How many replications? '
      READ(*, *) nrepl
      WRITE(11, *) 'No. of replicates = ', nrepl
      DO i = 1, nrepl
        IF (i == 1) THEN
          WRITE(*, '(" ", a, i6)') 'Replicate number ', i
        ELSE
          WRITE(*, '("+", a, i6, a, g14.6)')  &
                'Replicate number ', i, '  last RSS =', rss(nv+i0)
        END IF
        CALL random_pick(first, last, nv+i0)
        CALL replace2(first, last, nv+i0)
        WRITE(11, '(a, i4, a, g14.6)')  &
                  'Replicate ', nrepl, '  RSS = ', rss(nv+i0)
      END DO

    CASE ('d', 'D')          ! Display best found
      DO i = first, max_size
        WRITE(*, '(" ", a, i3, a)') 'Best subsets found of', i-i0, ' variables'
        WRITE(11, '(" ", a, i3, a)') 'Best subsets found of', i-i0, ' variables'
        WRITE(*, *) '   R.S.S.     Variable numbers'
        WRITE(11, *) '   R.S.S.     Variable numbers'
        pos = (i-1)*i/2 + 1
        DO j = 1, nbest
          WRITE(*, 940) ress(i,j), lopt(pos:pos+i-1, j)
          WRITE(11, 940) ress(i,j), lopt(pos:pos+i-1, j)
          940 FORMAT(' ', g13.5, '  ', 15i4: (/ 16(' '), 15i4))
        END DO
      END DO

    CASE ('i', 'I')          ! Force variables IN
      CALL current_status()
      ALLOCATE( list(nvar) )
      CALL get_numbers(list, nv)
      ier = 0
      IF (nv > 0) THEN
                             ! Check that none of the variables to be forced IN
                             ! is amongst any currently forced OUT.
        DO i = last+1, ncol
          IF (ANY(list(1:nv) == vorder(i))) THEN
            WRITE(*, '(a, i4, a)') ' Variable ', vorder(i),  &
                                   ' is currently forced OUT'
            nv = 0
            EXIT
          END IF
        END DO
        IF (fit_const) THEN
          CALL reordr(list, nv, 2, ier)
        ELSE
          CALL reordr(list, nv, 1, ier)
        END IF
      END IF
      IF (ier > 0 .AND. ier /= 4) THEN
        WRITE(*, *)'** Error in list of variable numbers **'
        WRITE(*, *)'** Variable not found or on list twice **'
      ELSE
        first = nv + 1
        IF (fit_const) first = nv + 2
        WRITE(11, *)
        WRITE(11, '(75a1)') ('-', i=1,75)
        WRITE(11, '(a, i4)') 'No. of variables to be forced in = ', nv
        IF (nv > 0) THEN
          WRITE(11, *) 'Variables:'
          WRITE(11, '(5(i4, 1x, a8, 2x))') (list(i), vname(list(i)), i=1,nv)
        END IF
        CALL init_subsets(nvar_max, fit_const)
      END IF
      DEALLOCATE( list )

    CASE ('o', 'O')          ! Force variables OUT
      CALL current_status()
      ALLOCATE( list(nvar) )
      CALL get_numbers(list, nv)
      ier = 0
      IF (nv > 0) THEN
                             ! Check that none of the variables to be forced OUT
                             ! is amongst any currently forced IN.
        DO i = 1, first-1
          IF (ANY(list(1:nv) == vorder(i))) THEN
            WRITE(*, '(a, i4, a)') ' Variable ', vorder(i),  &
                                   ' is currently forced IN'
            nv = 0
            EXIT
          END IF
        END DO
        DO i = 1, nv                   ! Lower variable to end positions
          DO j = first, ncol-i
            IF (vorder(j) == list(i)) THEN
              CALL vmove(j, ncol+1-i, ier)
              EXIT
            END IF
          END DO
        END DO
      END IF
      IF (ier > 0 .AND. ier /= 4) THEN
        WRITE(*, *)'** Error in list of variable numbers **'
        WRITE(*, *)'** Variable not found or on list twice **'
      ELSE
        last = ncol - nv
        WRITE(11, *)
        WRITE(11, '(75a1)') ('-', i=1,75)
        WRITE(11, '(a, i4)') 'No. of variables to be forced out = ', nv
        IF (nv > 0) THEN
          WRITE(11, *) 'Variables:'
          WRITE(11, '(5(i4, 1x, a8, 2x))') (list(i), vname(list(i)), i=1,nv)
        END IF
        CALL init_subsets(nvar_max, fit_const)
      END IF
      DEALLOCATE( list )

    CASE ('l', 'L')          ! Least-squares regression coeffs. and R^2
      WRITE(*, *)'WARNING: These estimates may be seriously biassed'
      WRITE(11, *)'WARNING: These estimates may be seriously biassed'
      CALL current_status()
      ALLOCATE( list(nvar) )
      CALL get_numbers(list, nv)
      IF (nv > 0) THEN
                             ! If variables are being forced in and/or out,
                             ! this could move them out of position.
                             ! Protect against this by copying the list of
                             ! variables and restoring at the end.
        ALLOCATE( order_copy(ncol) )
        order_copy = vorder
        CALL reordr(list, nv, i0+1, ier)
        IF (ier == 0) THEN
          IF (fit_const) THEN
            ALLOCATE( beta(0:nv) )
          ELSE
            ALLOCATE( beta(1:nv) )
          END IF
          CALL regcf(beta, nv+i0, ier)
          IF (ier == 0) THEN
            WRITE(*, *)'LS regression coefficients'
            WRITE(11, *)'LS regression coefficients'
            WRITE(*, '(3(1x, a8, 1x, g13.5, " | "))')   &
                          (vname(vorder(i+i0)), beta(i), i=1-i0, nv)
            WRITE(11, '(3(a8, 1x, g13.5, " |  "))')      &
                          (vname(vorder(i+i0)), beta(i), i=1-i0, nv)
            r2 = one - rss(nv+i0) / total_sumsq
            WRITE(*, '(1x, a, f8.4)') 'R^2 for this model = ', r2
            WRITE(11, '(1x, a, f8.4)') 'R^2 for this model = ', r2
          END IF
          DEALLOCATE( beta )
        END IF
        IF (first > i0+1) CALL reordr(order_copy, first-1, 1, ier)
        IF (last < ncol) CALL reordr(order_copy, last, 1, ier)
        DEALLOCATE( order_copy )
      END IF
      DEALLOCATE( list )

    CASE ('v', 'V')          ! Show variable names
      WRITE(*, 950) (i, vname(i), i=1, nvar)
      950 FORMAT(6(' ', i3, ' ', a8))

    CASE ('m', 'M')          ! Mallows Cp for best subsets
      IF (ncol >= nobs) THEN
        WRITE(*, *) 'No degrees of freedom available for residual variance'
        CYCLE
      END IF

      var = sserr / (nobs - ncol)
      best_size = 1
      WRITE(*, *) 'Mallows Cp - with Gilmours correction'
      WRITE(11, *) 'Mallows Cp - with Gilmours correction'
      IF (fit_const) THEN
        WRITE(*, *) 'N.B. The no. of variables INCLUDES the constant term'
        WRITE(11, *) 'N.B. The no. of variables INCLUDES the constant term'
      END IF
      DO i = first, max_size
        Cp = ress(i,1)/var - nobs + 2*i - 2*(ncol-1)/DBLE(nobs-ncol-2)
        WRITE(*, '(" ", i6, "     ", f10.2)') i, Cp
        WRITE(11, '(i6, "     ", f10.2)') i, Cp
        IF (i > first) THEN
          fmax = Cp_last - Cp + 2*(nobs-ncol-1)/DBLE(nobs-ncol-2)
          IF (fmax > zero) THEN
            CALL f1max(nobs-ncol, ncol+1-i, f1, f5, f10, ier)
            IF (fmax > f1) THEN
              WRITE(*, '(20x, a)') 'Reduction significant at the 1% level'
              WRITE(11, '(20x, a)') 'Reduction significant at the 1% level'
              best_size = i
            ELSE IF (fmax > f5) THEN
              WRITE(*, '(20(" "), a)') 'Reduction significant at the 5% level'
              WRITE(11, '(20(" "), a)') 'Reduction significant at the 5% level'
              best_size = i
            ELSE IF (fmax > f10) THEN
              WRITE(*, '(20(" "), a)') 'Reduction significant at the 10% level'
              WRITE(11, '(20(" "), a)') 'Reduction significant at the 10% level'
            ELSE
              WRITE(*, '(20(" "), a)') 'Reduction is not significant'
              WRITE(11, '(20(" "), a)') 'Reduction is not significant'
            END IF
          END IF
        END IF
        Cp_last = Cp
      END DO
      pos = best_size * (best_size - 1) / 2
      WRITE(*, *) 'A good subset for prediction appears to be:'
      WRITE(*, '(5(1x, i3, 1x, a8, " |"))')  &
                   (lopt(pos+i,1), vname(lopt(pos+i,1)), i=1, best_size)
      WRITE(11, *) 'A good subset for prediction appears to be:'
      WRITE(11, '(5(1x, i3, 1x, a8, " |"))')  &
                    (lopt(pos+i,1), vname(lopt(pos+i,1)), i=1, best_size)
      WRITE(11, *)

    CASE ('p', 'P')          ! Calculate the PRESS statistic for best subsets
      WRITE(*, *) 'PRESS statistic for best-fitting subsets'
      WRITE(11, *) 'PRESS statistic for best-fitting subsets'
      WRITE(*, *) 'Size  PRESS statistic  Mean sq. error  Variable nos.'
      WRITE(11, *) 'Size  PRESS statistic  Mean sq. error  Variable nos.'
      INQUIRE(10, OPENED=OK)
      IF (.NOT. OK) OPEN(10, FILE=fname_dat, STATUS='OLD')
      ALLOCATE( x(0:nvar), beta(ncol), xcopy(ncol) )
      x(0) = one

      DO i = first, max_size
        pos = i*(i-1)/2 + 1
        CALL reordr(lopt(pos:,1), i, 1, ier)
        CALL regcf(beta, i, ier)
        REWIND (10)
        DO j = 1, line1-1              ! Skip to line1 in data file
          READ(10, *)
        END DO
        press = zero
        DO
          IF (ypos > nvar) THEN
            READ(10, *, IOSTAT=iostatus) x(1:nvar), y
          ELSE IF (ypos .EQ. 1) THEN
            READ(10, *, IOSTAT=iostatus) y, x(1:nvar)
          ELSE
            READ(10, *, IOSTAT=iostatus) x(1:ypos-1), y, x(ypos:nvar)
          END IF

          IF (iostatus > 0) CYCLE                   ! Error in data
          IF (iostatus < 0) EXIT                    ! End of file

          e = y
          DO j = 1, i
            xcopy(j) = x(vorder(j))
            e = e - beta(j) * xcopy(j)
          END DO
          CALL hdiag(xcopy, i, h, ier)
          press = press + (e/(one - h))**2
        END DO
        WRITE(*, '(i4, 4x, 2g14.6, 2x, 10i4/ (38x, 10i4))')    &
                                   i-i0, press, press/nobs, vorder(1:i)
        WRITE(11, '(i4, 4x, 2g14.6, 2x, 10i4/ (38x, 10i4))')    &
                                   i-i0, press, press/nobs, vorder(1:i)
      END DO
      DEALLOCATE( x, beta, xcopy )

    CASE ('q', 'Q')          ! Quit this menu
      EXIT

    CASE DEFAULT
      WRITE(*, *) bel, '** Unknown option - ', ans, '  Try again! **'

  END SELECT

END DO

WRITE(*, *) 'Do you want to use cross-validation? (Y/N): '
READ(*, '(a)') ans
IF (ans == 'Y' .OR. ans == 'y') THEN
  DO
    WRITE(*, *) 'Choose search method'
    WRITE(*, *) '1. Stepwise (Efroymson)      2. Sequential replacement'
    WRITE(*, *) '3. 2-at-a-time replacement   4. Best subsets (exhaustive search)'
    WRITE(*, *) 'Enter number of your choice: '
    READ(*, *) search_method
    IF (search_method < 1 .OR. search_method > 4) THEN
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  IF (search_method > 1) THEN
    DO
      WRITE(*, *) 'Choose criterion for deciding size of subset'
      WRITE(*, *) '1. AIC (Akaike)           2. BIC (Bayesian)'
      WRITE(*, *) '3. Mallows Cp             4. Hannan-Quinn'
      WRITE(*, *) '5. F-ratio = 4.0'
      WRITE(*, *) 'Enter number of your choice: '
      READ(*, *) criterion
      IF (criterion < 1 .OR. criterion > 5) THEN
        CYCLE
      ELSE
        EXIT
      END IF
    END DO
  ELSE
    criterion = 5
  END IF

  WRITE(*, *) 'How many complete sets of 10 x 10% cross-validations: '
  READ(*, *) nrepl

  CALL RANDOM_SEED(size=i)
  WRITE(*, '(1x, a, i4, a)') 'Enter ', i, ' integers as random number seed(s): '
  ALLOCATE( seed(i) )
  READ(*, *) seed

  INQUIRE(10, OPENED=OK)
  IF (.NOT. OK) OPEN(10, FILE=fname_dat, STATUS='OLD')

  WRITE(11, '(/ a)') 'Using 10% cross-validation'
  WRITE(11, '("Search method ", a, "  Stopping criterion ", a/)')   &
                      method(search_method), crit_name(criterion)
  WRITE(11, *) 'Random number seeds: ', seed
  CALL cross_validation(10, line1, ypos, fit_const, nvar, first, last,      &
                        search_method, criterion, nrepl, seed, 11, msep, ier)
END IF   ! Want cross-validation?
STOP

CONTAINS


SUBROUTINE start()

!     This is the starting routine for the SUBSETS package of programs.
!     If a QR reduction has not been formed, then it forms the
!     upper-triangular Banachiewicz factorization of the input data.
!     Free-format input is assumed, i.e. with data fields separated by
!     spaces, CR's, tabs or commas.   N.B. Some Fortran compilers will
!     not accept tabs and/or commas as delimiters.

!     Latest revision - 10 May 2001

CHARACTER (LEN=20)    :: response
REAL (dp)             :: eps
REAL (dp), PARAMETER  :: pt0001 = 0.0001
INTEGER               :: i

!     Does the user already have a QR factorization?

WRITE(*, '(a)', ADVANCE='NO') ' Do you want to use an existing QR-factorization? (Y/N): '
READ(*, '(a)') ans
DO
  SELECT CASE (ans)
    CASE ('y', 'Y')                              ! Use existing .red file
      CALL get_file_name(fname_red, fname_red, 'red')
      OPEN(9, FILE=fname_red, STATUS='OLD', ACCESS='SEQUENTIAL',          &
           FORM='UNFORMATTED')
      CALL read_red_file()
      EXIT
    CASE ('n', 'N')
      CALL get_file_name(fname_dat, fname_dat, 'dat')
      CALL form_red_file()
      EXIT
    CASE DEFAULT
      WRITE(*, *) bel, '** You must enter Y or N **'
      CYCLE
  END SELECT
END DO

!     Set tolerances and test for singularities

eps = EPSILON(pt0001) ** 0.667
WRITE(*, '(a, g12.3)') ' Default tolerance for singularity test = ', eps
WRITE(*, *) 'This may be too small'
WRITE(*, *) 'e.g. if data recorded to say 5 significant decimals, use 1.E-5'
WRITE(*, '(a)', ADVANCE='NO')  &
            ' Enter tolerance or press RETURN to accept default: '
READ(*, '(a)') response
IF (LEN_TRIM(response) /= 0) THEN
  READ(response, '(f20.0)') eps
  eps = MAX( ABS(eps), 1.E-14)
  eps = MIN( eps, pt0001)
  WRITE(11, '(a, e12.3)') 'Tolerance used in singularity testing = ', eps
END IF
CALL tolset(eps)
IF (ALLOCATED(lindep)) DEALLOCATE(lindep)
IF (fit_const) THEN
  ALLOCATE ( lindep(0:nvar) )
ELSE
  ALLOCATE ( lindep(1:nvar) )
END IF
CALL sing(lindep, ier)
rank_deficit = -ier
IF (ANY(lindep)) THEN
  WRITE(*, '(i5, a)') -ier, ' singularities detected in predictor variables'
  WRITE(*, *) 'These variables are linearly related to earlier ones:'
  WRITE(11, '(i5, a)') -ier, ' singularities detected in predictor variables'
  WRITE(11, *) 'These variables are linearly related to earlier ones:'
  DO i = 1, nvar
    IF (lindep(i)) THEN
      WRITE(*, *) vname(i)
      WRITE(11, *) vname(i)
    END IF
  END DO
  WRITE(*, *)
  WRITE(11, *)
END IF

!     Write brief data summary to the screen and to the report file.

WRITE(*, *) 'QR-factorization is in file: ', fname_red
WRITE(11, *) 'QR-factorization is in file: ', fname_red

WRITE(*, '(" Dependent variable is: ", a8)') yname
WRITE(11, '(" Dependent variable is: ", a8)') yname
WRITE(*, 900) nvar, nobs, sserr
WRITE(11, 900) nvar, nobs, sserr
900 FORMAT(' No. of variables = ', i4, '     No. of observations = ', i5,   &
           / ' Residual sum of squares = ', g13.5)
IF (fit_const) THEN
  WRITE(*, *)'All models include a constant (intercept) term'
  WRITE(11, *)'All models include a constant (intercept) term'
  WRITE(*, 960) vname(0:nvar)
  WRITE(11, 960) vname(0:nvar)
  960 FORMAT(' Variable names: '/ (8(' ', a8)))
ELSE
  WRITE(*, *)'Note: Models do not include an intercept'
  WRITE(11, *)'Note: Models do not include an intercept'
  WRITE(*, 960) vname(1:nvar)
  WRITE(11, 960) vname(1:nvar)
END IF

IF (nobs > ncol) THEN
  ndf = nobs - ncol + rank_deficit
  WRITE(*, 910) sserr / ndf, ndf
  WRITE(11, 910) sserr / ndf, ndf
  910 FORMAT(' Residual variance estimate = ', g13.5, ' with ', i5,    &
             ' deg. of freedom'/)
END IF

!     Set default values for first & last
!     first = first variable which may an be omitted from subsets
!     last  = last variable which may be included in subsets

IF (fit_const) THEN
  first = 2
ELSE
  first = 1
END IF
last = ncol

RETURN
END SUBROUTINE start



SUBROUTINE get_file_name(fname1, fname2, extn)
!     Ask for details of the file name.   Add extension, if needed.

CHARACTER (LEN = 40), INTENT(OUT) :: fname1, fname2
CHARACTER (LEN =  3), INTENT(IN)  :: extn

!     Local variables

INTEGER              :: ipos, int_value(8)
LOGICAL              :: ok
CHARACTER (LEN = 3)  :: month(12) = (/'Jan', 'Feb', 'Mar', 'Apr', 'May',   &
                        'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/)

DO
  WRITE(*, 900)
  900 FORMAT(' Name of data file = ? ')
  READ(*, *) fname1
  fname2 = fname1

!     Add extension if none has been entered,
!     detected by the lack of a '.'

  IF (INDEX(fname1, '.') .EQ. 0) THEN
    ipos = INDEX(fname1, ' ')
    fname2 = fname1(1:ipos-1) // '.' // extn
  END IF

!     Check that file exists.

  INQUIRE(FILE=fname2, EXIST=ok)
  IF (.NOT. ok) THEN
    WRITE(*, 910) bel, fname2
    910 FORMAT(' ', a, ' ** File not found - ', a, ' **')
    CYCLE
  ELSE
    EXIT
  END IF
END DO

!     Open report file and write time and date

ipos = INDEX(fname2, '.')
fname_rpt = fname2(1:ipos) // 'rpt'
INQUIRE(FILE=fname_rpt, EXIST=ok)
IF (.NOT. ok) THEN
  OPEN(11, FILE=fname_rpt)
ELSE
  DO
    WRITE(*, *) bel, 'File: ', fname_rpt, ' already exists'
    WRITE(*, *) 'Append (A) or overwrite (O)?: '
    READ(*, '(a)') ans
    SELECT CASE (ans)
      CASE ('a', 'A')        ! Append
        OPEN(11, FILE=fname_rpt, POSITION='APPEND')
        WRITE(11, '(//)')
        EXIT
      CASE ('o', 'O')        ! Overwrite
        OPEN(11, FILE=fname_rpt, STATUS='REPLACE')
        EXIT
      CASE DEFAULT
        WRITE(*, *) bel, 'Answer must be A or O -- try again!'
    END SELECT
  END DO
END IF

WRITE(11, '(a)') version
CALL DATE_AND_TIME(values = int_value)
WRITE(11, 950) int_value(5), int_value(6), int_value(7), int_value(3),    &
               month(int_value(2)), int_value(1)
950 FORMAT('Time ', i2, ':', i2.2, ':', i2.2, '     Date ', i2, ' ', a3,  &
           ' ', i4)
RETURN
END SUBROUTINE get_file_name



SUBROUTINE form_red_file()
!     Read in data from an ASCII file and form the QR-factorization

REAL (dp)               :: weight = 1.0_dp, y
REAL (dp), ALLOCATABLE  :: x(:)
CHARACTER (LEN = 79)    :: text
INTEGER                 :: ipos, i, iostatus, ier
LOGICAL                 :: ok
REAL, ALLOCATABLE       :: xmin(:), xmax(:), xmean(:), xstdev(:)
REAL                    :: ymin, ymax, ymean, ystdev

!     Display first part of file.

OPEN(10, FILE=fname_dat, STATUS='OLD')
WRITE(*, *)'Start of your data file follows'
DO i = 1, 12
  READ(10, '(a)') text
  WRITE(*, '(" ", a)') text
END DO
REWIND (10)

WRITE(*, 920, ADVANCE='NO')
920 FORMAT(' How many X-variables ? ')
READ(*, *) nvar
WRITE(*, 930, ADVANCE='NO')
930 FORMAT(' Do you want a constant in the model ? ')
READ(*, *) ans
fit_const = (ans .EQ. 'y' .OR. ans .EQ. 'Y')
IF (fit_const) THEN
  i0 = 1
ELSE
  i0 = 0
END IF

!     Allocate memory for arrays X and VNAME.

IF (fit_const) THEN
  ALLOCATE ( x(0:nvar), vname(0:nvar) )
ELSE
  ALLOCATE ( x(nvar), vname(nvar) )
END IF

!     Get position of dependant variable.

WRITE(*, '(a)', ADVANCE='NO') ' Is dependant variable at end ? (Y/N): '
READ(*, *) ans
IF (ans == 'Y' .OR. ans == 'y') THEN
  ypos = nvar + 1
ELSE
  WRITE(*, '(a)', ADVANCE='NO') ' Enter no. of position of dependant variable: '
  READ(*, *) ypos
  IF (ypos < 1) ypos = 1
  IF (ypos > nvar) ypos = nvar + 1
END IF

!     Enter variable names, read them from file, or set defaults.

IF (fit_const) vname(0) = 'Constant'
WRITE(*, *)'Are variable names in data file ? (Y/N): '
READ(*, *) ans
IF (ans == 'Y' .OR. ans == 'y') THEN
  WRITE(*, *)'Which line do names start on ? '
  READ(*, *) line1
  IF (line1 > 1) THEN
    DO i = 1, line1-1
      READ(10, *)
    END DO
  END IF
  IF (ypos > nvar) THEN
    READ(10, *) vname(1:nvar), yname
  ELSE IF (ypos == 1) THEN
    READ(10, *) yname, vname(1:nvar)
  ELSE
    READ(10, *) vname(1:ypos-1), yname, vname(ypos:nvar)
  END IF
  REWIND (10)
ELSE
  WRITE(*, *)'Do you want to name variables ? (Y/N): '
  READ(*, '(a)') ans
  IF (ans == 'Y' .OR. ans == 'y') THEN
    WRITE(*, *)'Variable names may contain up to 8 characters'
    WRITE(*, *)'Name for dependant (Y) variable = ? '
    READ(*, '(a)') yname
    DO i = 1, nvar
      WRITE(*, '(a)', ADVANCE='NO') ' Name for variable', i, ' = ? '
      READ(*, '(a)') vname(i)
    END DO
  ELSE
    DO i = 1, nvar
      IF (i < 10) THEN
        WRITE(vname(i), '("X(", i1, ")    ")') i
      ELSE IF (i < 100) THEN
        WRITE(vname(i), '("X(", i2, ")   ")') i
      ELSE IF (i < 1000) THEN
        WRITE(vname(i), '("X(", i3, ")  ")') i
      ELSE
        WRITE(vname(i), '("X(", i4, ") ")') i
      END IF
    END DO
    yname = 'Dept.var'
  END IF
END IF

WRITE(*, '(a)', ADVANCE='NO') ' Which line does the data start on ? '
READ(*, *) line1
IF (line1 > 1) THEN
  DO i = 1, line1-1
    READ(10, *)
  END DO
END IF

CALL startup(nvar, fit_const)          ! From LSQ module

!     Allocate arrays for storing descriptive statistics
ALLOCATE( xmin(nvar), xmax(nvar), xmean(nvar), xstdev(nvar) )

!     Read in data and form the upper-triangular factorization.

IF (fit_const) x(0) = one

!     Case is skipped if spurious characters are found (e.g. for missing values).

DO
  IF (ypos > nvar) THEN
    READ(10, *, IOSTAT=iostatus) x(1:nvar), y
  ELSE IF (ypos == 1) THEN
    READ(10, *, IOSTAT=iostatus) y, x(1:nvar)
  ELSE
    READ(10, *, IOSTAT=iostatus) x(1:ypos-1), y, x(ypos:nvar)
  END IF

  IF (iostatus > 0) CYCLE                   ! Error in data
  IF (iostatus < 0) EXIT                    ! End of file

! Update descriptive statistics
  DO i = 1, nvar
    CALL update(REAL(x(i)), nobs+1, xmin(i), xmax(i), xmean(i), xstdev(i))
  END DO
  CALL update(REAL(y), nobs+1, ymin, ymax, ymean, ystdev)

  CALL includ(weight, x, y)
END DO

!     Output descriptive statistics
WRITE(11, '("Variable      MinVal       MaxVal       Mean        Std.devn.   Range/stdev")')
xstdev = SQRT(xstdev/(nobs-1))
ystdev = SQRT(ystdev/(nobs-1))
DO i = 1, nvar
  WRITE(11, '(a8, 2x, 4d13.4, f9.1)') vname(i), xmin(i), xmax(i), xmean(i),  &
                                xstdev(i), (xmax(i) - xmin(i))/xstdev(i)
END DO
WRITE(11, '(a8, 2x, 4d13.4, f9.1)') yname, ymin, ymax, ymean, ystdev,   &
                                    (ymax - ymin)/ystdev
WRITE(11, *)

!     Test for singularities & form initial array (rss) of sums of squares
!     of residuals

ALLOCATE ( lindep(ncol) )
CALL sing(lindep, ier)
IF (ier .NE. 0) THEN
  WRITE(*, *)'** Rank deficiency in data = ', -ier, ' **'
  WRITE(11, *)'** Rank deficiency in data = ', -ier, ' **'
  DO i = i0+1, ncol
    IF (lindep(i)) THEN
      WRITE(*, '(1x, a, a8, a)') 'Variable', vname(vorder(i)-i0),  &
                   ' is linearly dependent on previous variables or is constant'
    END IF
  END DO
END IF
CALL ss()

!     Change extension to .red for output file.

ipos = INDEX(fname_dat, '.')
IF (ipos == 0) ipos = LEN_TRIM(fname_dat) + 1
fname_red = fname_dat(1:ipos) // 'red'

!     Find if a file with this name already exists

INQUIRE(FILE=fname_red, EXIST=ok)
IF (ok) THEN
  WRITE(*, *)'File: ', fname_red, ' already exists.  Overwrite? (Y/N): '
  READ(*, *) ans
  ok = (ans == 'Y' .OR. ans == 'y')
ELSE
  ok = .true.
END IF

!     Write R, EL, RHS & the residual sum of squares (sserr) to disk.

IF (ok) THEN
  OPEN(9, FILE=fname_red, ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
  WRITE(9) fname_dat, ypos, line1
  WRITE(9) nvar, ncol, nobs, r_dim, lsel
  IF (.NOT. fit_const) THEN
    WRITE(9) yname, vname(1:nvar)
  ELSE
    WRITE(9) yname, vname(0:nvar)
  END IF

  WRITE(9) r, d, rhs, rss, sserr
  CLOSE(9)
END IF

DEALLOCATE ( x, lindep, xmin, xmax, xmean, xstdev )

RETURN
END SUBROUTINE form_red_file



SUBROUTINE read_red_file()
!     Read R, EL, RHS & the residual sum of squares (sserr) from disk.

IMPLICIT NONE
INTEGER :: ncases

OPEN(9, FILE=fname_red, STATUS='OLD', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
REWIND (9)
READ(9) fname_dat, ypos, line1
READ(9) nvar, ncol, ncases, r_dim, lsel
fit_const = (ncol > nvar)
CALL startup(nvar, fit_const)          ! From LSQ module
nobs = ncases

!     Allocate memory for array VNAME.

ALLOCATE ( vname(0:nvar) )

IF (.NOT. fit_const) THEN
  READ(9) yname, vname(1:nvar)
ELSE
  READ(9) yname, vname(0:nvar)
END IF

READ(9) r, d, rhs, rss, sserr
CLOSE (9)

WRITE(*, *)'File: ', fname_red, ' successfully read'

RETURN
END SUBROUTINE read_red_file



SUBROUTINE print_correlations()

!     Local variables

INTEGER              :: i1, i2, row, pos1, pos2
CHARACTER (LEN = 80) :: text

WRITE(*, *)'Do you want correlations amongst the predictor variables? (Y/N): '
READ(*, *) yesno
IF (yesno == 'y' .OR. yesno == 'Y') THEN
  IF (nv == 0) THEN
    WRITE(*, *)'Correlations amongst predictors'
    WRITE(11, *)'Correlations amongst predictors'
  ELSE
    WRITE(*, *)'Partial correlations amongst predictors'
    WRITE(*, '(1x, a/ (1x, 8a9))')'Variables forced in: ', vname(list)
    WRITE(11, *)'Partial correlations amongst predictors'
    WRITE(11, '(1x, a/ (1x, 8a9))')'Variables forced in: ', vname(list)
  END IF

!     i1 & i2 are the first & last column numbers
  i2 = in + 1
  DO
    i1 = i2 + 1
    i2 = MIN(i2 + 7, last)
    WRITE(*, '(11x, 7a9)') vname(vorder(i1:i2))
    WRITE(11, '(11x, 7a9)') vname(vorder(i1:i2))
    DO row = in+1, i2-1
      text = ' '
                     ! Find position in CORMAT of 1st off-diagonal correlation
                     ! in current row.
      pos = (row-in-1)*(last-in) - (row-in)*(row-in-1)/2 + 1
      pos2 = pos + i2 - row - 1
      WRITE(text, '(1x, a8)') vname(vorder(row))
      IF (row .GE. i1) THEN
        j = 20 + 9*(row-i1)
        pos1 = pos
      ELSE
        j = 11
        pos1 = pos + i1 - row - 1
      END IF
      WRITE(text(j:), '(7f9.4)') cormat(pos1:pos2)
      WRITE(*, '(a)') text
      WRITE(11, '(a)') text
    END DO
    IF (i2 .GE. last) EXIT
  END DO
END IF

IF (nv == 0) THEN
  WRITE(*, *)'Correlations with ', yname
  WRITE(11, *)'Correlations with ', yname
ELSE
  WRITE(*, *)'Partial correlations with ', yname
  WRITE(*, '(1x, a/ (1x, 8a9))')'Variables forced in: ', vname(list)
  WRITE(11, *)'Partial correlations with ', yname
  WRITE(11, '(1x, a/ (1x, 8a9))')'Variables forced in: ', vname(list)
END IF
i2 = in
DO
  i1 = i2 + 1
  i2 = MIN(i2 + 8, last)
  WRITE(*, '(2x, 8a9)') vname(vorder(i1:i2))
  WRITE(*, '(1x, 8f9.4)') ycorr(i1:i2)
  WRITE(11, '(2x, 8a9)') vname(vorder(i1:i2))
  WRITE(11, '(1x, 8f9.4)') ycorr(i1:i2)
  IF (i2 .GE. last) EXIT
END DO

RETURN
END SUBROUTINE print_correlations



SUBROUTINE current_status()
!     Print variable names & numbers and variables currently being forced
!     in or out (if any).

!     Local variables
CHARACTER (LEN=79) :: text
INTEGER            :: i

WRITE(*, 950) (i, vname(i), i=1, nvar)
950 FORMAT(6(' ', i3, ' ', a8))

text = 'Variables currently forced IN:  '
IF (first == i0+1) THEN
  text(33:36) = 'None'
ELSE
  DO i = i0+1, first-1
    WRITE(text(4*(i-i0)+29:4*(i-i0)+32), '(i4)') vorder(i)
  END DO
END IF
WRITE(*, '(" ", a)') text

text = 'Variables currently forced OUT: '
IF (last == ncol) THEN
  text(33:36) = 'None'
ELSE
  DO i = last+1, ncol
    WRITE(text(4*(i-last)+29:4*(i-last)+32), '(i4)') vorder(i)
  END DO
END IF
WRITE(*, '(1x, a)') text

RETURN
END SUBROUTINE current_status



SUBROUTINE f1max(ndf, nf, f1, f5, f10, ier)
! Calculates approximations to the 1%, 5% & 10% points of the distribution
! of the maximum of NF values of an F-ratio with 1 d.f. for the numerator
! and NDF d.f. for the denominator.
! An approximation is used to the values given in table 2 of:
!     Gilmour, S.G. (1996) `The interpretation of Mallows's Cp-statistic',
!     The Statistician, vol.45, pp.49-56

IMPLICIT NONE
INTEGER, INTENT(IN)  :: ndf, nf
REAL, INTENT(OUT)    :: f1, f5, f10
INTEGER, INTENT(OUT) :: ier

! Local variables
REAL :: a1(6) = (/ 1.67649, 6.94330,  1.22627, 0.25319,  0.06136, -2.41097 /)
REAL :: a5(6) = (/ 1.28152, 4.93268, -0.29583, 0.28518, -0.23566, -1.60581 /)
REAL :: a10(6) = (/1.06642, 3.96276, -0.62483, 0.30228, -0.52843, -1.04499 /)

REAL :: log_nf, one = 1.0

IF (ndf < 4 .OR. nf < 1) THEN
  ier = 1
  RETURN
END IF
ier = 0

log_nf = LOG(DBLE(nf))

f1 = one + EXP(a1(1) + (a1(3)/ndf + a1(2))/ndf + a1(4)*log_nf          &
               + (a1(6)/ndf + a1(5))/nf)

f5 = one + EXP(a5(1) + (a5(3)/ndf + a5(2))/ndf + a5(4)*log_nf          &
               + (a5(6)/ndf + a5(5))/nf)

f10 = one + EXP(a10(1) + (a10(3)/ndf + a10(2))/ndf + a10(4)*log_nf     &
                + (a10(6)/ndf + a10(5))/nf)

RETURN
END SUBROUTINE f1max



SUBROUTINE update(x, n, xmin, xmax, xmean, sxx)
! Update statistics for x
! n = observation number ( = 1 for the first call)
! (sxx is the sum of squares of deviations from the mean)

REAL, INTENT(IN)      :: x
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: xmin, xmax, xmean, sxx

! Local variables
REAL  :: dev

IF (n == 1) THEN
  xmin = x
  xmax = x
  xmean = x
  sxx = 0.0
  RETURN
END IF

IF (x < xmin) THEN
  xmin = x
ELSE IF (x > xmax) THEN
  xmax = x
END IF
dev = x - xmean
xmean = xmean + dev/n
sxx = sxx + dev * (x - xmean)
RETURN
END SUBROUTINE update



SUBROUTINE calc_residuals(unit_data, unit_resid, fname_dat, nvar, vname,   &
                          ypos, nreq, beta, list, stdev, autoc1)

INTEGER, INTENT(IN)              :: unit_data, unit_resid, nvar, ypos, nreq,  &
                                    list(:)
CHARACTER (LEN = 40), INTENT(IN) :: fname_dat
CHARACTER (LEN = 8), INTENT(IN)  :: vname(0:)
REAL (dp), INTENT(IN)            :: beta(:), stdev
REAL (dp), INTENT(OUT)           :: autoc1

! Local variables
INTEGER    :: i, ier, iostatus
REAL (dp)  :: fit, lsq_resid, h, std_resid, last_resid = 0.0, x(0:nvar), &
              y, xrow(nreq), one = 1.0, zero = 0.0, sumsq

WRITE(unit_resid) 'Data file name: ', fname_dat
WRITE(unit_resid, 900) (list(i), vname(list(i)), i=1,nreq)
900 FORMAT('Variables in model:' / (5(i3, ' ', a8, ' | ')))
WRITE(unit_resid, 910) beta(1:nreq)
910 FORMAT('Regression coefficients used:' / (5g15.6))
WRITE(unit_resid) 'Actual Y  Fitted Y  LS-residual  Std-residual  Leverage'
x(0) = one
last_resid = zero
autoc1 = zero
sumsq = zero

DO
  READ(unit_data, *, IOSTAT=iostatus) x(1:ypos-1), y, x(ypos:nvar)
  IF (iostatus > 0) THEN
    CYCLE
  ELSE IF (iostatus < 0) THEN
    EXIT
  END IF
  fit = zero
  DO i = 1, nreq
    xrow(i) = x(list(i))
    fit = fit + beta(i) * xrow(i)
  END DO
  lsq_resid = y - fit
  CALL hdiag(xrow, nreq, h, ier)
  std_resid = lsq_resid / (stdev * SQRT(one - h))
  WRITE(unit_resid, 920) y, fit, lsq_resid, std_resid, h
  920 FORMAT(3g12.4, f9.2, f10.3)
  sumsq = sumsq + lsq_resid**2
  autoc1 = autoc1 + lsq_resid * last_resid
  last_resid = lsq_resid

END DO

autoc1 = autoc1 / sumsq
WRITE(*, '(a, f9.3)') 'Lag 1 auto-correlation of residuals = ', autoc1
WRITE(unit_resid, '(a, f9.3)') 'Lag 1 auto-correlation of residuals = ', autoc1

RETURN
END SUBROUTINE calc_residuals



SUBROUTINE cross_validation(unit_data, line1, ypos, fit_const, nvar, first,  &
			    last, search_method, criterion, nrepl, seed,   &
			    unit_rpt, msep, ier)
! Cross-validation routine excluding 10% at a time

! Search_method
! 1  Efroymson stepwise (F-to-enter = F-to-delete = 4.0)
! 2  Sequential replacement
! 3  Two-at-a-time sequential replacement
! 4  Best subsets (exhaustive search with branch-and-bound)

! Criterion (set to 5 for Efroymson search)
! 1  AIC
! 2  BIC
! 3  Mallows Cp
! 4  Hannan-Quinn
! 5  F-ratio = 4.0

! Other arguments:
! unit_data Unit number from which to read the data
! line1     Number of the first line of data in the input file (in case file
!           contains variable names & other header information)
! ypos      Position of the dependent variable in each line of data
! fit_const .true. if a constant is to be included in the model
! nvar      Number of predictor variables, excluding any constant
! first     Position of the first variable available for selection.   Any
!           variables in earlier positions are to be forced into all subsets.
! last      Position of the last variable available for selection.   Variables
!           to be forced out of subsets (if any) should be after position last.
! nrepl     Number of complete replications of 10 subsets of 10% of the data
! seed()    Array of random number seeds (dimension depends upon compiler)
! unit_rpt  Unit number for output of report
! msep      Average mean squared error of the predictions
! ier       Error indicator
!           = 0 if no error detected

! This version - 15 August 2002

INTEGER, INTENT(IN)      :: unit_data, line1, ypos, nvar, first, last,  &
                            search_method, nrepl, seed(:), unit_rpt
INTEGER, INTENT(IN OUT)  :: criterion
LOGICAL, INTENT(IN)      :: fit_const
INTEGER, INTENT(OUT)     :: ier
REAL (dp), INTENT(OUT)   :: msep

! Local variables
INTEGER              :: replicate, i, nobs_full, percentile, i1, i2, case,   &
                        iostatus, num_seeds, nvar_max, nsize, ipos, lout = 6, j
INTEGER, ALLOCATABLE :: order(:), seeds(:), list(:), vorder_cpy(:),  &
                        init_vorder(:)
REAL (dp)            :: sumsq_ls, total_LS, zero = 0.0_dp, minus1 = -1.0_dp, y,   &
			weight, fin = 4.0_dp, fout = 4.0_dp, variance, one = 1.0_dp, &
                        fit_LS, penalty, crit_val, min_crit_val
REAL (dp), ALLOCATABLE :: x(:), beta_LS(:)

IF (search_method == 1) criterion = 5

ALLOCATE( order(nobs), x(0:nvar), init_vorder(ncol) )
nobs_full = nobs
msep = zero
init_vorder = vorder

ALLOCATE( list(1:nvar) )
DO i = 1, nvar
  list(i) = i
END DO

CALL RANDOM_SEED(size=num_seeds)
ALLOCATE( seeds(num_seeds) )
seeds = seed(1:num_seeds)
CALL RANDOM_SEED(put=seeds)

! Return the QR factorization back to the order of variables in the data set.
DO i = 1, nvar
  CALL reordr(list, i, 2, ier)
END DO

DO replicate = 1, nrepl
! Choose a random permutation of the integers 1, 2, ..., nobs

  DO i = 1, nobs_full
    order(i) = i
  END DO

  CALL ranord(order, nobs_full)

  total_LS = zero

! ------------------- Cycle through subsets of the data ----------------

! Delete 10% of the observations using array `order'

  DO percentile = 1, 10
    IF (percentile == 1) THEN
      i1 = 1
      i2 = 0.1 * REAL(nobs_full) + 0.5
    ELSE
      i1 = i2 + 1
      IF (percentile == 10) THEN
	i2 = nobs_full
      ELSE
        i2 = 0.1 * REAL(nobs_full) * percentile + 0.5
      END IF
    END IF

    WRITE(*, '(a, i4, a, i4)')' Replicate: ', replicate,  &
                              '  Percentile: ', percentile

    sumsq_LS = zero

    REWIND(unit_data)
    IF (line1 > 1) THEN
      DO i = 1, line1-1
        READ(unit_data, *)
      END DO
    END IF
    case = 1
    weight = minus1
    x(0) = one
    DO
      IF (ypos > nvar) THEN
        READ(unit_data, *, IOSTAT=iostatus) x(1:nvar), y
      ELSE IF (ypos == 1) THEN
        READ(unit_data, *, IOSTAT=iostatus) y, x(1:nvar)
      ELSE
        READ(unit_data, *, IOSTAT=iostatus) x(1:ypos-1), y, x(ypos:nvar)
      END IF

      IF (iostatus > 0) CYCLE                   ! Error in data
      IF (iostatus < 0) EXIT                    ! End of file

      IF(ANY(case == order(i1:i2))) THEN        ! Delete case if in this 10%
        CALL includ(weight, x, y)
      END IF
      case = case + 1
    END DO

! N.B. Subroutine INCLUD increases nobs even when weights are negative.
! Calculate correct value for the present nobs.

    nobs = nobs_full - (i2 + 1 - i1)

! Find subsets which fit well

    IF (fit_const) THEN
      nvar_max = max_size - 1
    ELSE
      nvar_max = max_size
    END IF
    CALL init_subsets(nvar_max, fit_const)
                             ! Re-order the QR-factorization if variables are
                             ! to be forced in or out.
    IF (first > 1) THEN
      CALL reordr(init_vorder, first-1, 1, ier)
    END IF
    IF (last < ncol) THEN
      CALL reordr(init_vorder, last, 1, ier)
    END IF
    SELECT CASE(search_method)
      CASE(1)
        CALL efroym(first, last, fin, fout, nsize, ier, lout)
      CASE(2)
        CALL seqrep(first, last, ier)
      CASE(3)
        CALL seq2(first, last, ier)
      CASE(4)
        CALL seq2(first, last, ier)
        CALL xhaust(first, last, ier)
    END SELECT

! Pick best subset

    min_crit_val = HUGE(one)
    IF (search_method > 1) THEN
      IF (criterion == 3) variance = sserr / (nobs - ncol)
      nsize = 0
      DO i = first, max_size
        IF (criterion /= 3) variance = ress(i,1) / (nobs - i)
        CALL calc_penalty(i, first, variance, criterion, penalty)
        crit_val = ress(i,1) + penalty
        IF (nsize == 0 .OR. crit_val < min_crit_val) THEN
          min_crit_val = crit_val
          nsize = i
        END IF
      END DO
    END IF

    ipos = ((nsize-1)*nsize)/2 + 1
    CALL reordr(lopt(ipos:,1), nsize, 1, ier)

! Estimate the regression coefficients using the LS-projections

    ALLOCATE( beta_LS(1:nsize), vorder_cpy(1:nsize) )
    vorder_cpy = vorder(1:nsize)
    CALL shell(vorder_cpy, nsize)                 ! Shell sort from find_sub
    WRITE(unit_rpt, 970) percentile, vorder_cpy(1:nsize)
    970 FORMAT('Percentile no.', i3, '   Selected variables:'/(' ', 15i5))
    CALL regcf(beta_LS, nsize, ier)
!    WRITE(unit_rpt, 980) beta_LS
!    980 FORMAT('LS regression coefficients:'/(' ', 6g13.5))
    vorder_cpy = vorder(1:nsize)

! Return the order of variables to the original order, as in the data set

    DO i = 1, nvar
      CALL reordr(list, i, 2, ier)
    END DO

! Estimate the 10% omitted, and re-instate the deleted cases

    REWIND(unit_data)
    IF (line1 > 1) THEN
      DO i = 1, line1-1
        READ(unit_data, *)
      END DO
    END IF
    case = 1
    weight = one
    x(0) = one
    DO
      IF (ypos > nvar) THEN
        READ(unit_data, *, IOSTAT=iostatus) x(1:nvar), y
      ELSE IF (ypos == 1) THEN
        READ(unit_data, *, IOSTAT=iostatus) y, x(1:nvar)
      ELSE
        READ(unit_data, *, IOSTAT=iostatus) x(1:ypos-1), y, x(ypos:nvar)
      END IF

      IF (iostatus > 0) CYCLE                   ! Error in data
      IF (iostatus < 0) EXIT                    ! End of file

      IF(ANY(case == order(i1:i2))) THEN        ! Restore case if in this 10%
        fit_LS = zero
        DO i = 1, nsize
          j = vorder_cpy(i)
          fit_LS = fit_LS + beta_LS(i) * x(j)
        END DO
        sumsq_LS = sumsq_LS + (y - fit_LS)**2
        CALL includ(weight, x, y)                ! INCLUD destroys x
      END IF
      case = case + 1
    END DO
    WRITE(unit_rpt,1000) sumsq_LS
    1000 FORMAT('Sums of sq. (LS) = ', g12.4)
    total_LS = total_LS + sumsq_LS

!   CALL print_QR
    DEALLOCATE( beta_LS, vorder_cpy )

  END DO             ! percentile = 1, 10

  WRITE(unit_rpt, 1030) total_LS
  WRITE(*, 1030) total_LS
  1030 FORMAT(/' Total sum of squares (LS) = ', g13.5)
  WRITE(unit_rpt, '(/)')
  msep = msep + total_LS

END DO             ! replicate = 1, nrepl

msep = msep / (nrepl * nobs_full)
WRITE(*, 900) msep
WRITE(unit_rpt, 900) msep
900 FORMAT(' Overall mean squared error of prediction = ', g13.5)
WRITE(*, 910) SQRT(msep)
WRITE(unit_rpt, 910) SQRT(msep)
910 FORMAT(' RMS (prediction error) = ', g13.5)

DEALLOCATE( order, x, list, seeds )
STOP
END SUBROUTINE cross_validation


SUBROUTINE ranord(order, n)

! Generate a random ordering of the integers 1 ... n.

INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(OUT) :: order(:)

! Local variables

REAL      :: wk(n)
INTEGER   :: i

DO i = 1, n
  order(i) = i
END DO
CALL RANDOM_NUMBER(wk)

CALL sqsort(wk, n, order)

RETURN
END SUBROUTINE ranord




SUBROUTINE sqsort(a, n, t)

!   NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!   BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

!   SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN OUT) :: t(:)
REAL, INTENT(IN OUT)    :: a(:)

! Local variables

INTEGER :: i, j, k, l, r, s, stackl(15), stackr(15), ww
REAL    :: w, x

s = 1
stackl(1) = 1

! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

stackr(1) = n

10 l = stackl(s)
r = stackr(s)

! KEEP SPLITTING A(L), ... , A(R) UNTIL L>=R.

s = s - 1

20 i = l
j = r
k = (l + r)/2

! REPEAT UNTIL I > J.

x = a(k)

30 IF(a(i) >= x) GO TO 40
i = i + 1
GO TO 30

40 IF(x >= a(j)) GO TO 50
j = j - 1
GO TO 40

50 IF(i > j) GO TO 60
w = a(i)
ww = t(i)
a(i) = a(j)
t(i) = t(j)
a(j) = w
t(j) = ww
i = i + 1
j = j - 1
IF(i <= j) GO TO 30

60 IF(j - l < r - i) GO TO 75
IF(l >= j) GO TO 65
s = s + 1
stackl(s) = l
stackr(s) = j

65 l = i
GO TO 90

75 IF(i >= r) GO TO 80
s = s + 1
stackl(s) = i
stackr(s) = r

80 r = j

90 IF(l < r) GO TO 20

IF(s /= 0) GO TO 10

RETURN
END SUBROUTINE sqsort



SUBROUTINE calc_penalty(size1, size2, var, penalty_num, penalty_val)

! Calculate value of penalty for size of subset.
! Currently the penalties available are:
! Number Name
!    1   AIC
!    2   BIC
!    3   Mallows' Cp
!    4   Hannan-Quinn
!    5   Efroymson (F-to-delete = F-to-add = 4.0)

INTEGER, INTENT(IN)    :: size1, size2, penalty_num
REAL(dp), INTENT(IN)   :: var
REAL(dp), INTENT(OUT)  :: penalty_val

! Local variables
REAL(dp)  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, four = 4.0_dp

IF (size1 == size2) THEN
  penalty_val = zero
  RETURN
END IF

IF (penalty_num < 1 .OR. penalty_num > 5) THEN
  penalty_val = zero
  RETURN
END IF

SELECT CASE(penalty_num)
  CASE(1)
    penalty_val = rss(size1) * (one - exp(two * (size2 - size1) / nobs))
  CASE(2)
    penalty_val = rss(size1) *    &
                  (one - exp((size2 - size1) * LOG(REAL(nobs)) / nobs))
  CASE(3)
    penalty_val = two * var * (size1 - size2)
  CASE(4)
    penalty_val = rss(size1) *    &
                  (one - exp(two * (size2 - size1) * LOG(LOG(REAL(nobs))) / nobs))
  CASE(5)
    penalty_val = four * var * (size1 - size2)
END SELECT

RETURN
END SUBROUTINE calc_penalty



SUBROUTINE get_numbers(list, n)
!     Read in a list of numbers which may be separated by commas, blanks
!     or either `..' or `-'.

INTEGER, DIMENSION(:), INTENT(OUT) :: list
INTEGER, INTENT(OUT)               :: n

!     Local variables
CHARACTER (LEN=4)   :: delimiters = ' ,-.'
CHARACTER (LEN=100) :: text
INTEGER             :: nmax, length, i1, i2, iostatus, i, number
LOGICAL             :: sequence

nmax = SIZE(list)

start: DO
  WRITE(*, *) 'Enter variable numbers on one line'
  WRITE(*, *) 'e.g. 1-5 8 11 .. 15  use commas or blanks as separators'
  WRITE(*, *) ': '
  READ(*, '(a)') text
  text = ADJUSTL(text)
  length = LEN_TRIM(text)
  IF (length == 0) THEN
    n = 0
    RETURN
  END IF

  n = 1
  i1 = 1
  sequence = .FALSE.
  DO
    i2 = SCAN( text(i1:), delimiters )
    IF (i2 == 0) THEN
      i2 = length
    ELSE
      i2 = i2 + i1 - 2
    END IF
    READ(text(i1:i2), *, IOSTAT=iostatus) number
    IF (iostatus /= 0) THEN
      WRITE(*, *) '** Error: numeric data expected **'
      WRITE(*, '(1x, a)') text(1:length)
      text = ' '
      DO i = i1, i2
        text(i:i) = '^'
      END DO
      WRITE(*, '(1x, a)') text(1:i2)
      CYCLE start
    END IF

    IF (sequence) THEN
      IF (number <= list(n-1)) THEN
        WRITE(*, *) 'Variable numbers not increasing'
        WRITE(*, '(1x, a)') text(1:length)
        text = ' '
        DO i = i1, i2
          text(i:i) = '^'
        END DO
        WRITE(*, '(1x, a)') text(1:i2)
        CYCLE start
      END IF
      DO
        list(n) = list(n-1) + 1
        IF (list(n) >= number) EXIT
        n = n + 1
      END DO
    ELSE
      list(n) = number
    END IF
    IF (i2 == length) RETURN

    i1 = i2 + 1
    sequence = .FALSE.

                                       ! Find end of delimiters
    DO
      IF ( SCAN( text(i1:i1), delimiters ) > 0) THEN
        IF (text(i1:i1) == '-' .OR. text(i1:i1+1) == '..') sequence = .TRUE.
        i1 = i1 + 1
      ELSE
        EXIT
      END IF
    END DO
    n = n + 1
    IF (n > nmax) THEN
      WRITE(*, *) '** Too many numbers entered - list truncated **'
      n = nmax
      RETURN
    END IF
  END DO
END DO start

RETURN
END SUBROUTINE get_numbers



SUBROUTINE set_seed()

INTEGER, ALLOCATABLE :: seed(:)
INTEGER              :: k

!     Set the random number seed.

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
CALL RANDOM_SEED(get=seed)
WRITE(*, *)'Old random number seeds: ', seed

WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
WRITE(11, '(a/ 10(" ", i10))') 'New random number seeds:', seed
CALL RANDOM_SEED(put=seed)

RETURN
END SUBROUTINE set_seed

END PROGRAM subset
