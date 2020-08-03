PROGRAM demo

!     Program to demonstrate the use of module LSQ for unconstrained
!     linear least-squares regression.
!     The data on fuel consumption are from:
!     Sanford Weisberg `Applied Linear Regression', 2nd edition, 1985,
!     pages 35-36.   Publisher: Wiley

!     The program is designed so that users can easily edit it for their
!     problems.   If you are likely to have more than 500 cases then increase
!     the value of maxcases below.   Similarly, if you may have more than
!     30 predictor variables, change maxvar below.

!     The subscript 0 is used to relate to the constant in the model.
!     Module LSQ treats the constant in the model, if there is one, just as
!     any other variable.   The numbering of all arrays starts from 1 in LSQ.

USE lsq                      ! Must be the first declaration in the program
IMPLICIT NONE

INTEGER, PARAMETER :: maxcases = 500, maxvar = 30,                            &
                      max_cdim = maxvar*(maxvar+1)/2, lout = 6
INTEGER            :: case, nvar, iostatus, nreq, ifault, i, list(maxvar), j, &
                      in, i1, i2
REAL (dp)          :: xx(0:maxvar), yy, wt, one = 1.D0, beta(0:maxvar),       &
                      var, covmat(max_cdim), sterr(0:maxvar), hii,            &
                      cormat(max_cdim), ycorr(maxvar)
REAL (KIND(0.E0))  :: x(maxcases, maxvar), y(maxcases), t(0:maxvar), tmin,    &
                      r2, resid(maxcases), std_resid(maxcases), std_err_pred, &
                      fitted, stdev_res
CHARACTER (LEN=80) :: heading, output
CHARACTER (LEN= 2) :: state(50)
CHARACTER (LEN= 8) :: vname(0:maxvar), y_name
LOGICAL            :: fit_const, lindep(0:maxvar)

!     Interfaces

INTERFACE
  SUBROUTINE separate_text(heading, vname, nvar)
    IMPLICIT NONE
    CHARACTER (LEN = *), INTENT(IN OUT)            :: heading
    CHARACTER (LEN = *), DIMENSION(:), INTENT(OUT) :: vname
    INTEGER, INTENT(OUT)                           :: nvar
  END SUBROUTINE separate_text
END INTERFACE

INTERFACE
  SUBROUTINE printc(in, cormat, dimc, ycorr, vname, yname, iopt, lout, ier)
    USE lsq
    IMPLICIT NONE
    INTEGER, INTENT(IN)                          :: in, dimc, iopt, lout
    INTEGER, INTENT(OUT)                         :: ier
    REAL (dp), DIMENSION(:), INTENT(IN)          :: cormat, ycorr
    CHARACTER (LEN=*), DIMENSION(0:), INTENT(IN) :: vname
    CHARACTER (LEN=*), INTENT(IN)                :: yname
  END SUBROUTINE printc
END INTERFACE


WRITE(*, *)'The data on fuel consumption are from:'
WRITE(*, *)'Sanford Weisberg "Applied Linear Regression", 2nd edition, 1985,'
WRITE(*, *)'pages 35-36.   Publisher: Wiley   ISBN: 0-471-87957-6'
WRITE(*, *)

OPEN(8, file='fuelcons.dat', status='old')

!     The first line of this file contains the names of the variables.

READ(8, '(a)') heading

vname(0) = 'Constant'

!     Read the heading to get the variable names & the number of variables

CALL separate_text(heading, vname, nvar)

nvar = nvar - 2              ! nvar is the number of variables.
                             ! 1 is subtracted as the first field in this
                             ! file contains the abbreviated state name,
                             ! & 1 is subtracted for the Y-variable.
vname(1:nvar) = vname(2:nvar+1)
y_name = vname(nvar+2)

WRITE(*, *)'No. of variables =', nvar
WRITE(*, *)'Predictor variables are:'
WRITE(*, '(" ", 9a9)') vname(1:nvar)
WRITE(*, *)'Dependent variable is: ', y_name

fit_const = .true.           ! Change to .false. if fitting a model without
                             ! a constant.

CALL startup(nvar, fit_const)          ! Initializes the QR-factorization

!     Read in the data, one line at a time, and progressively update the
!     QR-factorization.

wt = one
case = 1

DO
  READ(8, *, IOSTAT=iostatus) state(case), x(case, 1:nvar), y(case)
  IF (iostatus > 0) CYCLE              ! Error in data
  IF (iostatus < 0) EXIT               ! End of file

  xx(0) = one                          ! A one is inserted as the first
                                       ! variable if a constant is being fitted.
  xx(1:nvar) = x(case, 1:nvar)         ! New variables and transformed variables
                                       ! will often be generated here.
  yy = y(case)
  CALL includ(wt, xx, yy)
  case = case + 1
END DO

WRITE(*, *)'No. of observations =', nobs

CALL sing(lindep, ifault)              ! Checks for singularities

IF (ifault == 0) THEN
  WRITE(*, *)'QR-factorization is not singular'
ELSE
  DO i = 1, nvar
    IF (lindep(i)) THEN
      WRITE(*, *) vname(i), ' is exactly linearly related to earlier variables'
    END IF
  END DO ! i = 1, nvar
END IF ! (ifault == 0)
WRITE(*, *)

!     Show correlations (IN = 1 for `usual' correlations)

in = 1
CALL partial_corr(in, cormat, max_cdim, ycorr, ifault)
CALL printc(in, cormat, max_cdim, ycorr, vname, y_name, 1, lout, ifault)
WRITE(*, *)
WRITE(*, *)'Press ENTER to continue'
READ(*, *)

!     Weisberg only uses variables TAX, INC, ROAD and DLIC.   These are
!     currently in positions 2, 4, 5 & 7.
!     We could have set NVAR = 4 and copied only these variables from X to
!     XX, but we will show how regressions can be performed for a subset
!     of the variables in the QR-factorization.   Routine REORDR will be used
!     to re-order the variables.

list(1:4) = (/ 2, 4, 5, 7/)
CALL reordr(list, 4, 2, ifault)        ! Re-order so that the first 4 variables
                                       ! appear in positions 2,3,4 & 5.
                                       ! N.B. Though variables # 2, 4, 5 & 7
                                       !      will occupy positions 2-5, they
                                       !      may be in any order.
                                       ! The constant remains in position 1.
WRITE(*, 910) (vname(vorder(1:ncol)))
910 FORMAT(' Current order of variables:'/' ', 8a9/)

!     Calculate regression coefficients of Y against the first variable, which
!     was just a constant = 1.0, and the next 4 predictors.

CALL tolset()                          ! Calculate tolerances before calling
                                       ! subroutine regcf.
nreq = 5                               ! i.e. Const, TAX, INC, ROAD, DLIC
CALL regcf(beta, nreq, ifault)

CALL ss()                              ! Calculate residual sums of squares

!     Calculate covariance matrix of the regression coefficients & their
!     standard errors.

var = rss(nreq) / (nobs - nreq)
CALL cov(nreq, var, covmat, max_cdim, sterr, ifault)

!     Calculate t-values

t(0:nreq-1) = beta(0:nreq-1) / sterr(0:nreq-1)

!     Output regression table, residual sums of squares, and R-squared.

WRITE(*, *)
WRITE(*, *)'Variable   Regn.coeff.   Std.error  t-value   Res.sum of sq.'
DO i = 0, nreq-1
  WRITE(*, 900) vname(vorder(i+1)), beta(i), sterr(i), t(i), rss(i+1)
  900 FORMAT(' ', a8, '  ', g12.4, '  ', g11.4, ' ', f7.2, '  ', g14.6)
END DO
WRITE(*, *)

!     Output correlations of the parameter estimates

WRITE(*, *) 'Covariances of parameter estimates'
i2 = nreq
DO i = 1, nreq-1
  i1 = i2 + 1
  i2 = i2 + nreq - i
  WRITE(output, '(" ", a8)') vname(vorder(i+1))
  WRITE(output(10*i:), '(7f10.3)') covmat(i1:i2)
  WRITE(*, '(a)') output
END DO
WRITE(*, *)

!     Now delete the variable with the smallest t-value by moving it
!     to position 5 and then repeating the calculations for the constant
!     and the next 3 variables.

j = 1
tmin = ABS(t(1))
DO i = 2, nreq-1
  IF (ABS(t(i)) < tmin) THEN
    j = i
    tmin = ABS(t(i))
  END IF
END DO
j = j + 1                    ! Add 1 as the t-array started at subscript 0

WRITE(*, *)'Removing variable in position', j
CALL vmove(j, nreq, ifault)
nreq = nreq - 1
CALL regcf(beta, nreq, ifault)
CALL ss()
CALL cov(nreq, var, covmat, max_cdim, sterr, ifault)
t(0:nreq-1) = beta(0:nreq-1) / sterr(0:nreq-1)

!     Output regression table, residual sums of squares, and R-squared.

WRITE(*, *)
WRITE(*, *)'Variable   Regn.coeff.   Std.error  t-value   Res.sum of sq.'
DO i = 0, nreq-1
  WRITE(*, 900) vname(vorder(i+1)), beta(i), sterr(i), t(i), rss(i+1)
END DO
WRITE(*, *)

var = rss(nreq)/(nobs - nreq)
stdev_res = SQRT( var )
r2 = one - rss(nreq) / rss(1)          ! RSS(1) is the rss if only the constant
                                       ! is fitted; RSS(nreq) is the rss after
                                       ! fitting the requested NREQ variables.
WRITE(*, '(" R^2 =", f8.4, "     Std. devn. of residuals =", g12.4/)')      &
      r2, stdev_res
WRITE(*, *)'N.B. Some statistical packages wrongly call the standard deviation'
WRITE(*, *)'     of the residuals the standard error of prediction'
WRITE(*, *)

!     Calculate residuals, hii, standardized residuals & standard errors of
!     prediction.

WRITE(*, *)'Press ENTER to continue'
READ(*, *)

WRITE(*, *)'State     Actual    Fitted  Residual  Std.resid.  SE(prediction)'
DO i = 1, nobs
  xx(0) = one
  DO j = 1, nreq-1
    xx(j) = x(i, vorder(j+1))          ! N.B. Regression coefficient j is for
                                       !      the variable vorder(j+1)
  END DO
  fitted = DOT_PRODUCT( beta(0:nreq-1), xx(0:nreq-1) )
  resid(i) = y(i) - fitted
  CALL hdiag(xx, nreq, hii, ifault)
  std_resid(i) = resid(i) / SQRT(var*(one - hii))
!     The sqrt was omitted from the initial version of this demo in line below.
  std_err_pred = sqrt( varprd(xx, nreq) )
  WRITE(*, 920) state(i), y(i), fitted, resid(i), std_resid(i), std_err_pred
  920 FORMAT('   ', a2, '  ', 3f10.1, f9.2, '     ', f10.0)

  IF (i .EQ. 24) THEN
    WRITE(*, *)'Press ENTER to continue'
    READ(*, *)
    WRITE(*, *)'State     Actual    Fitted  Residual  Std.resid.  SE(prediction)'
  END IF
END DO ! i = 1, nobs

STOP
END PROGRAM demo



SUBROUTINE separate_text(text, name, number)

!     Takes the character string in `text' and separates it into `names'.
!     This version only allows spaces as delimiters.

IMPLICIT NONE
CHARACTER (LEN = *), INTENT(IN OUT)            :: text
CHARACTER (LEN = *), DIMENSION(:), INTENT(OUT) :: name
INTEGER, INTENT(OUT)                           :: number

!     Local variables

INTEGER  :: length, len_name, pos1, pos2, nchars

length = LEN_TRIM(text)
number = 0
IF (length == 0) RETURN

len_name = LEN(name(1))
IF (len_name < 1) RETURN

pos1 = 1
DO
  DO                                   ! Remove any leading blanks
    IF (text(pos1:pos1) == ' ') THEN
      text(pos1:) = text(pos1+1:)
      length = length - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO
  pos2 = INDEX(text(pos1:length), ' ') ! Find position of next blank
  IF (pos2 > 0) THEN
    pos2 = pos2 + pos1 - 1
  ELSE                                 ! Last name, no blank found
    pos2 = length + 1
  END IF
  number = number + 1
  nchars = MIN(len_name, pos2-pos1)
  name(number) = text(pos1:pos1+nchars-1)
  pos1 = pos2 + 1                      ! Move past the blank
  IF (pos1 > length) RETURN
END DO

RETURN
END SUBROUTINE separate_text



SUBROUTINE printc(in, cormat, dimc, ycorr, vname, yname, iopt, lout, ier)

!     Print (partial) correlations calculated using partial_corr to unit lout.
!     If IOPT = 0, print correlations with the Y-variable only.

USE lsq

IMPLICIT NONE
INTEGER, INTENT(IN)                          :: in, dimc, iopt, lout
INTEGER, INTENT(OUT)                         :: ier
REAL (dp), DIMENSION(:), INTENT(IN)          :: cormat, ycorr
CHARACTER (LEN=*), DIMENSION(0:), INTENT(IN) :: vname
CHARACTER (LEN=*), INTENT(IN)                :: yname

!     Local variables.

INTEGER            :: nrows, j1, j2, i1, i2, row, upos, tpos, last
CHARACTER (LEN=74) :: text
CHARACTER (LEN= 9) :: char1 = ' 1.0     '

!     Check validity of arguments

ier = 0
IF (in .GE. ncol) ier = 1
IF (ncol .LE. 1) ier = ier + 2
nrows = ncol - in
IF (dimc .LE. nrows*(nrows-1)/2) ier = ier + 4
IF (ier .NE. 0) RETURN

!     If iopt.NE.0 output heading

IF (iopt .EQ. 0) GO TO 30
WRITE(lout, 900)
900 FORMAT(/'     ', 'Correlation matrix')
j1 = in + 1

DO
  j2 = MIN(j1+6, ncol)
  i1 = j1 - in
  i2 = j2 - in
  WRITE(lout, 910) vname(vorder(j1:j2))
  910 FORMAT('           ', 7(a8, ' '))

!     Print correlations for rows 1 to i2, columns i1 to i2.

  DO row = 1, i2
    text = ' ' // vname(vorder(row+in))
    IF (i1 .GT. row) THEN
      upos = (row-1) * (nrows+nrows-row) /2 + (i1-row)
      last = upos + i2 - i1
      WRITE(text(12:74), '(7(F8.5, " "))') cormat(upos:last)
    ELSE
      upos = (row-1) * (nrows+nrows-row) /2 + 1
      tpos = 12 + 9*(row-i1)
      text(tpos:tpos+8) = char1
      last = upos + i2 - row - 1
      IF (row .LT. i2) WRITE(text(tpos+9:74), '(6(F8.5, " "))') cormat(upos:last)
    END IF
    WRITE(lout, '(a)') text
  END DO ! row = 1, i2

!     Move onto the next block of columns.

  j1 = j2 + 1
  IF (j1 .GT. ncol) EXIT
END DO

!     Correlations with the Y-variable.

30 WRITE(lout, 920) yname
   920 FORMAT(/'     Correlations with the dependent variable: ', a)
   j1 = in + 1
DO
  j2 = MIN(j1+7, ncol)
  WRITE(lout, 930) vname(vorder(j1:j2))
  930 FORMAT(/' ', 8(A8, ' '))
  WRITE(lout, 940) ycorr(j1:j2)
  940 FORMAT(' ', 8(F8.5, ' '))
  j1 = j2 + 1
  IF (j1 .GT. ncol) EXIT
END DO

RETURN
END SUBROUTINE printc
