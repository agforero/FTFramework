PROGRAM t_logistic
! Use the data from the file `surgical.dat'
! The printed versions of this data all have an error.
! We first test that this program detects the error.

USE Logistic_Regression
IMPLICIT NONE

INTEGER             :: ngroups = 40, n(40), s(40), count, i, pos, ndf, ier
REAL (dp)           :: x(40,3), chisq, devnce, beta(0:3), se_beta(0:3),  &
                       p(40), stdres(40), tail_prob
CHARACTER (LEN=80)  :: text
CHARACTER (LEN=40)  :: head
CHARACTER (LEN= 5)  :: age
CHARACTER (LEN= 1)  :: sex
CHARACTER (LEN= 6)  :: vname(0:3) = (/ 'Const.', 'Area  ', 'Age   ', 'Sex   ' /)

! Read the data

OPEN(UNIT=8, FILE='surgical.dat', STATUS='OLD')
count = 0
DO
  READ(8, '(a)') text
  IF (LEN_TRIM(text) == 0) CYCLE
  IF (text(1:1) == '!') THEN
    WRITE(*, '(" ",a)') TRIM(text)
  ELSE IF (text(1:1) == 'A') THEN
    head = TRIM(text)
  ELSE
    count = count + 1
    BACKSPACE(8)
    READ(8, '(f2.0, a7, a4, i9, i8)') x(count,1), age, sex, n(count), s(count)
    pos = INDEX(age, '-')
    READ(age(pos+1:pos+2), '(i2)') i
    IF (i == 0) THEN
      x(count,2) = 88.0_dp
    ELSE IF (i == 4) THEN
      x(count,2) = 2.0_dp
    ELSE
      x(count,2) = i - 4
    END IF
    IF (sex == 'M') THEN
      x(count,3) = 0
    ELSE
      x(count,3) = 1
    END IF
  END IF
  IF (count == 40) EXIT
END DO

WRITE(*, *)
CALL logistic(ngroups, x, 3, s, n, chisq, devnce, ndf, beta, se_beta, ier)
WRITE(*, *) '*** IER =', ier
WRITE(*, *)
WRITE(*, *) 'Data being corrected'
s(39) = 50
WRITE(*, *)

! Fit full model and report

CALL logistic(ngroups, x, 3, s, n, chisq, devnce, ndf, beta, se_beta, ier,  &
              fit=p, stdres=stdres)
IF (ier /= 0) THEN
  WRITE(*, *) 'Error number', ier
ELSE
  WRITE(*, '(a, f9.3, a, f9.3, a, i4, a)')  &
        ' Deviance = ', devnce, '   Chi-squared = ', chisq,   &
        ' with ', ndf, ' deg. of freedom'
  tail_prob = 1.0_dp - chi_squared(ndf, chisq)
  IF (tail_prob < 0.01_dp) THEN
    WRITE(*, *) '*** Significantly bad fit at the 1% level ***'
  END IF
  WRITE(*, '(a, g12.4/)') ' Chi-squared tail prob. = ', tail_prob
  WRITE(*, *) '        Coefficient   Std.error'
  DO i = 0, 3
    WRITE(*, '(" ", a6, "  ", 2g13.5)') vname(i), beta(i), se_beta(i)
  END DO

! Show fitted values
  WRITE(*, '(/a)') ' Area  Age  Sex  Patients  Died   Fitted  Std.Resid.'
  DO i = 1, ngroups
    WRITE(*, '(3f5.0, 2i8, 2f9.1)') x(i,1:3), n(i), s(i), n(i)*p(i), stdres(i)
  END DO
END IF

WRITE(*, *)
WRITE(*, *) 'From the standardized residuals, the model obviously does not'
WRITE(*, *) 'well for the youngest age group.'
WRITE(*, *)
WRITE(*, *) 'Leaving out the first age group'
WRITE(*, *)

! Removing cases 1, 2, 21 & 22.

x(1:18,:) = x(3:20,:)
n(1:18) = n(3:20)
s(1:18) = s(3:20)
x(19:36,:) = x(23:40,:)
n(19:36) = n(23:40)
s(19:36) = s(23:40)
ngroups = ngroups - 4

CALL logistic(ngroups, x, 3, s, n, chisq, devnce, ndf, beta, se_beta, ier,  &
              fit=p, stdres=stdres)
IF (ier /= 0) THEN
  WRITE(*, *) 'Error number', ier
ELSE
  WRITE(*, '(a, f9.3, a, f9.3, a, i4, a)')  &
        ' Deviance = ', devnce, '   Chi-squared = ', chisq,   &
        ' with ', ndf, ' deg. of freedom'
  tail_prob = 1.0_dp - chi_squared(ndf, chisq)
  IF (tail_prob < 0.01_dp) THEN
    WRITE(*, *) '*** Significantly bad fit at the 1% level ***'
  END IF
  WRITE(*, '(a, g12.4/)') ' Chi-squared tail prob. = ', tail_prob
  WRITE(*, *) '        Coefficient   Std.error'
  DO i = 0, 3
    WRITE(*, '(" ", a6, "  ", 2g13.5)') vname(i), beta(i), se_beta(i)
  END DO

! Show fitted values
  WRITE(*, '(/a)') ' Area  Age  Sex  Patients  Died   Fitted  Std.Resid.'
  DO i = 1, ngroups
    WRITE(*, '(3f5.0, 2i8, 2f9.1)') x(i,1:3), n(i), s(i), n(i)*p(i), stdres(i)
  END DO
END IF

WRITE(*, *)
WRITE(*, *) 'The model still does not fit well, though it fits better now.'
WRITE(*, *) 'There appears to be no difference between areas, but the'
WRITE(*, *) 'difference between sexes is significant'

STOP
END PROGRAM t_logistic

