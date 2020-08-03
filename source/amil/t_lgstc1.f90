PROGRAM t_logistic
! Use the low infant birth weight data from Hosmer & Lemeshow

USE Logistic_Regression
IMPLICIT NONE

CHARACTER (LEN= 6)  :: vname(0:10) = (/  &
                       'Const.', 'Age   ', 'LWT   ', 'Race_2', 'Race_3', &
                       'Smoke?', 'PTL   ', 'HT    ', 'UI    ', 'FTV   ', &
                       'BWT   ' /)
INTEGER    :: iostatus, id, low, age, lwt, race, smoke, ptl, ht, ui, ftv,  &
              bwt, ncases, n(200), s(200), ndf, ier, i

REAL (dp)  :: x(200, 9), chisq, devnce, beta(0:9), se_beta(0:9), tail_prob

! Read the data file

OPEN(UNIT=8, FILE='birthwt.dat', STATUS='OLD')
ncases = 0
DO
  READ(8, *, IOSTAT=iostatus) id, low, age, lwt, race, smoke, ptl, ht, ui,  &
                              ftv, bwt
  IF (iostatus < 0) EXIT
  IF (iostatus > 0) CYCLE
  ncases = ncases + 1
  n(ncases) = 1
  s(ncases) = low
  x(ncases, 1) = age
  x(ncases, 2) = lwt
  x(ncases, 3:4) = 0.0_dp
  IF (race > 1) x(ncases, race+1) = 1.0_dp
  x(ncases, 5) = smoke
  x(ncases, 6) = ptl
  x(ncases, 7) = ht
  x(ncases, 8) = ui
  x(ncases, 9) = ftv
END DO

WRITE(*, *) 'NCASES = ', ncases
WRITE(*, *)

! Fit full model and report

CALL logistic(ncases, x, 9, s, n, chisq, devnce, ndf, beta, se_beta, ier)
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
  DO i = 0, 9
    WRITE(*, '(" ", a6, "  ", 2g13.5)') vname(i), beta(i), se_beta(i)
  END DO

END IF

! Omit the last variable for which the coefficient is much smaller
! than its standard error.

WRITE(*, *)
CALL logistic(ncases, x, 8, s, n, chisq, devnce, ndf, beta, se_beta, ier)
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
  DO i = 0, 8
    WRITE(*, '(" ", a6, "  ", 2g13.5)') vname(i), beta(i), se_beta(i)
  END DO

END IF

! Now drop variable 1 (age)

x(:,1:7) = x(:,2:8)
vname(1:7) = vname(2:8)
WRITE(*, *)
CALL logistic(ncases, x, 7, s, n, chisq, devnce, ndf, beta, se_beta, ier)
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
  DO i = 0, 7
    WRITE(*, '(" ", a6, "  ", 2g13.5)') vname(i), beta(i), se_beta(i)
  END DO

END IF

STOP
END PROGRAM t_logistic
