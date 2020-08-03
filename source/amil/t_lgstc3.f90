PROGRAM t_logistic
! Illustration of the problem of a boundary between two regions.
! In one region, every case falls into one category, in the other
! region, every case is in the other category.

USE Logistic_Regression
IMPLICIT NONE

REAL (dp)  :: x(16,2), chisq, devnce, beta(0:2), se_beta(0:2), scale
INTEGER    :: n(16), s(16), i, iostatus, ndf, ier

! Read in the artificial data

OPEN(UNIT=8, FILE='clearcut.dat', STATUS='OLD')

i = 1
DO
  READ(8, *, IOSTAT=iostatus) x(i,1), x(i,2), s(i), n(i)
  IF (iostatus < 0) EXIT
  IF (iostatus > 0) CYCLE
  i = i + 1
  IF (i > 16) EXIT
END DO

CALL logistic(16, x, 2, s, n, chisq, devnce, ndf, beta, se_beta, ier)
WRITE(*, *) 'IER =', ier
WRITE(*, *)
WRITE(*, *) 'If IER = 7, it means that there is a linear boundary'
WRITE(*, *) 'between cases in which s = 0 and s = 1'

IF (ier == 0) WRITE(*, '(a, f8.3, a, f8.3, a, i3, a)')  &
                       ' Deviance = ', devnce, '  Chi-square = ', chisq,  &
                       '  with ', ndf, ' deg. of freedom'
WRITE(*, *)
WRITE(*, *) '   Coefficient'
DO i = 0, 2
  WRITE(*, '(i3, g12.4)') i, beta(i)
END DO
WRITE(*, *)

scale = SQRT(ABS(beta(1))*ABS(beta(2)))
beta = beta / scale
WRITE(*, *) 'Boundary is approx:'
WRITE(*, '(f9.4, a, f9.4, a, f9.4, a)')  &
         beta(0), ' + ', beta(1), ' * X1 + ', beta(2), ' * X2 = 0.0'

STOP
END PROGRAM t_logistic

