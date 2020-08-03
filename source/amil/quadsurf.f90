PROGRAM quad_surf
! Fit a quadratic surface to a set of data

USE lsq
IMPLICIT NONE
INTEGER              :: iostatus, n, ier, row, col
INTEGER, PARAMETER   :: nmax = 100
REAL (dp), PARAMETER :: wt = 1.0_dp
REAL (dp)            :: rh, temp, yval, x(nmax,5), y(nmax), xrow(0:5), b(0:5), &
                        yfit(3)

OPEN(UNIT=8, FILE='dry_bulb.dat', STATUS='OLD')
n = 0
CALL startup(5, .TRUE.)
DO
  READ(8, *, IOSTAT=iostatus) rh, temp, yval
  IF (iostatus < 0) EXIT
  IF (iostatus > 0) CYCLE
  n = n + 1
  rh = rh - 70._dp
  temp = temp - 95._dp
  y(n) = yval
  x(n,1) = rh
  x(n,2) = temp
  x(n,3) = rh**2
  x(n,4) = rh*temp
  x(n,5) = temp**2
  xrow(1:5) = x(n,:)
  xrow(0) = 1.0_dp
  CALL includ(wt, xrow, y(n))
END DO

CALL regcf(b, 6, ier)
WRITE(*, *) 'Fitted quadratic surface:'
WRITE(*, '(a, 2f8.3, a, f8.4, a, f8.4, a, f8.4, a, f8.4, a)')  &
         ' Y = ', b(0), b(1), '*X1 ', b(2), '*X2 ', b(3), '*X1^2 ',  &
         b(4), '*X1*X2 ', b(5), '*X2^2'
WRITE(*, *)
WRITE(*, *) 'Fitted values:'
WRITE(*, *) '   %RH          90       95      100'
WRITE(*, *) '-------------------------------------'
DO row = 1, 4
  xrow(0) = 1.0_dp
  xrow(1) = 20.0 - 10.0*row
  xrow(3) = xrow(1)**2
  DO col = 1, 3
    xrow(2) = -10.0 + 5.0*col
    xrow(4) = xrow(1)*xrow(2)
    xrow(5) = xrow(2)**2
    yfit(col) = DOT_PRODUCT( b, xrow )
  END DO
  WRITE(*, '(f8.0, "   ", 3f9.1)') xrow(1) + 70., yfit
END DO

STOP
END PROGRAM quad_surf
