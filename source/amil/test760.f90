PROGRAM tprg3p
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-11  Time: 10:12:29

! Test program for the RGBI3P/RGSF3P subroutine package

! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 1995/08

! This program calls the RGBI3P and RGSF3P subroutines.

! This program requires no input data files.

! This program creates the TPRG3P data file.  All elements of
! the DZI array in the data file are expected to be zero.


! Specification statements
USE Grid_Interpolation
IMPLICIT NONE

!     .. Parameters ..
INTEGER, PARAMETER  :: newpg = 210000000
INTEGER, PARAMETER  :: nxd = 9, nyd = 11, nxi = 19, nyi = 23
REAL, PARAMETER     :: ximn = -0.5, ximx = 8.5, yimn = -0.5, yimx = 10.5
!     ..
!     .. Local Scalars ..
REAL     :: anxim1, anyim1, dxi, dyi
INTEGER  :: ier, isec, ixd, ixi, iximn, iximx, iyd, iydr, iyi, iyir, md, nydo2
CHARACTER (LEN=6)  :: nmpr = 'TPRG3P', nmwf = 'WFRG3P'
!     ..
!     .. Local Arrays ..
REAL     :: dzi(nxi,nyi), xi(nxi), yi(nyi), zi(nxi,nyi)
CHARACTER (LEN=6)   :: nmsr(2) = (/ 'RGBI3P', 'RGSF3P' /)
CHARACTER (LEN=20)  :: lbl(2) = (/ 'Calculated ZI Values',  &
                                   'Differences         ' /)
!     ..
!     .. External Subroutines ..
! EXTERNAL rgbi3p, rgsf3p
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC MOD, REAL
!     ..
! Data statements
REAL  :: xd(nxd) = (/ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 /)
REAL  :: yd(nyd) = (/ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 /)
REAL  :: zd(nxd,nyd) = RESHAPE(  &
                    (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
                       3.2, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
                        7.4,  4.8,  1.4,  0.1,  0.0, 0.0, 0.0, 0.0, 0.0,   &
                       12.0,  8.0,  5.3,  2.9,  0.6, 0.0, 0.0, 0.0, 0.0,   &
                       16.8, 14.4,  8.1,  6.9,  6.2,  0.6,  0.1, 0.0,  0.0, &
                       21.8, 20.5, 12.8, 17.6,  5.8,  7.6,  0.8,  0.6,  0.6, &
                       22.4, 22.5, 14.6, 22.5,  4.7,  7.2,  1.8,  2.1,  2.1, &
                       37.2, 40.0, 27.0, 41.3, 14.1, 24.5, 17.3, 20.2, 20.8, &
                       58.2, 61.5, 47.9, 62.3, 34.6, 45.5, 38.2, 41.2, 41.7 /), &
                       (/ 9, 11 /) )
REAL  :: zie(nxi,23) = RESHAPE(  &
 (/ -.847, -.533, -.274, -.117, -.031, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, .401, .250, .119,   &
  .043, .011, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -.665, -.376, -.143, -.033, -.007, 0.0, 0.0,  &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  &
  0.0, 2.449, 1.368, .537, .149, .025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.083, 3.200, 1.642, .700, .187, 0.0,  &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.588,   &
  5.234, 3.878, 2.542, 1.188, .253, .026, .026, .007, 0.0, 0.0, 0.0, 0.0,   &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.017, 7.400, 6.400, 4.800, 2.963, 1.400,   &
  .457, .100, .027, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11.055, &
  9.670, 8.083, 6.305, 4.786, 3.421, 2.043, 1.112, .565, .131, -.019, 0.0,  &
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 14.492, 12.000, 9.746, 8.000, 6.594,   &
  5.300, 4.081, 2.900, 1.697, .600, .059, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  &
  0.0, 15.999, 14.376, 12.657, 10.774, 8.620, 6.659, 5.291, 4.392, 3.926,   &
  3.005, 1.223, .139, .051, .025, .009, 0.0, 0.0, 0.0, -.005, 15.525, 16.800, &
  16.749, 14.400, 10.956, 8.100, 6.735, 6.900, 7.298, 6.200, 3.010, .600,   &
  .248, .100, .024, 0.0, .006, 0.0, -.025, 15.876, 19.280, 20.563, 17.856,   &
  13.242, 10.219, 10.577, 11.999, 10.170, 7.053, 5.198, 3.543, 1.831, .350,  &
  -.130, .168, .408, .168, -.224, 17.700, 21.800, 23.531, 20.500, 15.087,   &
  12.800, 15.817, 17.600, 11.477, 5.800, 6.988, 7.600, 4.410, .800, -.392,  &
  .600, 1.261, .600, -.417, 17.913, 22.788, 24.944, 21.881, 16.302, 14.382,  &
  18.557, 20.807, 11.916, 4.561, 7.327, 8.518, 5.133, 1.284, -.013, 1.201,   &
  1.998, 1.200, -.065, 16.383, 22.400, 25.330, 22.500, 16.796, 14.600,   &
  19.172, 22.500, 13.159, 4.700, 6.689, 7.200, 4.392, 1.800, 1.150, 2.100,   &
  2.734, 2.100, 1.025, 18.109, 26.756, 31.311, 28.143, 21.004, 18.237,   &
  24.236, 28.979, 17.970, 7.469, 10.467, 11.985, 9.022, 6.833, 6.901, 8.292,  &
  9.186, 8.524, 7.101, 24.667, 37.200, 44.007, 40.000, 30.508, 27.000,   &
  34.974, 41.300, 27.136, 14.100, 20.473, 24.500, 20.557, 17.300, 17.639,  &
  20.200, 21.826, 20.800, 18.458, 33.414, 48.009, 56.017, 51.561, 40.817,  &
  36.922, 45.856, 52.860, 37.376, 23.200, 30.839, 36.192, 31.969, 28.037,  &
  28.437, 31.604, 33.579, 32.332, 29.561, 44.842, 58.200, 65.537, 61.500,  &
  51.657, 47.900, 55.899, 62.300, 47.891, 34.600, 41.239, 45.500, 41.479,  &
  38.200, 38.591, 41.200, 42.823, 41.700, 39.192, 58.284, 68.917, 74.644,  &
  71.333, 63.413, 60.125, 66.293, 71.400, 59.129, 47.725, 52.451, 54.592,  &
  50.842, 48.483, 48.639, 50.142, 51.089, 50.200, 48.268 /), (/ 19, 23 /) )

!     ..
! Calculation
! Opens the output file and writes the input data.
OPEN (6,FILE=nmwf)
nydo2 = nyd / 2
WRITE (6,FMT=5000) nmpr
WRITE (6,FMT=5100) xd
DO  iydr = 1, nyd
  IF (MOD(iydr-1,nydo2) <= 1) WRITE (6,FMT='(1X)')
  iyd = nyd + 1 - iydr
  WRITE (6,FMT=5200) yd(iyd), (zd(ixd,iyd),ixd = 1,nxd)
END DO

! Program check for the RGBI3P subroutine
! - Performs interpolation and calculates the differences.
dxi = ximx - ximn
anxim1 = nxi - 1
DO  ixi = 1, nxi
  xi(ixi) = ximn + dxi*REAL(ixi-1) / anxim1
END DO
dyi = yimx - yimn
anyim1 = nyi - 1
DO  iyi = 1, nyi
  yi(iyi) = yimn + dyi*REAL(iyi-1) / anyim1
END DO
DO  iyi = 1, nyi
  DO  ixi = 1, nxi
    IF (ixi == 1.AND.iyi == 1) THEN
      md = 1
    ELSE
      md = 2
    END IF
    CALL rgbi3p(md, nxd, nyd, xd, yd, zd, 1, xi(ixi), yi(iyi), zi(ixi,iyi), ier)
    IF (ier > 0) STOP
    dzi(ixi,iyi) = zi(ixi,iyi) - zie(ixi,iyi)
  END DO
END DO

! - Writes the calculated results.
WRITE (6,FMT=5300) newpg, nmpr, nmsr(1), lbl(1)
DO  isec = 1, 2
  IF (isec == 1) THEN
    iximn = 1
    iximx = 11
  ELSE
    iximn = 9
    iximx = nxi
  END IF
  WRITE (6,FMT=5400) (xi(ixi),ixi = iximn,iximx)
  DO  iyir = 1, nyi
    iyi = nyi + 1 - iyir
    WRITE (6,FMT=5500) yi(iyi), (zi(ixi,iyi),ixi = iximn,iximx)
  END DO
END DO

! - Writes the differences.
WRITE (6,FMT=5300) newpg, nmpr, nmsr(1), lbl(2)
DO  isec = 1, 2
  IF (isec == 1) THEN
    iximn = 1
    iximx = 11
  ELSE
    iximn = 9
    iximx = nxi
  END IF
  WRITE (6,FMT=5600) (xi(ixi),ixi = iximn,iximx)
  DO  iyir = 1, nyi
    iyi = nyi + 1 - iyir
    WRITE (6,FMT=5500) yi(iyi), (dzi(ixi,iyi),ixi = iximn,iximx)
  END DO
END DO

! Program check for the RGSF3P subroutine
! - Performs surface fitting and calculates the differences.
md = 1
CALL rgsf3p(md, nxd, nyd, xd, yd, zd, nxi, xi, nyi, yi, zi, ier)
IF (ier > 0) STOP
DO  iyi = 1, nyi
  DO  ixi = 1, nxi
    dzi(ixi,iyi) = zi(ixi,iyi) - zie(ixi,iyi)
  END DO
END DO

! - Writes the calculated results.
WRITE (6,FMT=5300) newpg, nmpr, nmsr(2), lbl(1)
DO  isec = 1, 2
  IF (isec == 1) THEN
    iximn = 1
    iximx = 11
  ELSE
    iximn = 9
    iximx = nxi
  END IF
  WRITE (6,FMT=5400) (xi(ixi),ixi = iximn,iximx)
  DO  iyir = 1, nyi
    iyi = nyi + 1 - iyir
    WRITE (6,FMT=5500) yi(iyi), (zi(ixi,iyi),ixi = iximn,iximx)
  END DO
END DO

! - Writes the differences.
WRITE (6,FMT=5300) newpg, nmpr, nmsr(2), lbl(2)
DO  isec = 1, 2
  IF (isec == 1) THEN
    iximn = 1
    iximx = 11
  ELSE
    iximn = 9
    iximx = nxi
  END IF
  WRITE (6,FMT=5600) (xi(ixi),ixi = iximn,iximx)
  DO  iyir = 1, nyi
    iyi = nyi + 1 - iyir
    WRITE (6,FMT=5500) yi(iyi), (dzi(ixi,iyi),ixi = iximn,iximx)
  END DO
END DO
STOP

! Format statements
5000 FORMAT (a6, t17, 'Original Data'//// t36, 'ZD(XD,YD)')
5100 FORMAT ('    YD    XD='/ t8, f8.1, 2(' ', 3F6.1, f7.1),/)
5200 FORMAT (' ', f6.1, f8.1, 2(' ', 3F6.1, f7.1))
5300 FORMAT (a1, a6, t14, 'Program Check for ', a6, t43, a20)
5400 FORMAT (/t39, 'ZI(XI,YI)'/ '  YI   XI='/ t6, 3F7.2, 2F6.2, 2f7.2,  &
             2F6.2, 2F7.2/)
5500 FORMAT (f5.2, 3F7.2, 2F6.2, 2F7.2, 2F6.2, 2F7.2)
5600 FORMAT (/t39, 'DZI(XI,YI)'/ '  YI   XI='/ t6, 3F7.2, 2F6.2,  &
             2F7.2, 2F6.2, 2F7.2,/)
END PROGRAM tprg3p
