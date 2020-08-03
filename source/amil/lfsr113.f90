MODULE Lin_Feedback_Shift_Reg
! L'Ecuyer's 1999 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! N.B. This version is compatible with Lahey's ELF90
! This version requires that the default integer type is of 32-bits
! http://www.ozemail.com.au/~milleraj
! http://users.bigpond.net.au/amiller/
! Latest revision - 12 January 2001

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)
! These are unsigned integers in the C version
INTEGER, SAVE      :: s1 = 153587801, s2 = -759022222, s3 = 1288503317, &
                      s4 = -1718083407

CONTAINS

SUBROUTINE init_seeds(i1, i2, i3, i4)
IMPLICIT NONE

INTEGER, INTENT(IN) :: i1, i2, i3, i4

s1 = i1
s2 = i2
s3 = i3
s4 = i4
IF (IAND(s1,  -2) == 0) s1 = i1 - 1023
IF (IAND(s2,  -8) == 0) s2 = i2 - 1023
IF (IAND(s3, -16) == 0) s3 = i3 - 1023
IF (IAND(s4,-128) == 0) s4 = i4 - 1023

RETURN
END SUBROUTINE init_seeds



FUNCTION lfsr113() RESULT(random_numb)
! Generates a random number between 0 and 1.  Translated from C function in:
! Reference:
! L'Ecuyer, P. (1999) `Tables of maximally equidistributed combined LFSR
! generators', Math. of Comput., 68, 261-269.

! The cycle length is claimed to be about 2^(113) or about 10^(34).
! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).(2^25 - 1)

IMPLICIT NONE
REAL (dp) :: random_numb

INTEGER   :: b

! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
!      to the left if j > 0, otherwise to the right.

b  = ISHFT( IEOR( ISHFT(s1,6), s1), -13)
s1 = IEOR( ISHFT( IAND(s1,-2), 18), b)
b  = ISHFT( IEOR( ISHFT(s2,2), s2), -27)
s2 = IEOR( ISHFT( IAND(s2,-8), 2), b)
b  = ISHFT( IEOR( ISHFT(s3,13), s3), -21)
s3 = IEOR( ISHFT( IAND(s3,-16), 7), b)
b  = ISHFT( IEOR( ISHFT(s4,3), s4), -12)
s4 = IEOR( ISHFT( IAND(s4,-128), 13), b)

! The constant below is the reciprocal of (2^32 - 1)
random_numb = IEOR( IEOR( IEOR(s1,s2), s3), s4) * 2.3283064365E-10_dp + 0.5_dp

RETURN
END FUNCTION lfsr113

END MODULE Lin_Feedback_Shift_Reg



PROGRAM t_lfsr113
USE Lin_Feedback_Shift_Reg
IMPLICIT NONE
INTEGER              :: i1, i2, i3, i4, i, j, k, l, m, freq(0:15,0:15,0:15,0:15)
INTEGER, ALLOCATABLE :: seed(:)
REAL (dp)            :: x, chi_sq, expctd, deg_freedom, stdev, lower, upper

WRITE(*, *)'Enter 4 integers as seeds: '
READ(*, *) i1, i2, i3, i4
CALL init_seeds(i1, i2, i3, i4)

x = lfsr113()
i = 16*x
x = lfsr113()
j = 16*x
x = lfsr113()
k = 16*x
freq = 0
DO m= 1, 10000000
  x = lfsr113()
  l = 16*x
  freq(i,j,k,l) = freq(i,j,k,l) + 1
  i = j
  j = k
  k = l
END DO

expctd = REAL( SUM(freq) ) / (16. * 16. * 16. * 16.)
chi_sq = 0.0_dp
DO i = 0, 15
  DO j = 0, 15
    DO k = 0, 15
      chi_sq = chi_sq + SUM( (REAL(freq(i,j,k,:)) - expctd)**2 ) / expctd
    END DO
  END DO
END DO

deg_freedom = (16. * 16. * 16. * 16.) - 1.
WRITE(*, '(a, f10.1, a, f6.0, a)') ' Chi-squared = ', chi_sq, ' with ',  &
                                  deg_freedom, ' deg. of freedom'

! Now repeat the exercise with the compiler's random number generator.

CALL RANDOM_SEED(size=k)
ALLOCATE( seed(k) )
IF (i1 /= 0) THEN
  seed(1) = i1
ELSE
  seed(1) = 1234567
END IF
IF (k >= 2) seed(2) = i2
IF (k >= 3) seed(3) = i3
DO i = 4, k
  seed(i) = seed(i-3)
END DO
CALL RANDOM_SEED(put=seed)

CALL RANDOM_NUMBER(x)
i = 16*x
CALL RANDOM_NUMBER(x)
j = 16*x
CALL RANDOM_NUMBER(x)
k = 16*x
freq = 0
DO m= 1, 10000000
  CALL RANDOM_NUMBER(x)
  l = 16*x
  freq(i,j,k,l) = freq(i,j,k,l) + 1
  i = j
  j = k
  k = l
END DO

expctd = REAL( SUM(freq) ) / (16. * 16. * 16. * 16.)
chi_sq = 0.0_dp
DO i = 0, 15
  DO j = 0, 15
    DO k = 0, 15
      chi_sq = chi_sq + SUM( (REAL(freq(i,j,k,:)) - expctd)**2 ) / expctd
    END DO
  END DO
END DO

WRITE(*, '(a, f10.1, a, f6.0, a)') ' Chi-squared = ', chi_sq, ' with ',  &
                                  deg_freedom, ' deg. of freedom'

! Calculate rough limits for chi-squared based upon a normal approximation.
stdev = SQRT(2. * deg_freedom)
upper = deg_freedom + 2.*stdev
lower = deg_freedom - 2.*stdev
WRITE(*, '(a, f8.1, a, f8.1)') ' Approx. 95% limits for chi-squared: ',  &
            lower, ' - ', upper

STOP
END PROGRAM t_lfsr113
