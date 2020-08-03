PROGRAM print_constants
! Print some of the constants in the quad. precision package in decimal

USE quadruple_precision
IMPLICIT NONE

CHARACTER (LEN=36)      :: string
INTEGER                 :: ier
TYPE (quad), PARAMETER  :: pt1 = quad( 0.1_dp, -0.5551115123125784E-17_dp)
REAL (dp)               :: error

CALL quad_string(pi, string, ier)
WRITE(*, '(a, t15, a)') ' PI = ', string

CALL quad_string(rad2deg, string, ier)
WRITE(*, '(a, t15, a)') ' rad2deg = ', string

CALL quad_string(ln10, string, ier)
WRITE(*, '(a, t15, a)') ' ln10 = ', string

CALL quad_string(ln2, string, ier)
WRITE(*, '(a, t15, a)') ' ln2 = ', string

CALL quad_string(log10e, string, ier)
WRITE(*, '(a, t15, a)') ' log10e = ', string

CALL quad_string(e, string, ier)
WRITE(*, '(a, t15, a)') ' e = ', string

CALL quad_string(euler, string, ier)
WRITE(*, '(a, t15, a)') ' euler = ', string

CALL quad_string(sqrt2, string, ier)
WRITE(*, '(a, t15, a)') ' sqrt2 = ', string

CALL quad_string(sqrt10, string, ier)
WRITE(*, '(a, t15, a)') ' sqrt10 = ', string
WRITE(*, *)

! And for anyone who still thinks the constants are wrong:

CALL quad_string(pt1, string, ier)
WRITE(*, '(a, t15, a)') ' pt1 = ', string
error = pt1 - 1.0_dp / quad( 10.0_dp, 0.0_dp )
WRITE(*, '(a, g13.5)') ' Error = ', error

STOP
END PROGRAM print_constants
