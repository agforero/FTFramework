MODULE quad_prec_complex
! Quadruple-precision complex arithmetic.
! This should compile and run correctly using any compiler on which the
! module quadruple_precision runs correctly.   There should be no need
! for different versions for different compilers.

! Programmer : Alan Miller, CSIRO Mathematical & Information Sciences
! e-mail: alan @ vic.cmis.csiro.au   URL:  www.ozemail.com.au/~milleraj
! Fax: (+61) 3-9545-8080

! Latest revision - 11 November 1999

USE quadruple_precision
IMPLICIT NONE

TYPE  :: qc
  TYPE (quad) :: qr, qi                ! quad_real, quad_imaginary
END TYPE qc

INTERFACE OPERATOR (+)
  MODULE PROCEDURE qc_add
END INTERFACE

INTERFACE OPERATOR (-)
  MODULE PROCEDURE qc_sub
  MODULE PROCEDURE negate_qc
END INTERFACE

INTERFACE OPERATOR (*)
  MODULE PROCEDURE qc_mult
END INTERFACE

INTERFACE OPERATOR (/)
  MODULE PROCEDURE qc_div
END INTERFACE

INTERFACE ABS
  MODULE PROCEDURE qc_abs
END INTERFACE

INTERFACE CMPLX
  MODULE PROCEDURE qc_cmplx
END INTERFACE

INTERFACE CONJG
  MODULE PROCEDURE qc_conjg
END INTERFACE

INTERFACE AIMAG
  MODULE PROCEDURE qc_aimag
END INTERFACE

INTERFACE SQRT
  MODULE PROCEDURE qc_sqrt
END INTERFACE

INTERFACE LOG
  MODULE PROCEDURE qc_log
END INTERFACE

INTERFACE EXP
  MODULE PROCEDURE qc_exp
END INTERFACE

CONTAINS


FUNCTION qc_add(x, y) RESULT(z)
! Quadruple-precision complex addition

TYPE (qc), INTENT(IN) :: x, y
TYPE (qc)             :: z

z % qr = x % qr + y % qr
z % qi = x % qi + y % qi

RETURN
END FUNCTION qc_add



FUNCTION qc_sub(x, y) RESULT(z)
! Quadruple-precision complex subtraction

TYPE (qc), INTENT(IN) :: x, y
TYPE (qc)             :: z

z%qr = x%qr - y%qr
z%qi = x%qi - y%qi

RETURN
END FUNCTION qc_sub



FUNCTION negate_qc(x) RESULT(z)
! Quadruple-precision complex negate sign

TYPE (qc), INTENT(IN) :: x
TYPE (qc)             :: z

z%qr = - x%qr
z%qi = - x%qi

RETURN
END FUNCTION negate_qc



FUNCTION qc_mult(x, y) RESULT(z)
! Quadruple-precision complex multiplication

TYPE (qc), INTENT(IN) :: x, y
TYPE (qc)             :: z

z%qr = x%qr * y%qr - x%qi * y%qi
z%qi = x%qi * y%qr + x%qr * y%qi

RETURN
END FUNCTION qc_mult



FUNCTION qc_div(x, y) RESULT(z)
! Quadruple-precision complex division

TYPE (qc), INTENT(IN) :: x, y
TYPE (qc)             :: z

TYPE (quad) :: t, d

IF ( ABS(y%qr%hi) > ABS(y%qi%hi) ) THEN
  t = y%qi / y%qr
  d = y%qr + t * y%qi
  z%qr = (x%qr + t * x%qi) / d
  z%qi = (x%qi - t * x%qr) / d
ELSE
  t = y%qr / y%qi
  d = y%qi + t * y%qr
  z%qr = (t * x%qr + x%qi) / d
  z%qi = (t * x%qi - x%qr) / d
END IF

RETURN
END FUNCTION qc_div



FUNCTION qc_cmplx(xr, xi) RESULT(z)
! Convert pair of quad.precision numbers to qc.

TYPE (quad), INTENT(IN) :: xr, xi
TYPE (qc)               :: z

z%qr = xr
z%qi = xi

RETURN
END FUNCTION qc_cmplx



FUNCTION qc_aimag(x) RESULT(z)
! Convert pair of quad.precision numbers to qc.

TYPE (qc), INTENT(IN)   :: x
TYPE (quad)             :: z

z = x%qi

RETURN
END FUNCTION qc_aimag



FUNCTION qc_conjg(x) RESULT(z)
! Quadruple-precision complex square root

TYPE (qc), INTENT(IN) :: x
TYPE (qc)             :: z

z%qr = x%qr
z%qi = - x%qi

RETURN
END FUNCTION qc_conjg



FUNCTION qc_sqrt(x) RESULT(z)
! Quadruple-precision complex square root

TYPE (qc), INTENT(IN) :: x
TYPE (qc)             :: z

! Local variables
TYPE (quad)           :: r
REAL (dp), PARAMETER  :: vsmall = TINY(1.0_dp)

r = SQRT(x%qr **2 + x%qi **2)
IF (r%hi < vsmall) THEN
  z%qr = 0.0_dp
  z%qi = 0.0_dp
  RETURN
END IF

z%qr = SQRT( SCALE( r + ABS(x%qr), -1 ) )
z%qi = ABS( SCALE( x%qi, -1 ) ) / z%qr
IF (x%qr%hi < 0.0_dp) THEN
  r = z%qr
  z%qr = z%qi
  z%qi = r
END IF
IF (x%qi%hi < 0.0_dp) z%qi = - z%qi

RETURN
END FUNCTION qc_sqrt



FUNCTION qc_abs(x) RESULT(z)
! Quadruple-precision complex absolute (modulus)
! SQRT( x%qr **2 + x%qi **2 ) avoiding overflow
! Using either:   X^2 + Y^2 = X^2.(1 + (Y/X)^2)
!           or:             = Y^2.(1 + (X/Y)^2)

TYPE (qc), INTENT(IN) :: x
TYPE (quad)           :: z

TYPE (quad) :: a

IF ( ABS(x%qr%hi) > ABS(x%qi%hi) ) THEN
  a = x%qi / x%qr
  z = ABS( x%qr ) * SQRT( quad(1.0_dp, 0.0_dp) + a**2 )
ELSE IF ( x%qi%hi == 0.0_dp ) THEN
  z = 0.0_dp
ELSE
  a = x%qr / x%qi
  z = ABS( x%qi ) * SQRT( quad(1.0_dp, 0.0_dp) + a**2 )
END IF

RETURN
END FUNCTION qc_abs



FUNCTION qc_exp(x) RESULT(z)
! Quadruple-precision complex exponential

TYPE (qc), INTENT(IN) :: x
TYPE (qc)             :: z

TYPE (quad)           :: expx

expx = EXP(x%qr)
z%qr = expx * COS(x%qi)
z%qi = expx * SIN(x%qi)

RETURN
END FUNCTION qc_exp



FUNCTION qc_log(x) RESULT(z)
! Quadruple-precision complex logarithm

TYPE (qc), INTENT(IN) :: x
TYPE (qc)             :: z

z%qr = LOG( ABS(x) )
z%qi = ATAN2( x%qi, x%qr )

RETURN
END FUNCTION qc_log


END MODULE quad_prec_complex
