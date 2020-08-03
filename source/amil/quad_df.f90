MODULE quadruple_precision_df
!
! This version is for Compaq Fortran, Lahey/Fujitsu LF95 and
! Absoft ProFortran 6.0.
! N.B. The -o0 option (no optimization) must be used with LF95.
!
! N.B. This is quadruple precision implemented in SOFTWARE, hence it is SLOW.
! This package has NOT been tested with other Fortran 90 compilers.
! The operations +-*/ & sqrt will probably work correctly with other compilers
! for PCs.   Functions and routines which initialize quadruple-precision
! constants, particularly the trigonometric and exponential functions will
! only work if the compiler initializes double-precision constants from
! decimal code in exactly the same way as the Lahey compilers.

! The basic algorithms come from:
! Linnainmma, S.  Software for doubled-precision floating-point computations,
!    ACM Trans, Math. Software, vol. 7, pp.272-283, 1981.

! Programmer : Alan Miller
! e-mail: amiller @ bigpond.net.au   URL:  http://users.bigpond.net.au/amiller

! Latest revision - 18 September 2002

! 4 Aug 97.  Added new algorithm for exponential.  Takes half the time of the
!            previous Taylor series, but errors 2.5 times larger.   The old
!            exponential is included here under the name exp_taylor for anyone
!            who needs the extra accuracy.  To use it instead of the new
!            algorithm, change the module procedure name in INTERFACE exp from
!            longexp to exp_taylor.
! 5 Aug 97.  Found way to reduce cancellation errors in yesterday's algorithm
!            for the exponential.   Removed exp_taylor.
! 18 Aug 97. Added table of quadruple-precision constants.
! 8 Sept 97. Added str_quad which reads a character string and converts it to
!            quadruple-precision form.
! 15 Oct 97. Added quad_str which takes a quadruple-precision value and outputs
!            a character string containing its decimal value to 30 significant
!            digits.
!            Added overlays to the ** operator for quadruple-precision values.
! 15 Jan 98. Added ATAN2.
! 19 Jan 98. Added <, <=, ==, >= and >.
! 27 Dec 98. Added quad/real, quad/integer, integer/quad, real/quad, dp/quad,
!            int+quad, real+quad, dp+quad, int-quad, real-quad, dp-quad,
!            SUM, DOT_PRODUCT & MATMUL.
! 29 Dec 98. Correction to routine string_quad for strings of < 5 characters.
! 10 May 99. Added EPSILON for quad type.
! 5 Oct 99.  log10e corrected.
! 15 Oct 99. Corrected function quad_pow_int.
! 18 Oct 99. Rewrote quad_pow_int to use the binary power method.
! 11 Nov 99. Added overlaying of assignments, e.g. quad = int, etc.
! 21 Jan 00. Added inreface for EPSILON.
! 17 Feb 02. Corrected ACOS & ASIN for arguments close to +1 or -1.
! 18 Feb 02. Further improvement to ACOS.
! 18 Sep 02. Added reference to Linnainmaa in comments at beginning.

IMPLICIT NONE

INTEGER, PARAMETER   :: dp = SELECTED_REAL_KIND(12, 100)
INTEGER, PARAMETER   :: nbits = DIGITS(1.0_dp)
REAL (dp), PARAMETER :: constant = 2._dp**(nbits - nbits/2) + 1._dp
                        !
                        ! Special for Digital Visual Fortran 5.0
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp

PRIVATE :: zero, one, constant, nbits

TYPE    :: quad
  REAL (dp) :: hi, lo
END TYPE quad

INTERFACE OPERATOR (*)
  MODULE PROCEDURE longmul
  MODULE PROCEDURE mult_quad_int
  MODULE PROCEDURE mult_int_quad
  MODULE PROCEDURE mult_quad_dp
  MODULE PROCEDURE mult_dp_quad
  MODULE PROCEDURE mult_quad_real
  MODULE PROCEDURE mult_real_quad
END INTERFACE

INTERFACE OPERATOR (/)
  MODULE PROCEDURE longdiv
  MODULE PROCEDURE div_quad_dp
  MODULE PROCEDURE div_quad_real
  MODULE PROCEDURE div_quad_int
  MODULE PROCEDURE div_int_quad
  MODULE PROCEDURE div_real_quad
  MODULE PROCEDURE div_dp_quad
END INTERFACE

INTERFACE OPERATOR (+)
  MODULE PROCEDURE longadd
  MODULE PROCEDURE quad_add_int
  MODULE PROCEDURE quad_add_real
  MODULE PROCEDURE quad_add_dp
  MODULE PROCEDURE int_add_quad
  MODULE PROCEDURE real_add_quad
  MODULE PROCEDURE dp_add_quad
END INTERFACE

INTERFACE OPERATOR (-)
  MODULE PROCEDURE longsub
  MODULE PROCEDURE quad_sub_int
  MODULE PROCEDURE quad_sub_real
  MODULE PROCEDURE quad_sub_dp
  MODULE PROCEDURE int_sub_quad
  MODULE PROCEDURE real_sub_quad
  MODULE PROCEDURE dp_sub_quad
  MODULE PROCEDURE negate_quad
END INTERFACE

INTERFACE ASSIGNMENT (=)
  MODULE PROCEDURE quad_eq_int
  MODULE PROCEDURE quad_eq_real
  MODULE PROCEDURE quad_eq_dp
  MODULE PROCEDURE int_eq_quad
  MODULE PROCEDURE real_eq_quad
  MODULE PROCEDURE dp_eq_quad
END INTERFACE

INTERFACE OPERATOR (**)
  MODULE PROCEDURE quad_pow_int
  MODULE PROCEDURE quad_pow_real
  MODULE PROCEDURE quad_pow_dp
  MODULE PROCEDURE quad_pow_quad
END INTERFACE

INTERFACE OPERATOR (<)
  MODULE PROCEDURE quad_lt
END INTERFACE

INTERFACE OPERATOR (<=)
  MODULE PROCEDURE quad_le
END INTERFACE

INTERFACE OPERATOR (==)
  MODULE PROCEDURE quad_eq
END INTERFACE

INTERFACE OPERATOR (>=)
  MODULE PROCEDURE quad_ge
END INTERFACE

INTERFACE OPERATOR (>)
  MODULE PROCEDURE quad_gt
END INTERFACE

INTERFACE SCALE
  MODULE PROCEDURE qscale
END INTERFACE

INTERFACE ABS
  MODULE PROCEDURE qabs
END INTERFACE

INTERFACE SQRT
  MODULE PROCEDURE longsqrt
END INTERFACE

INTERFACE LOG
  MODULE PROCEDURE longlog
END INTERFACE

INTERFACE EXP
  MODULE PROCEDURE longexp
END INTERFACE

INTERFACE SIN
  MODULE PROCEDURE longsin
END INTERFACE

INTERFACE COS
  MODULE PROCEDURE longcos
END INTERFACE

INTERFACE TAN
  MODULE PROCEDURE longtan
END INTERFACE

INTERFACE ASIN
  MODULE PROCEDURE longasin
END INTERFACE

INTERFACE ACOS
  MODULE PROCEDURE longacos
END INTERFACE

INTERFACE ATAN
  MODULE PROCEDURE longatan
END INTERFACE

INTERFACE ATAN2
  MODULE PROCEDURE qatan2
END INTERFACE

INTERFACE SUM
  MODULE PROCEDURE quad_sum
END INTERFACE

INTERFACE DOT_PRODUCT
  MODULE PROCEDURE quad_dot_product
END INTERFACE

INTERFACE MATMUL
  MODULE PROCEDURE q_matmul12
  MODULE PROCEDURE q_matmul21
  MODULE PROCEDURE q_matmul22
END INTERFACE

INTERFACE EPSILON
  MODULE PROCEDURE q_epsilon
END INTERFACE

TYPE (quad), PARAMETER ::  &
             pi =  quad( 0.3141592653589793D+01, 0.1224646799147353D-15 ),  &
          piby2 =  quad( 0.1570796326794897D+01, -.3828568698926950D-15 ),  &
          piby3 =  quad( 0.1047197551196598D+01, -.3292527815701405D-15 ),  &
          piby4 =  quad( 0.7853981633974484D+00, -.8040613248383182D-16 ),  &
          piby6 =  quad( 0.5235987755982990D+00, -.1646263907850702D-15 ),  &
          twopi =  quad( 0.6283185307179586D+01, 0.2449293598294707D-15 ),  &
          ln_pi =  quad( 0.1144729885849400D+01, 0.2323105560877391D-15 ),  &
         sqrtpi =  quad( 0.1772453850905516D+01, -.7666586499825800D-16 ),  &
       fact_pt5 =  quad( 0.8862269254527582D+00, -.1493552349616447D-15 ),  &
        sqrt2pi =  quad( 0.2506628274631001D+01, -.6273750096546544D-15 ),  &
      lnsqrt2pi =  quad( 0.9189385332046728D+00, -.3878294158067242D-16 ),  &
      one_on2pi =  quad( 0.1591549430918953D+00, 0.4567181289366658D-16 ),  &
    two_on_rtpi =  quad( 0.1128379167095513D+01, -.4287537502368968D-15 ),  &
        deg2rad =  quad( 0.1745329251994330D-01, -.3174581724866598D-17 ),  &
        rad2deg =  quad( 0.5729577951308232D+02, -.1987849567057628D-14 ),  &
            ln2 =  quad( 0.6931471805599454D+00, -.8783183432405266D-16 ),  &
           ln10 =  quad( 0.2302585092994046D+01, -.2170756223382249D-15 ),  &
          log2e =  quad( 0.1442695040888964D+01, -.6457785410341630D-15 ),  &
         log10e =  quad( 0.4342944819032519D+00, -.5552037773430574D-16 ),  &
        log2_10 =  quad( 0.3321928094887362D+01, 0.1661617516973592D-15 ),  &
        log10_2 =  quad( 0.3010299956639812D+00, -.2803728127785171D-17 ),  &
          euler =  quad( 0.5772156649015330D+00, -.1159652176149463D-15 ),  &
              e =  quad( 0.2718281828459045D+01, 0.1445646891729250D-15 ),  &
          sqrt2 =  quad( 0.1414213562373095D+01, 0.1253716717905022D-15 ),  &
          sqrt3 =  quad( 0.1732050807568877D+01, 0.3223954471431004D-15 ),  &
         sqrt10 =  quad( 0.3162277660168380D+01, -.6348773795572286D-15 )

CONTAINS


FUNCTION exactmul2(a, c) RESULT(ac)
!  Procedure exactmul2, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

REAL (dp), INTENT(IN) :: a, c
TYPE (quad)           :: ac

!     Local variables
REAL (dp) :: a1, a2, c1, c2, t

t = constant * a
a1 = (a - t) + t             ! Lahey's optimization removes the brackets
                             ! and sets a1 = a which defeats the whole point.
a2 = a - a1
t = constant * c
c1 = (c - t) + t
c2 = c - c1
ac%hi = a * c
ac%lo = (((a1 * c1 - ac%hi) + a1 * c2) + c1 * a2) + c2 * a2

RETURN
END FUNCTION exactmul2



FUNCTION longmul(a, c) RESULT(ac)
!  Procedure longmul, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

TYPE (quad), INTENT(IN) :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: zz
TYPE (quad) :: z

z = exactmul2(a%hi, c%hi)
zz = ((a%hi + a%lo) * c%lo + a%lo * c%hi) + z%lo
ac%hi = z%hi + zz
ac%lo = (z%hi - ac%hi) + zz

RETURN
END FUNCTION longmul



FUNCTION mult_quad_int(a, i) RESULT(c)
!     Multiply quadruple-precision number (a) by an integer (i).

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: i
TYPE (quad)             :: c

IF (i == 0) THEN
  c = quad(zero, zero)
ELSE IF (i == 1) THEN
  c = a
ELSE IF (i == -1) THEN
  c = -a
ELSE
  c = exactmul2(a%hi, DBLE(i)) + exactmul2(a%lo, DBLE(i))
END IF

RETURN
END FUNCTION mult_quad_int



FUNCTION mult_int_quad(i, a) RESULT(c)
!     Multiply quadruple-precision number (a) by an integer (i).

INTEGER, INTENT(IN)     :: i
TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: c

IF (i == 0) THEN
  c = quad(zero, zero)
ELSE IF (i == 1) THEN
  c = a
ELSE IF (i == -1) THEN
  c = -a
ELSE
  c = exactmul2(a%hi, DBLE(i)) + exactmul2(a%lo, DBLE(i))
END IF

RETURN
END FUNCTION mult_int_quad



FUNCTION mult_quad_dp(a, b) RESULT(c)
!  Multiply a quadruple-precision number (a) by a double-precision number (b).

TYPE (quad), INTENT(IN) :: a
REAL (dp), INTENT(IN)   :: b
TYPE (quad)             :: c

!     Local variables
TYPE (quad)             :: z
REAL (dp)               :: zz

z = exactmul2(a%hi, b)
zz = a%lo * b + z%lo
c%hi = z%hi + zz
c%lo = (z%hi - c%hi) + zz

RETURN
END FUNCTION mult_quad_dp



FUNCTION mult_quad_real(a, b) RESULT(c)
!  Multiply a quadruple-precision number (a) by a real number (b).

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: b
TYPE (quad)             :: c

!     Local variables
TYPE (quad)             :: z
REAL (dp)               :: zz

z = exactmul2(a%hi, DBLE(b))
zz = a%lo * b + z%lo
c%hi = z%hi + zz
c%lo = (z%hi - c%hi) + zz

RETURN
END FUNCTION mult_quad_real



FUNCTION mult_dp_quad(b, a) RESULT(c)
!  Multiply a quadruple-precision number (a) by a double-precision number (b).

TYPE (quad), INTENT(IN) :: a
REAL (dp), INTENT(IN)   :: b
TYPE (quad)             :: c

!     Local variables
TYPE (quad)             :: z
REAL (dp)               :: zz

z = exactmul2(a%hi, b)
zz = a%lo * b + z%lo
c%hi = z%hi + zz
c%lo = (z%hi - c%hi) + zz

RETURN
END FUNCTION mult_dp_quad



FUNCTION mult_real_quad(a, b) RESULT(c)
!  Multiply a quadruple-precision number (a) by a double-precision number (b).

REAL, INTENT(IN)        :: a
TYPE (quad), INTENT(IN) :: b
TYPE (quad)             :: c

!     Local variables
TYPE (quad)             :: z
REAL (dp)               :: zz

z = exactmul2(b%hi, DBLE(a))
zz = b%lo * a + z%lo
c%hi = z%hi + zz
c%lo = (z%hi - c%hi) + zz

RETURN
END FUNCTION mult_real_quad



FUNCTION longdiv(a, c) RESULT(ac)
!  Procedure longdiv, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

TYPE (quad), INTENT(IN) :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, zz
TYPE (quad) :: q

z = a%hi / c%hi
q = exactmul2(c%hi, z)
zz = ((((a%hi - q%hi) - q%lo) + a%lo) - z*c%lo) / (c%hi + c%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION longdiv



FUNCTION div_quad_dp(a, b) RESULT(c)
!     Divide a quadruple-precision number (a) by a double-precision number (b)

TYPE (quad), INTENT(IN) :: a
REAL (dp), INTENT(IN)   :: b
TYPE (quad)             :: c

!     Local variables
REAL (dp  ) :: z, zz
TYPE (quad) :: q

z = a%hi / b
q = exactmul2(b, z)
zz = (((a%hi - q%hi) - q%lo) + a%lo) / b
c%hi = z + zz
c%lo = (z - c%hi) + zz

RETURN
END FUNCTION div_quad_dp


FUNCTION div_quad_real(a, b) RESULT(c)
!     Divide a quadruple-precision number (a) by a real number (b)

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: b
TYPE (quad)             :: c

!     Local variables
REAL (dp  ) :: z, zz
TYPE (quad) :: q

z = a%hi / b
q = exactmul2(DBLE(b), z)
zz = (((a%hi - q%hi) - q%lo) + a%lo) / b
c%hi = z + zz
c%lo = (z - c%hi) + zz

RETURN
END FUNCTION div_quad_real


FUNCTION div_quad_int(a, b) RESULT(c)
!     Divide a quadruple-precision number (a) by an integer (b)

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: b
TYPE (quad)             :: c

!     Local variables
REAL (dp  ) :: z, zz
TYPE (quad) :: q

z = a%hi / b
q = exactmul2(DBLE(b), z)
zz = (((a%hi - q%hi) - q%lo) + a%lo) / b
c%hi = z + zz
c%lo = (z - c%hi) + zz

RETURN
END FUNCTION div_quad_int



FUNCTION div_int_quad(a, c) RESULT(ac)
!  Procedure longdiv, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

INTEGER, INTENT(IN)     :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, zz
TYPE (quad) :: q

z = DBLE(a) / c%hi
q = exactmul2(c%hi, z)
zz = (((a - q%hi) - q%lo) - z*c%lo) / (c%hi + c%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION div_int_quad



FUNCTION div_real_quad(a, c) RESULT(ac)
!  Procedure longdiv, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

REAL, INTENT(IN)        :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, zz
TYPE (quad) :: q

z = DBLE(a) / c%hi
q = exactmul2(c%hi, z)
zz = (((a - q%hi) - q%lo) - z*c%lo) / (c%hi + c%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION div_real_quad



FUNCTION div_dp_quad(a, c) RESULT(ac)
!  Procedure longdiv, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

REAL (dp), INTENT(IN)   :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, zz
TYPE (quad) :: q

z = a / c%hi
q = exactmul2(c%hi, z)
zz = (((a - q%hi) - q%lo) - z*c%lo) / (c%hi + c%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION div_dp_quad



FUNCTION longadd(a, c) RESULT(ac)
!  Procedure longadd, translated from Pascal, from:
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

TYPE (quad), INTENT(IN) :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi + c%hi
q = a%hi - z
zz = (((q + c%hi) + (a%hi - (q + z))) + a%lo) + c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION longadd



FUNCTION quad_add_int(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_add_int



FUNCTION quad_add_real(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_add_real



FUNCTION quad_add_dp(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
REAL (dp), INTENT(IN)   :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_add_dp



FUNCTION int_add_quad(c, a) RESULT(ac)

INTEGER, INTENT(IN)     :: c
TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION int_add_quad



FUNCTION real_add_quad(c, a) RESULT(ac)

REAL, INTENT(IN)        :: c
TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION real_add_quad



FUNCTION dp_add_quad(c, a) RESULT(ac)

REAL (dp), INTENT(IN)   :: c
TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi + c
q = a%hi - z
zz = (((q + c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION dp_add_quad



FUNCTION longsub(a, c) RESULT(ac)
!  Adapted from longadd by changing signs of c.
!  Linnainmaa, Seppo (1981).   Software for doubled-precision floating-point
!  computations.   ACM Trans. on Math. Software (TOMS), 7, 272-283.

TYPE (quad), INTENT(IN) :: a, c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi - c%hi
q = a%hi - z
zz = (((q - c%hi) + (a%hi - (q + z))) + a%lo) - c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION longsub



FUNCTION quad_sub_int(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi - c
q = a%hi - z
zz = (((q - c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_sub_int



FUNCTION quad_sub_real(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi - c
q = a%hi - z
zz = (((q - c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_sub_real



FUNCTION quad_sub_dp(a, c) RESULT(ac)

TYPE (quad), INTENT(IN) :: a
REAL (dp), INTENT(IN)   :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a%hi - c
q = a%hi - z
zz = (((q - c) + (a%hi - (q + z))) + a%lo)
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION quad_sub_dp



FUNCTION int_sub_quad(a, c) RESULT(ac)

INTEGER, INTENT(IN)     :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = DBLE(a) - c%hi
q = DBLE(a) - z
zz = ((q - c%hi) + (DBLE(a) - (q + z))) - c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION int_sub_quad



FUNCTION real_sub_quad(a, c) RESULT(ac)

REAL, INTENT(IN)        :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a - c%hi
q = a - z
zz = ((q - c%hi) + (a - (q + z))) - c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION real_sub_quad



FUNCTION dp_sub_quad(a, c) RESULT(ac)

REAL (dp), INTENT(IN)   :: a
TYPE (quad), INTENT(IN) :: c
TYPE (quad)             :: ac

!     Local variables
REAL (dp  ) :: z, q, zz

z = a - c%hi
q = a - z
zz = ((q - c%hi) + (a - (q + z))) - c%lo
ac%hi = z + zz
ac%lo = (z - ac%hi) + zz

RETURN
END FUNCTION dp_sub_quad



FUNCTION negate_quad(a) RESULT(b)
!     Change the sign of a quadruple-precision number.
!     In many cases, a & b will occupy the same locations.

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

b%hi = -a%hi
b%lo = -a%lo

RETURN
END FUNCTION negate_quad



SUBROUTINE quad_eq_int(a, i)
!     Assignment

TYPE (quad), INTENT(OUT) :: a
INTEGER, INTENT(IN)      :: i

a%hi = i
a%lo = 0

RETURN
END SUBROUTINE quad_eq_int



SUBROUTINE quad_eq_real(a, r)
!     Assignment

TYPE (quad), INTENT(OUT) :: a
REAL, INTENT(IN)         :: r

a%hi = r
a%lo = zero

RETURN
END SUBROUTINE quad_eq_real



SUBROUTINE quad_eq_dp(a, d)
!     Assignment

TYPE (quad), INTENT(OUT) :: a
REAL (dp), INTENT(IN)    :: d

a%hi = d
a%lo = zero

RETURN
END SUBROUTINE quad_eq_dp



SUBROUTINE int_eq_quad(i, a)
!     Assignment

INTEGER, INTENT(OUT)     :: i
TYPE (quad), INTENT(IN)  :: a

i = a%hi

RETURN
END SUBROUTINE int_eq_quad



SUBROUTINE real_eq_quad(r, a)
!     Assignment

REAL, INTENT(OUT)        :: r
TYPE (quad), INTENT(IN)  :: a

r = a%hi

RETURN
END SUBROUTINE real_eq_quad



SUBROUTINE dp_eq_quad(d, a)
!     Assignment

REAL (dp), INTENT(OUT)   :: d
TYPE (quad), INTENT(IN)  :: a

d = a%hi

RETURN
END SUBROUTINE dp_eq_quad



FUNCTION quad_lt(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = (diff%hi < zero)

RETURN
END FUNCTION quad_lt



FUNCTION quad_le(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = .NOT. (diff%hi > zero)

RETURN
END FUNCTION quad_le



FUNCTION quad_eq(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = (diff%hi == zero)

RETURN
END FUNCTION quad_eq



FUNCTION quad_ge(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = .NOT. (diff%hi < zero)

RETURN
END FUNCTION quad_ge



FUNCTION quad_gt(x, y) RESULT(is_it)
!     Comparison of 2 logical numbers

TYPE (quad), INTENT(IN) :: x, y
LOGICAL                 :: is_it

! Local variable
TYPE (quad) :: diff

diff = x - y
is_it = (diff%hi > zero)

RETURN
END FUNCTION quad_gt



FUNCTION quad_pow_int(a, i) RESULT(b)
!     Raise a quadruple=precision number (a) to a power.

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: i
TYPE (quad)             :: b

! Local variables

INTEGER     :: ia, j, emax
TYPE (quad) :: power
LOGICAL     :: first

SELECT CASE (i)
  CASE (0)
    b = quad(one, zero)
  CASE (-1023:-1, 1:1023)
    ia = ABS(i)
    first = .TRUE.
    power = a
    emax = MAXEXPONENT(one)
    DO j = 0, 9
      IF (BTEST(ia, j)) THEN
        IF (first) THEN
          b = power
          first = .FALSE.
        ELSE
          IF (EXPONENT(power%hi) + EXPONENT(b%hi) > emax) THEN
            WRITE(*, *) '** Exponential overflow in routine QUAD_POW_INT **'
            RETURN
          END IF
          b = b * power
        END IF
        ia = IBCLR(ia, j)
        IF (ia == 0) EXIT
      END IF
      IF (EXPONENT(power%hi) > emax/2) THEN
        WRITE(*, *) '** Exponential overflow in routine QUAD_POW_INT **'
        RETURN
      END IF
      power = power * power
    END DO
    IF (i < 0) b = quad(one, zero) / b
  CASE DEFAULT
    IF (a%hi < zero) THEN
      IF (i * LOG(-a%hi) > LOG( HUGE(one) )) THEN
        WRITE(*, *) '** Exponential overflow in routine QUAD_POW_INT **'
        RETURN
      ELSE IF (i * LOG(-a%hi) < LOG( TINY(one) )) THEN
        b = quad(zero, zero)
      ELSE
        ia = ABS(i)
        IF (BTEST(ia,0)) THEN
          b = - EXP( LOG(-a)*i )         ! i is odd (test least sign. bit)
        ELSE
          b = EXP( LOG(-a)*i )           ! i is even
        END IF
      END IF
    ELSE IF (a%hi > zero) THEN
      IF (i * LOG(a%hi) > LOG( HUGE(one) )) THEN
        WRITE(*, *) '** Exponential overflow in routine QUAD_POW_INT **'
        RETURN
      ELSE IF (i * LOG(a%hi) < LOG( TINY(one) )) THEN
        b = quad(zero, zero)
      ELSE
        b = EXP( LOG(a)*i )
      END IF
    ELSE
      b = quad(zero, zero)
    END IF
END SELECT

RETURN
END FUNCTION quad_pow_int



FUNCTION quad_pow_real(a, r) RESULT(b)
!     Raise a quadruple=precision number (a) to a power.

TYPE (quad), INTENT(IN) :: a
REAL, INTENT(IN)        :: r
TYPE (quad)             :: b

IF (a%hi < zero) THEN
  WRITE(*, *)   &
  ' *** Error: attempt to raise negative quad. prec. number to a real power ***'
  b = quad(zero, zero)
ELSE IF (a%hi > zero) THEN
  b = EXP( LOG(a)*r )
ELSE
  b = quad(zero, zero)
END IF

RETURN
END FUNCTION quad_pow_real



FUNCTION quad_pow_dp(a, d) RESULT(b)
!     Raise a quadruple=precision number (a) to a power.

TYPE (quad), INTENT(IN) :: a
REAL (dp), INTENT(IN)   :: d
TYPE (quad)             :: b

IF (a%hi < zero) THEN
  WRITE(*, *)   &
  ' *** Error: attempt to raise negative quad. prec. number to a real power ***'
  b = quad(zero, zero)
ELSE IF (a%hi > zero) THEN
  b = EXP( LOG(a)*d )
ELSE
  b = quad(zero, zero)
END IF

RETURN
END FUNCTION quad_pow_dp



FUNCTION quad_pow_quad(a, q) RESULT(b)
!     Raise a quadruple=precision number (a) to a power.

TYPE (quad), INTENT(IN) :: a, q
TYPE (quad)             :: b

IF (a%hi < zero) THEN
  WRITE(*, *)   &
  ' *** Error: attempt to raise negative quad. prec. number to a real power ***'
  b = quad(zero, zero)
ELSE IF (a%hi > zero) THEN
  b = EXP( LOG(a)*q )
ELSE
  b = quad(zero, zero)
END IF

RETURN
END FUNCTION quad_pow_quad



FUNCTION qscale(a, i) RESULT(b)
!     Multiply a by 2^i

TYPE (quad), INTENT(IN) :: a
INTEGER, INTENT(IN)     :: i
TYPE (quad)             :: b

b%hi = SCALE(a%hi, i)
b%lo = SCALE(a%lo, i)

RETURN
END FUNCTION qscale



FUNCTION qabs(a) RESULT(b)
!     Absolute value of a quadruple-precision number

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

IF (a%hi < zero) THEN
  b%hi = - a%hi
  b%lo = - a%lo
ELSE
  b%hi = a%hi
  b%lo = a%lo
END IF

RETURN
END FUNCTION qabs



FUNCTION longsqrt(a) RESULT(b)
! This is modified from procedure sqrt2 of:
!    Dekker, T.J. (1971). 'A floating-point technique for extending the
!    available precision', Numer. Math., 18, 224-242.

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

!     Local variables
REAL (dp)   :: t, res
TYPE (quad) :: tt

! Check that ahi >= 0.

IF (a%hi < 0._dp) THEN
  WRITE(*, *) ' *** Negative argument for longsqrt ***'
  RETURN
ELSE IF (a%hi == 0.d0) THEN
  b%hi = 0._dp
  b%lo = 0._dp
  RETURN
END IF

! First approximation is  t = sqrt(a).

t = SQRT(a%hi)
tt = exactmul2(t, t)
res = t + (((a%hi - tt%hi) - tt%lo) + a%lo) * 0.5_dp / t
b%lo = (t - res) + (((a%hi - tt%hi) - tt%lo) + a%lo) * 0.5_dp / t
b%hi = res

RETURN
END FUNCTION longsqrt



FUNCTION longlog(x) RESULT(y)
!  Quadruple-precision logarithm to base e
!  Halley's algorithm using double-precision logarithm as starting value.
!  Solves:         y
!          f(y) = e  - x = 0

TYPE (quad), INTENT(IN) :: x
TYPE (quad)             :: y

!     Local variables
TYPE (quad) :: expy, f

y%hi = LOG(x%hi)
y%lo = 0._dp
expy = EXP(y)
f = expy - x
f = SCALE(f, 1)
y = y - f / (expy + x)

RETURN
END FUNCTION longlog



FUNCTION longexp(x) RESULT(y)
!  Calculate a quadruple-precision exponential
!  Method:
!   x    x.log2(e)    nint[x.log2(e)] + frac[x.log2(e)]
!  e  = 2          = 2
!
!                     iy    fy
!                  = 2   . 2
!  Then
!   fy    y.ln(2)
!  2   = e
!
!  Now y.ln(2) will be less than 0.3466 in absolute value.
!  This is halved and a Pade approximation is used to approximate e^x over
!  the region (-0.1733, +0.1733).   This approximation is then squared.

!  WARNING: No overflow checks!

TYPE (quad), INTENT(IN) :: x
TYPE (quad)             :: y

! Local variables
TYPE (quad)          :: temp, ysq, sum1, sum2
INTEGER              :: iy

y = x / ln2
iy = NINT(REAL(y%hi))
y = (y - DBLE(iy)) * ln2
y = SCALE(y, -1)

! The Pade series is:
!     p0 + p1.y + p2.y^2 + p3.y^3 + ... + p9.y^9
!     ------------------------------------------
!     p0 - p1.y + p2.y^2 - p3.y^3 + ... - p9.y^9
!
! sum1 is the sum of the odd powers, sum2 is the sum of the even powers

ysq = y * y
sum1 = y * ((((ysq + 3960.)*ysq + 2162160._dp)*ysq + 302702400._dp)*ysq +   &
               8821612800._dp)
sum2 = (((90.*ysq + 110880.)*ysq + 30270240._dp)*ysq + 2075673600._dp)*ysq +  &
       17643225600._dp

!                     sum2 + sum1         2.sum1
! Now approximation = ----------- = 1 + ----------- = 1 + 2.temp
!                     sum2 - sum1       sum2 - sum1
!
! Then (1 + 2.temp)^2 = 4.temp.(1 + temp) + 1

temp = sum1 / (sum2 - sum1)
y = temp * (temp + 1._dp)
y = SCALE(y, 2)
y = y + 1._dp
y = SCALE(y, iy)

RETURN
END FUNCTION longexp



SUBROUTINE longmodr(a, b, n, rem)

! Extended arithmetic calculation of the 'rounded' modulus:
!  a = n.b + rem
! where all quantities are in quadruple-precision, except the integer
! number of multiples, n.   The absolute value of the remainder (rem)
! is not greater than b/2.
! The result is exact.   rem may occupy the same location as either input.

! Programmer: Alan Miller

! Latest revision - 11 September 1986
! Fortran version - 4 December 1996

TYPE (quad), INTENT(IN)  :: a, b
INTEGER, INTENT(OUT)     :: n
TYPE (quad), INTENT(OUT) :: rem

! Local variables

TYPE (quad) :: temp

! Check that b%hi .ne. 0

IF (b%hi == 0._dp) THEN
  WRITE(*, *) ' *** Error in longmodr - 3rd argument zero ***'
  RETURN
END IF

! Calculate n.

temp = a / b
n = NINT(REAL(temp%hi))

! Calculate remainder preserving full accuracy.

temp = exactmul2(DBLE(n), b%hi)
rem%hi = a%hi
rem%lo = zero
temp = rem - temp
rem%hi = a%lo
rem%lo = zero
temp = rem + temp
rem = exactmul2(DBLE(n), b%lo)
rem = temp - rem

RETURN
END SUBROUTINE longmodr



! Extended accuracy arithmetic sine, cosine & tangent (about 31 decimals).
! Calculates  b = sin, cos or tan (a), where all quantities are in
! quadruple-precision, using table look-up and a Taylor series expansion.
! The result may occupy the same locations as the input value.
! Much of the code is common to all three functions, and this is in a
! subroutine longcst.


FUNCTION longsin(a) RESULT(b)

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables

LOGICAL :: sine, cosine, tangent

! Set logical variables for sine function.

sine = .true.
cosine = .false.
tangent = .false.
CALL longcst(a, b, sine, cosine, tangent)

RETURN
END FUNCTION longsin



FUNCTION longcos(a) RESULT(b)

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables

LOGICAL :: sine, cosine, tangent

! Set logical variables for sine function.

sine = .false.
cosine = .true.
tangent = .false.
CALL longcst(a, b, sine, cosine, tangent)

RETURN
END FUNCTION longcos



FUNCTION longtan(a) RESULT(b)

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables

LOGICAL :: sine, cosine, tangent

! Set logical variables for sine function.

sine = .false.
cosine = .false.
tangent = .true.
CALL longcst(a, b, sine, cosine, tangent)

RETURN
END FUNCTION longtan



SUBROUTINE longcst(a, b, sine, cosine, tangent)

TYPE (quad), INTENT(IN)  :: a
TYPE (quad), INTENT(OUT) :: b
LOGICAL, INTENT(IN)      :: sine, cosine, tangent

! Local variables

LOGICAL     :: pos
TYPE (quad) :: d, term, temp, angle, piby40, sum1, sum2, sin
INTEGER     :: npi, ipt, i
REAL (dp)   :: tol15 = 1.E-15_dp, tol30 = 1.E-30_dp

! sin(i.pi/40), i = 0(1)20
REAL (dp)   :: table(2, 0:20) = RESHAPE( (/    &
               0.0000000000000000E+00_dp,  0.0000000000000000E+00_dp,  &
               0.7845909572784494E-01_dp,  0.1464397249532491E-17_dp,  &
               0.1564344650402309E+00_dp,  -.2770509565052586E-16_dp,  &
               0.2334453638559054E+00_dp,  0.2058612230858154E-16_dp,  &
               0.3090169943749475E+00_dp,  -.8267172724967036E-16_dp,  &
               0.3826834323650898E+00_dp,  -.1005077269646159E-16_dp,  &
               0.4539904997395468E+00_dp,  -.1292033036231312E-16_dp,  &
               0.5224985647159488E+00_dp,  0.6606794454708078E-16_dp,  &
               0.5877852522924732E+00_dp,  -.1189570533007057E-15_dp,  &
               0.6494480483301838E+00_dp,  -.1134903961116171E-15_dp,  &
               0.7071067811865476E+00_dp,  -.4833646656726458E-16_dp,  &
               0.7604059656000310E+00_dp,  -.1036987135483477E-15_dp,  &
               0.8090169943749476E+00_dp,  -.1381828784809282E-15_dp,  &
               0.8526401643540922E+00_dp,  0.4331886637554353E-16_dp,  &
               0.8910065241883680E+00_dp,  -.1474714419679880E-15_dp,  &
               0.9238795325112868E+00_dp,  -.9337725537817898E-16_dp,  &
               0.9510565162951536E+00_dp,  -.7008780156242836E-16_dp,  &
               0.9723699203976766E+00_dp,  0.4478912629332321E-16_dp,  &
               0.9876883405951378E+00_dp,  -.4416018005989794E-16_dp,  &
               0.9969173337331280E+00_dp,  0.1235153006196267E-16_dp,  &
               0.1000000000000000E+01_dp,  0.0000000000000000E+00_dp  /),  &
               (/ 2, 21 /) )

! pi/40

piby40%hi = 0.7853981633974484E-01_dp
piby40%lo = -.1081617080994607E-16_dp

! Reduce angle to range (-pi/2, +pi/2) by subtracting an integer multiple of pi.

CALL longmodr(a, pi, npi, angle)

! Find nearest multiple of pi/40 to angle.

CALL longmodr(angle, piby40, ipt, d)

! Sum 1 = 1 - d**2/2! + d**4/4! - d**6/6! + ...
! Sum 2 = d - d**3/3! + d**5/5! - d**7/7! + ...

sum1%hi = zero
sum1%lo = zero
sum2%hi = zero
sum2%lo = zero
pos = .false.
term = d
i = 2
20 IF (ABS(term%hi) > tol15) THEN
  term = term * d                                ! Use quad. precision
  IF (i == 2 .OR. i == 4 .OR. i == 8) THEN
    term%hi = term%hi / i
    term%lo = term%lo / i
  ELSE
    term = term / DBLE(i)
  END IF
  IF (pos) THEN
    sum1 = sum1 + term
  ELSE
    sum1 = sum1 - term
  END IF
ELSE
  term%hi = term%hi * d%hi / DBLE(i)             ! Double prec. adequate
  IF (pos) THEN
    sum1%lo = sum1%lo + term%hi
  ELSE
    sum1%lo = sum1%lo - term%hi
  END IF
END IF

! Repeat for sum2

i = i + 1
IF (ABS(term%hi) > tol15) THEN
  term = term * d / DBLE(i)                      ! Use quad. precision
  IF (pos) THEN
    sum2 = sum2 + term
  ELSE
    sum2 = sum2 - term
  END IF
ELSE
  term%hi = term%hi * d%hi / DBLE(i)             ! Double prec. adequate
  IF (pos) THEN
    sum2%lo = sum2%lo + term%hi
  ELSE
    sum2%lo = sum2%lo - term%hi
  END IF
END IF

i = i + 1
pos = .NOT. pos
IF (ABS(term%hi) > tol30) GO TO 20

sum1 = sum1 + 1._dp                              ! Now add the 1st terms
sum2 = sum2 + d                                  ! for max. accuracy

! Construct sine, cosine or tangent.
! Sine first.    sin(angle + d) = sin(angle).cos(d) + cos(angle).sin(d)

IF (sine .OR. tangent) THEN
  IF (ipt >= 0) THEN
    temp%hi = table(1, ipt)
    temp%lo = table(2, ipt)
  ELSE
    temp%hi = - table(1, -ipt)
    temp%lo = - table(2, -ipt)
  END IF
  b = sum1 * temp
  IF (ipt >= 0) THEN
    temp%hi = table(1, 20-ipt)
    temp%lo = table(2, 20-ipt)
  ELSE
    temp%hi = table(1, 20+ipt)
    temp%lo = table(2, 20+ipt)
  END IF
  b = b + sum2 * temp
  IF (npi /= 2*(npi/2)) THEN
    b = -b
  END IF
  IF (tangent) THEN
    sin = b
  END IF
END IF

! Cosine or tangent.

IF (cosine .OR. tangent) THEN
  IF (ipt >= 0) THEN
    temp%hi = table(1, ipt)
    temp%lo = table(2, ipt)
  ELSE
    temp%hi = - table(1, -ipt)
    temp%lo = - table(2, -ipt)
  END IF
  b = sum2 * temp
  IF (ipt >= 0) THEN
    temp%hi = table(1, 20-ipt)
    temp%lo = table(2, 20-ipt)
  ELSE
    temp%hi = table(1, 20+ipt)
    temp%lo = table(2, 20+ipt)
  END IF
  b = sum1 * temp - b
  IF (npi /= 2*(npi/2)) THEN
    b = -b
  END IF
END IF

! Tangent.

IF (tangent) THEN

! Check that bhi .ne. 0

  IF (b%hi == 0.d0) THEN
    WRITE(*, *) ' *** Infinite tangent - routine longcst ***'
    b%hi = HUGE(one)
    b%lo = 0.d0
    RETURN
  END IF
  b = sin / b
END IF

RETURN
END SUBROUTINE longcst



FUNCTION longasin(a) RESULT(b)

! Quadratic-precision arc sine (about 31 decimals).
! One Newton-Raphson iteration to solve:  f(b) = sin(b) - a = 0,
! except when a close to -1 or +1.
! The result (b) may occupy the same location as the input values (a).
! Use ACOS when |a| is close to 1.

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables
TYPE (quad)  :: y, c

! Check that -1 <= a%hi <= +1.

IF (a%hi < -1._dp .OR. a%hi > 1._dp) THEN
  WRITE(*, *) ' *** Argument outside range for longasin ***'
  RETURN
END IF

IF (ABS(a%hi) < 0.866) THEN
  ! First approximation is  y = asin(a).
  ! Quadruple-precision result is  y - [sin(y) - a]/cos(y).

  y%hi = ASIN(a%hi)
  y%lo = zero
  b = y + (a - SIN(y)) / COS(y%hi)
ELSE
  ! Calculate acos(c) where c = sqrt(1 - a^2)
  c = SQRT(one - a*a)
  y%hi = ACOS(c%hi)
  y%lo = zero
  b = y + (COS(y) - c) / SIN(y%hi)
  IF (a%hi < zero) b = -b
END IF


RETURN
END FUNCTION longasin



FUNCTION longacos(a) RESULT(b)

! Quadratic-precision arc cosine (about 31 decimals).
! Newton-Raphson iteration to solve: f(b) = cos(b) - a = 0.
! The result (b) may occupy the same location as the input values (a).
! When |a| is near 1, use formula from p.175 of
! `Software Manual for the Elementary Functions' by W.J. Cody, Jr. &
! W. Waite, Prentice-Hall, 1980.

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables
TYPE (quad)  :: y, c

! Check that -1 <= a%hi <= +1.

IF (a%hi < -1._dp .OR. a%hi > 1._dp) THEN
  WRITE(*, *) ' *** Argument outside range for longacos ***'
  RETURN
END IF

IF (ABS(a%hi) < 0.866) THEN
  ! First approximation is  y = acos(a).
  ! Quadruple-precision result is  y + [cos(y) - a]/sin(y).

  y%hi = ACOS(a%hi)
  y%lo = 0._dp
  b = y + (COS(y) - a) / SIN(y%hi)
ELSE
  ! Calculate 2.asin(c) where c = sqrt([1 - |a|]/2)
  c = SQRT((one - ABS(a))/2)
  y%hi = ASIN(c%hi)
  y%lo = zero
  b = (y - (SIN(y) - c) / COS(y%hi))*2
  IF (a%hi < zero) b = pi - b
END IF

RETURN
END FUNCTION longacos



FUNCTION longatan(a) RESULT(b)

! Quadratic-precision arc tangent (about 31 decimals).
! Newton-Raphson iteration to solve: f(b) = tan(b) - a = 0.
! The result (b) may occupy the same location as the input values (a).

TYPE (quad), INTENT(IN) :: a
TYPE (quad)             :: b

! Local variables
TYPE (quad)  :: y

! First approximation is  y = atan(a).
! Quadruple-precision result is  y - [tan(y) - a] * cos(y)**2.

y%hi = ATAN(a%hi)
y%lo = 0._dp
b = y - (TAN(y) - a) * (COS(y%hi))**2

RETURN
END FUNCTION longatan



FUNCTION qatan2(y, x) RESULT(b)

! Quadratic-precision arc tangent (about 31 decimals).
! As for arc tangent (y/x) except that the result is in the range
!       -pi < ATAN2 <= pi.
! The signs of x and y determine the quadrant.

TYPE (quad), INTENT(IN) :: y, x
TYPE (quad)             :: b

! Local variables
TYPE (quad)  :: z

! First approximation is  z = atan2(y, x).
! Quadruple-precision result is  z - [tan(z) - (y/x)] * cos(z)**2.

z%hi = ATAN2(y%hi, x%hi)
z%lo = 0._dp
IF (x%hi == zero) THEN
  b = z
ELSE
  b = z - (TAN(z) - y/x) * (COS(z%hi))**2
END IF

RETURN
END FUNCTION qatan2



FUNCTION quad_sum(a) RESULT(s)

! Quadruple-precision SUM

TYPE (quad), INTENT(IN), DIMENSION(:) :: a
TYPE (quad)                           :: s

! Local variables

INTEGER :: i

s = quad(zero, zero)

DO i = 1, SIZE(a)
  s = s + a(i)
END DO

RETURN
END FUNCTION quad_sum



FUNCTION quad_dot_product(a, b) RESULT(ab)

! Quadruple-precision DOT_PRODUCT

TYPE (quad), INTENT(IN), DIMENSION(:) :: a, b
TYPE (quad)                           :: ab

! Local variables

INTEGER :: i, n

ab = quad(zero, zero)
n = SIZE(a)
IF (n /= SIZE(b)) THEN
  WRITE(*, *) ' ** Error invoking DOT_PRODUCT - different argument sizes **'
  WRITE(*, '(a, i10, a, i10)') ' Size of 1st argument = ', n,   &
                               '   Size of 2nd argument = ', SIZE(b)
  RETURN
END IF

DO i = 1, n
  ab = ab + a(i)*b(i)
END DO

RETURN
END FUNCTION quad_dot_product



FUNCTION q_matmul12(a, b) RESULT(ab)

! Quadruple-precision MATMUL
! Rank of A = 1, rank of B = 2

TYPE (quad), INTENT(IN), DIMENSION(:)   :: a
TYPE (quad), INTENT(IN), DIMENSION(:,:) :: b
TYPE (quad), DIMENSION( SIZE(b,2) )     :: ab

! Local variables

INTEGER :: j, na, nb1, nb2

! Check dimensions

na = SIZE(a)
nb1 = SIZE(b, 1)
nb2 = SIZE(b, 2)
IF (na /= nb1) THEN
  WRITE(*, *) ' ** Incompatible dimensions for quad-prec. MATMUL'
  RETURN
END IF

DO j = 1, nb2
  ab(j) = quad_dot_product( a, b(:,j) )
END DO

RETURN
END FUNCTION q_matmul12



FUNCTION q_matmul21(a, b) RESULT(ab)

! Quadruple-precision MATMUL
! Rank of A = 2, rank of B = 1

TYPE (quad), INTENT(IN), DIMENSION(:,:) :: a
TYPE (quad), INTENT(IN), DIMENSION(:)   :: b
TYPE (quad), DIMENSION( SIZE(a,1) )     :: ab

! Local variables

INTEGER :: i, na1, na2, nb

! Check dimensions

na1 = SIZE(a, 1)
na2 = SIZE(a, 2)
nb = SIZE(b)

IF (na2 /= nb) THEN
  WRITE(*, *) ' ** Incompatible dimensions for quad-prec. MATMUL'
  RETURN
END IF

DO i = 1, na1
  ab(i) = quad_dot_product( a(i,:), b )
END DO

RETURN
END FUNCTION q_matmul21



FUNCTION q_matmul22(a, b) RESULT(ab)

! Quadruple-precision MATMUL
! Rank of A = 2, rank of B = 2

TYPE (quad), INTENT(IN), DIMENSION(:,:)     :: a
TYPE (quad), INTENT(IN), DIMENSION(:,:)     :: b
TYPE (quad), DIMENSION(SIZE(a,1),SIZE(b,2)) :: ab

! Local variables

INTEGER :: i, j, na1, na2, nb1, nb2

! Check dimensions

na1 = SIZE(a, 1)
na2 = SIZE(a, 2)
nb1 = SIZE(b, 1)
nb2 = SIZE(b, 2)

IF (na2 /= nb1) THEN
  WRITE(*, *) ' ** Incompatible dimensions for quad-prec. MATMUL'
  RETURN
END IF

DO i = 1, na1
  DO j = 1, nb2
    ab(i,j) = quad_dot_product( a(i,:), b(:,j) )
  END DO
END DO

RETURN
END FUNCTION q_matmul22



SUBROUTINE string_quad(string, value, ier)
! Convert a character string to a quadruple-precision quantity.
! Error indicator ier = 0 if value OK
!                     = 1 if string has > 50 characters and no decimal point
!                            in about the first 45
!                     = 2 if a non-numeric character is read in the mantissa
!                            part

CHARACTER (LEN=*), INTENT(IN) :: string
TYPE (quad), INTENT(OUT)      :: value
INTEGER, INTENT(OUT)          :: ier

! Local variables
CHARACTER (LEN=50)            :: str
INTEGER                       :: length, pos, decpt, status, power10, i, i1, i2
REAL (dp)                     :: sgn, temp
CHARACTER (LEN= 5)            :: expnt

str = ADJUSTL(string)
length = LEN_TRIM( ADJUSTL(string) )
IF (length > 50) THEN
                             ! Truncate string to 50 characters preserving
                             ! the exponent, if present.
  pos = SCAN(string(length-4:length), 'DdEe')
  IF (pos > 0) THEN
    str(45+pos:50) = string(length-5+pos:length)
  END IF
                             ! Check that the string contains a '.'
  IF (INDEX(str, '.') == 0) THEN
    ier = 1
    value = quad(zero, zero)
    RETURN
  END IF
END IF

ier = 0

! Determine sign of mantissa

sgn = +1._dp
IF (str(1:1) == '-') THEN
  sgn = -1._dp
  str = str(2:)
ELSE IF (str(1:1) == '+') THEN
  str = str(2:)
END IF

! Separate the exponent, if there is one

length = LEN_TRIM(str)
IF (length > 4) THEN
  pos = SCAN(str(length-4:length), 'DdEe')
  IF (pos > 0) THEN
    expnt = str(length-5+pos:)
    str(length-5+pos:) = '     '
  ELSE
    expnt = ' '
  END IF
ELSE
  expnt = ' '
END IF

! decpt = position of decimal point

decpt = INDEX(str, '.')
IF (decpt > 0) str = str(:decpt-1) // str(decpt+1:)

! Read str in blocks of up to 9 digits at a time, as a large integer.

i1 = 1
length = LEN_TRIM(str)
value = quad(zero, zero)
DO
  i2 = MIN(i1+8, length)
  READ(str(i1:i2), '(f9.0)', IOSTAT=status) temp
  IF (status /= 0) THEN
    ier = 2
    RETURN
  END IF
  IF (i1 > 1) THEN
    value = value * (10._dp ** (i2+1-i1)) + temp
  ELSE
    value = quad(temp, zero)
  END IF
  i1 = i2 + 1
  IF (i1 > length) EXIT
END DO

IF (sgn < zero) value = -value

! Multiply by appropriate power of 10 for position of the decimal point,
! and the exponent.

IF (expnt == ' ') THEN
  power10 = 0
ELSE
  i = LEN_TRIM(expnt)
  READ(expnt(2:i), '(i4)') power10
END IF
IF (decpt > 0) power10 = power10 - (length + 1 - decpt)

! As 1.D+15 is represented exactly in double-precision IEEE arithmetic,
! multiply by multiples of 1.D+15 or 1.D-15.

IF (power10 == 0) RETURN
IF (power10 < 0) THEN
  DO i = 1, -power10/15
    value = value / 1.D+15
  END DO
  i = MOD(-power10, 15)
  IF (i /= 0) value = value / (10._dp ** i)
ELSE
  DO i = 1, power10/15
    value = value * 1.D+15
  END DO
  i = MOD(power10, 15)
  IF (i /= 0) value = value * (10._dp ** i)
END IF

RETURN
END SUBROUTINE string_quad


SUBROUTINE quad_string(value, string, ier)
! Convert a quadruple-precision quantity to a decimal character string.
! Error indicator ier = 0 if conversion OK
!                     = 1 if the length of the string < 36 characters.

TYPE (quad), INTENT(IN)        :: value
CHARACTER (LEN=*), INTENT(OUT) :: string
INTEGER, INTENT(OUT)           :: ier

! Local variables
CHARACTER (LEN= 1)   :: sgn
CHARACTER (LEN=17)   :: str1, str2
TYPE (quad)          :: val
INTEGER              :: dec_expnt, i
REAL (dp)            :: tmp

IF (LEN(string) < 36) THEN
  ier = 1
  string = '***'
  RETURN
END IF
ier = 0

! Check if value = zero.
IF (value%hi == zero) THEN
  string = ' 0.00'
  RETURN
END IF

IF (value%hi < zero) THEN
  sgn = '-'
  val = - value
ELSE
  sgn = ' '
  val = value
END IF

! Use LOG10 to set the exponent.
dec_expnt = FLOOR( LOG10(val%hi) )

! Get the first 15 decimal digits
IF (dec_expnt /= 14) THEN
  val = val * EXP( LOG(quad(10._dp, zero)) * REAL(14 - dec_expnt) )
END IF
WRITE(str1, '(f16.0)') val%hi

! Calculate the remainder
READ(str1, '(f16.0)') tmp
val = val - tmp

! If val is -ve, subtract 1 from the last digit of str1, and add 1 to val.
IF (val%hi < -0.5D-16) THEN
  tmp = tmp - one
  WRITE(str1, '(f16.0)') tmp
  val = val + one
END IF
val = val * 1.D15

! write the second 15 digits
WRITE(str2, '(f16.0)') val%hi

! If str2 consists of asterisks, add 1 in the last digit to str1.
! Set str2 to zeroes.
IF (str2(2:2) == '*') THEN
  tmp = tmp + one
  WRITE(str1, '(f17.0)') tmp
  IF (str1(1:1) /= ' ') THEN
    dec_expnt = dec_expnt + 1
  ELSE
    str1 = str1(2:17)
  END IF
  str2 = '000000000000000.'
END IF

! Replace leading blanks with zeroes
DO i = 1, 15
  IF (str2(i:i) /= ' ') EXIT
  str2(i:i) = '0'
END DO

! Combine str1 & str2, removing decimal points & adding exponent.
i = INDEX(str1, '.')
str1(i:i) = ' '
str2(16:16) = ' '
string = '.' // TRIM(ADJUSTL(str1)) // TRIM(ADJUSTL(str2)) // 'E'
WRITE(str1, '(i4.2)') dec_expnt+1
string = TRIM(string) // ADJUSTL(str1)

! Restore the sign.
IF (sgn == '-') THEN
  string = '-' // ADJUSTL(string)
ELSE
  string = ADJUSTL(string)
END IF

RETURN
END SUBROUTINE quad_string


FUNCTION q_epsilon(x) RESULT(eps)
! Returns the machine accuracy.
! eps is the smallest value such that x.(1 + eps) > x for all x.
! This value is machine dependent.   It has been checked only for the
! the following PC compilers:
! Lahey's ELF90
! Lahey/Fujitsi LF95
! Compaq (Digital) DVF5

TYPE (quad), INTENT(IN)  :: x
TYPE (quad)              :: eps

eps = quad(0.6163e-32_dp, 0.0_dp)

RETURN
END FUNCTION q_epsilon

END MODULE quadruple_precision_df
