! This file now contains TOMS algorithm 668 as a separate module.
! Latest revision - 24 July 2000
! amiller@bigpond.net.au
! http://users.bigpond.net.au/amiller/


MODULE random_hypergeometric
! Generate random deviates from the hypergeometric distribution.
! The notation used for the 2 x 2 table is:

!   x  r-x |  r
!          |
!  g-x     | N-r
! _________|_____
!          |
!   g  N-g |  N

! where N, r, g are fixed, and x is the random value to be generated.

! Two methods are included.

! r_hyperg is a simple method which is suitable when MIN(r,g) is small,
!          say <= 10.

! h_alias is much faster for larger values of MIN(r,g) but requires a table
!         of hypergeometric probabilities to be set up on the first call.
!         It is not suitable if the marginal totals (N, r, g) change.

! Latest revision - 24 July 2000

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), ALLOCATABLE, SAVE  :: p(:), f(:)
INTEGER, ALLOCATABLE, SAVE    :: l(:)
INTEGER, SAVE                 :: Nlast = 0, rlast = 0, glast = 0, xmax


CONTAINS


FUNCTION r_hyperg(N, r, g) RESULT(x)

! This is algorithm 3.1 from page 52 of:
!   Dagpunar, John (1988) `Principles of Random Variate Generation'
!   published by Clarendon Press, Oxford.   ISBN 0-19-852202-9

INTEGER, INTENT(IN)  :: N, r, g
INTEGER              :: x

! Local variables

INTEGER  :: gg, i, NN
REAL     :: p, u

gg = g
NN = N
i = 1
x = 0
DO
  IF (gg <= 0 .OR. i > r) EXIT
  p = REAL(gg) / NN
  CALL RANDOM_NUMBER(u)
  IF (u < p) THEN
    x = x + 1
    gg = gg - 1
  END IF
  i = i + 1
  NN = NN - 1
END DO

RETURN
END FUNCTION r_hyperg



FUNCTION h_alias(N, r, g) RESULT(x)

! This is algorithm HALIAS from:

! Kachitvichyanukul, V. & Schmeiser, B. (1985)
! `Computer Generation of Hypergeometric Random Variates',
! J. of Statistical Computation & Simulation, vol.22, 127-145.

INTEGER, INTENT(IN)  :: N, r, g
INTEGER              :: x

! Local variables
INTEGER    :: i, i0, nn(2,2)
REAL (dp)  :: u, p_last

! If this is the first call, or N, r or g has changed, call alias
! to set up new tables.

IF (N /= Nlast .OR. r /= rlast .OR. g /= glast) THEN
  Nlast = N
  rlast = r
  glast = g
  xmax = MIN(r,g)
  IF (ALLOCATED( p )) DEALLOCATE( p, f, l )
  ALLOCATE( p(0:xmax), f(0:xmax), l(0:xmax) )

! Calculate the set of hypergeometric probabilities
! First calculate the maximum probability which occurs near
! i = (r+1)(g+1)/(N+2)

  i0 = NINT( REAL((r+1)*(g+1)) / REAL(N+2) )
  i0 = MIN(i0, xmax)
  nn(1,1) = i0
  nn(1,2) = r - i0
  nn(2,1) = g - i0
  nn(2,2) = N - r - g + i0
  p(i0) = hyperg(nn)

  p_last = p(i0)
  DO i = i0-1, 0, -1
    p(i) = p_last * REAL( (i+1)*(N-r-g+i+1), KIND=dp ) /  &
                    REAL( (r-i)*(g-i), KIND=dp )
    p_last = p(i)
  END DO

  p_last = p(i0)
  DO i = i0+1, xmax
    p(i) = p_last * REAL( (r-i+1)*(g-i+1), KIND=dp ) /  &
                    REAL( i*(N-r-g+i), KIND=dp )
    p_last = p(i)
  END DO

! Call ALIAS to set up its tables of F and L.

  CALL alias(p, xmax+1, f, l)
END IF

! Generate an integer in the range 0 to XMAX
CALL RANDOM_NUMBER(u)
x = u * (xmax + 1)

! Accept it if it is <= f(x), otherwise use its alias.
CALL RANDOM_NUMBER(u)
IF (u > f(x)) x = l(x) - 1

RETURN
END FUNCTION h_alias



FUNCTION hyperg(n) RESULT(prob)

!  Calculate hypergeometric probability for a 2 x 2 table:
!      n(1,1)  n(1,2)
!      n(2,1)  n(2,2)

!  For large numbers, it would be faster to use Stirling's
!  approximation or similar for the factorials.

!  Programmer: Alan Miller
!  Latest revision - 9 January 1997

IMPLICIT NONE
INTEGER, INTENT(IN)  :: n(2,2)
REAL (dp)            :: prob

!     Local variables

INTEGER    :: ncell(1:4), ntot(1:4), nall
LOGICAL    :: ended
REAL (dp)  :: zero = 0.0_dp, one = 1.0_dp

!     Check that all elements are positive or zero.

IF (n(1,1) < 0 .OR. n(1,2) < 0 .OR. n(2,1) < 0 .OR. n(2,2) < 0) THEN
  WRITE(*, *) '** Negative cell frequency in argument to HYPERG **'
  prob = zero
  RETURN
END IF

!     Copy the cell frequencies into ncell and order them.
!     Calculate the row & column totals and order them.

ncell(1) = n(1,1)
ncell(2) = n(1,2)
ncell(3) = n(2,1)
ncell(4) = n(2,2)
CALL order4(ncell)

ntot(1) = n(1,1) + n(1,2)
ntot(2) = n(2,1) + n(2,2)
ntot(3) = n(1,1) + n(2,1)
ntot(4) = n(1,2) + n(2,2)
CALL order4(ntot)

nall = ncell(1) + ncell(2) + ncell(3) + ncell(4)
IF (nall == 0) THEN
  WRITE(*, *) '** Call to HYPERG with total cell frequency = 0 **'
  prob = zero
  RETURN
END IF

!     Calculate hyperg as:

!       ntot(1)!   ntot(2)!    ntot(3)!    ntot(4)!        1
!       -------- . --------- . --------- . --------- . ---------
!        nall!     ncell(1)!   ncell(2)!   ncell(3)!   ncell(4)!

!     The order of multiplications and divisions minimizes the chance
!     of overflows or underflows without using logarithms.

prob = one
DO
  ended = .true.
  IF (nall > ntot(1)) THEN
    prob = prob / nall
    nall = nall - 1
    ended = .false.
  END IF

  IF (ntot(2) > ncell(1)) THEN
    prob = prob * ntot(2)
    ntot(2) = ntot(2) - 1
    ended = .false.
  END IF

  IF (ntot(3) > ncell(2)) THEN
    prob = prob * ntot(3)
    ntot(3) = ntot(3) - 1
    ended = .false.
  END IF

  IF (ntot(4) > ncell(3)) THEN
    prob = prob * ntot(4)
    ntot(4) = ntot(4) - 1
    ended = .false.
  END IF

  IF (ncell(4) > 1) THEN
    prob = prob / ncell(4)
    ncell(4) = ncell(4) - 1
    ended = .false.
  END IF
  IF (ended) EXIT
END DO

RETURN
END FUNCTION hyperg



SUBROUTINE order4(n)

!     Order 4 integers; largest first.

IMPLICIT NONE
INTEGER, INTENT(IN OUT) :: n(4)

!     Local variable
INTEGER :: nn

IF (n(1) < n(2)) THEN
  nn = n(1)
  n(1) = n(2)
  n(2) = nn
END IF

IF (n(3) < n(4)) THEN
  nn = n(3)
  n(3) = n(4)
  n(4) = nn
END IF

IF (n(2) < n(3)) THEN
  nn = n(2)
  n(2) = n(3)
  n(3) = nn
  IF (n(3) < n(4)) THEN
    nn = n(3)
    n(3) = n(4)
    n(4) = nn
  END IF
  IF (n(1) < n(2)) THEN
    nn = n(1)
    n(1) = n(2)
    n(2) = nn
    IF (n(2) < n(3)) THEN
      nn = n(2)
      n(2) = n(3)
      n(3) = nn
    END IF
  END IF
END IF

RETURN
END SUBROUTINE order4



SUBROUTINE alias(p, n, f, l)

! Set up the aliases for the method described in:

! Kronmal, R.A. & Peterson, A.V. (1979)
! `On the Alias Method for Generating Random Variables from a Discrete
! Distribution', The American Statistician, vol.33(4), 214-218.

! The algorithm improves slightly on the basic algorithm given in:

! Walker, A.J. (1977)
! `An Efficient Method for Generating Discrete Random Variables with General
! Distributions', ACM Transactions on Mathematical Software, vol.3, 253-256.

! Arguments:
! p(1:n)  INPUT   Contains the discrete distribution
!                 N.B. It must have a finite range
! n       INPUT
! f(1:n)  OUTPUT  The cut-off values for each value of the random variable
! l(1:n)  OUTPUT  The aliases (at least one will be undefined, in which case
!                 it is not needed as the corresponding f(i) >= 1.

REAL (dp), INTENT(IN)   :: p(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: f(:)
INTEGER, INTENT(OUT)    :: l(:)

! Local variables

INTEGER    :: i, ig, is, gptr(n), j, k, sptr(n-1)
REAL (dp), PARAMETER  :: nearly_one = 1.0_dp - 100.0_dp*EPSILON(1.0_dp)

! Step 1
f(1:n) = n * p(1:n)

! Step 2
! G = set of i such that f(i) >= 1
! S = the others
ig = 0
is = 0
DO i = 1, n
  IF (f(i) < 1.0) THEN
    is = is + 1
    sptr(is) = i
  ELSE
    ig = ig + 1
    gptr(ig) = i
  END IF
END DO

! Step 3
! Main loop
DO
  IF (is == 0) EXIT

! Step 4
! Choose element k in G, and element j in S.
! We choose the last ones on each list.
  k = gptr(ig)
  j = sptr(is)
! This finalizes operations on element j, so we can remove it from S.
  is = is - 1

! Step 5
  l(j) = k

! Step 6
  f(k) = f(k) - (1.0 - f(j))

! Step 7
! f(k) may now be < 1.
! If so, remove it from G and add it to S.
  IF (f(k) < nearly_one) THEN
    ig = ig - 1
    is = is + 1
    sptr(is) = k
  END IF
END DO

RETURN
END SUBROUTINE alias

END MODULE random_hypergeometric



MODULE toms668
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, SAVE       :: iseed = 123456789

CONTAINS


SUBROUTINE h2pec(kk, nn1, nn2, iseed1, jx)

!  HYPERGEOMETRIC RANDOM VARIATE GENERATOR

!  METHOD
!     IF (MODE - MAX(0,KK-NN2) .LT. 10), USE THE INVERSE CDF.
!        OTHERWISE, USE ALGORITHM H2PE: ACCEPTANCE-REJECTION VIA
!        THREE REGION COMPOSITION.  THE THREE REGIONS ARE A
!        RECTANGLE, AND EXPONENTIAL LEFT AND RIGHT TAILS.
!     H2PE  REFERS TO HYPERGEOMETRIC-2 POINTS-EXPONENTIAL TAILS.
!     H2PEC REFERS TO H2PE AND "COMBINED."  THUS H2PE IS THE RESEARCH RESULT
!        AND H2PEC IS THE IMPLEMENTATION OF A COMPLETE USABLE ALGORITHM.

!  REFERENCE
!     VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,

!     "COMPUTER GENERATION OF HYPERGEOMETRIC RANDOM VARIATES,"
!     JOURNAL OF STATISTICAL COMPUTATION AND SIMULATION,
!     22(1985), 2, 1985, 127-145.

!  REQUIRED SUBPROGRAMS
!     AFC() : A DOUBLE-PRECISION FUNCTION TO EVALUATE
!                THE LOGARITHM OF THE FACTORIAL.
!     RAND(): A UNIFORM (0,1) RANDOM NUMBER GENERATOR.

!  ARGUMENTS
!     NN1   : NUMBER OF WHITE BALLS          (INPUT)
!     NN2   : NUMBER OF BLACK BALLS          (INPUT)
!     KK    : NUMBER OF BALLS TO BE DRAWN    (INPUT)
!     ISEED : RANDOM NUMBER SEED  (INPUT AND OUTPUT)
!     JX    : NUMBER OF WHITE BALLS DRAWN   (OUTPUT)

!  STRUCTURAL VARIABLES
!     REJECT: LOGICAL FLAG TO REJECT THE VARIATE GENERATE BY H2PE.
!     SETUP1: LOGICAL FLAG TO SETUP FOR NEW VALUES OF NN1 OR NN2.
!     SETUP2: LOGICAL FLAG TO SETUP FOR NEW VALUES OF KK.
!     IX    : INTEGER CANDIDATE VALUE.
!     M     : DISTRIBUTION MODE.
!     MINJX : DISTRIBUTION LOWER BOUND.
!     MAXJX : DISTRIBUTION UPPER BOUND.
!     KS    : SAVED VALUE OF KK FROM THE LAST CALL TO H2PEC.
!     N1S   : SAVED VALUE OF NN1 FROM THE LAST CALL TO H2PEC.
!     N2S   : SAVED VALUE OF NN2 FROM THE LAST CALL TO H2PEC.
!     K,N1,N2: ALTERNATE VARIABLES FOR KK, NN1, AND NN2
!                (ALWAYS (N1 .LE. N2) AND (K .LE. (N1+N2)/2)).
!     TN    : TOTAL NUMBER OF WHITE AND BLACK BALLS

!  INVERSE-TRANSFORMATION VARIABLES
!     CON   : NATURAL LOGARITHM  OF SCALE.
!     P     : CURRENT SCALED PROBABILITY FOR THE INVERSE CDF.
!     SCALE : A BIG CONSTANT (1.E25) USED TO SCALE THE
!                PROBABILITY TO AVOID NUMERICAL UNDERFLOW
!     U     : THE UNIFORM VARIATE BETWEEN (0, 1.E25).
!     W     : SCALED HYPERGEOMETRIC PROBABILITY OF MINJX.

!  H2PE VARIABLES
!     S     : DISTRIBUTION STANDARD DEVIATION.
!     D     : HALF THE AREA OF THE RECTANGLE.
!     XL    : LEFT END OF THE RECTANGLE.
!     XR    : RIGHT END OF THE RECTANGLE.
!     A     : A SCALING CONSTANT.
!     KL    : HIGHEST POINT OF THE LEFT-TAIL REGION.
!     KR    : HIGHEST POINT OF THE RIGHT-TAIL REGION.
!     LAMDL : RATE FOR THE LEFT EXPONENTIAL TAIL.
!     LAMDR : RATE FOR THE RIGHT EXPONENTIAL TAIL.
!     P1    : AREA OF THE RECTANGLE.
!     P2    : AREA OF THE LEFT EXPONENTIAL TAIL PLUS P1.
!     P3    : AREA OF THE RIGHT EXPONENTIAL TAIL PLUS P2.
!     U     : A UNIFORM (0,P3) RANDOM VARIATE USED FIRST TO SELECT
!                ONE OF THE THREE REGIONS AND THEN CONDITIONALLY TO
!                GENERATE A VALUE FROM THE REGION.
!     V     : U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM
!                VALUE OR TO ACCEPT OR REJECT THE CANDIDATE VALUE.
!     F     : THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
!                ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL.
!     I     : INDEX FOR EXPLICIT CALCULATION OF F FOR H2PE.

! THE FOLLOWING VARIABLES ARE TEMPORARY VARIABLES USED IN COMPUTING THE UPPER
! AND LOWER BOUNDS OF THE NATURAL LOGARITHM OF THE SCALED DENSITY.
! THE DETAILED DESCRIPTION IS GIVEN IN PROPOSITIONS 2 AND 3 OF THE APPENDIX IN
! THE REFERENCE.
!     Y, Y1, YM, YN, YK, NK, R, S, T, E, G, DG, GU, GL, XM, XN, XK, NM

!     Y     : PRELIMINARY CONTINUOUS CANDIDATE VALUE, FLOAT(IX)
!     UB    : UPPER BOUND FOR THE NATURAL LOGARITHM OF THE SCALED DENSITY.
!     ALV   : NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V.
!     DR, DS, DT, DE: ONE OF MANY TERMS SUBTRACTED FROM THE UPPER
!                BOUND TO OBTAIN THE LOWER BOUND ON THE NATURAL
!                LOGARITHM OF THE SCALED DENSITY.
!     DELTAU: A CONSTANT, THE VALUE 0.0034 IS OBTAINED BY SETTING
!                N1 = N2 = 200, K = 199, M = 100, AND Y = 50 IN
!                THE FUNCTION DELTA_U IN LEMMA 1 AND ROUNDING THE
!                VALUE TO FOUR DECIMAL PLACES.
!     DELTAL: A CONSTANT, THE VALUE 0.0078 IS OBTAINED BY SETTING
!                N1 = N2 = 200, K = 199, M = 100, AND Y = 50 IN
!                THE FUNCTION DELTA_L IN LEMMA 1 AND ROUNDING THE
!                VALUE TO FOUR DECIMAL PLACES.

INTEGER, INTENT(IN)      :: kk
INTEGER, INTENT(IN)      :: nn1
INTEGER, INTENT(IN)      :: nn2
INTEGER, INTENT(IN OUT)  :: iseed1
INTEGER, INTENT(OUT)     :: jx

REAL (dp)      :: p, u, w, a, xl, xr
REAL           :: alv, d, dg, de, dr, ds, dt, e, f, g, gl, gu, kl, kr,  &
                  lamdl, lamdr, nk, nm, p1, p2, p3, r, s, t, ub, v,  &
                  xk, xm, xn, y, y1, yk, ym, yn
REAL, SAVE     :: tn
INTEGER        :: i, ix, k, n1, n2
LOGICAL, SAVE  :: reject, setup1, setup2
INTEGER, SAVE  :: ks = -1, n1s = -1, n2s = -1, m, minjx, maxjx
REAL (dp), PARAMETER  :: con = 57.56462733_dp, deltal = 0.0078_dp,  &
                         deltau = 0.0034_dp, scale = 1.0e+25_dp

!*****CHECK PARAMETER VALIDITY

IF ( nn1 < 0  .OR. nn2 < 0 .OR. kk  < 0 .OR. kk > nn1 + nn2 ) THEN
  jx     = -1
  RETURN
END IF
iseed = iseed1

!*****IF NEW PARAMETER VALUES, INITIALIZE

reject = .true.
setup1 = .false.
setup2 = .false.
IF (nn1 /= n1s .OR. nn2 /= n2s)  THEN
  setup1 = .true.
  setup2 = .true.
ELSE IF (kk /= ks)  THEN
  setup2 = .true.
END IF

IF (setup1)  THEN
  n1s   = nn1
  n2s   = nn2
  tn    = nn1 + nn2
  IF (nn1 <= nn2)  THEN
    n1 = nn1
    n2 = nn2
  ELSE
    n1 = nn2
    n2 = nn1
  END IF
END IF

IF (setup2)  THEN
  ks    = kk
  IF (kk+kk >= tn)  THEN
    k  = tn - kk
  ELSE
    k  = kk
  END IF
END IF

IF (setup1 .OR. setup2)  THEN
  m     = (k+1.) * (n1+1.) / (tn+2.)
  minjx = MAX (0, k-n2)
  maxjx = MIN (n1, k)
END IF

!*****GENERATE RANDOM VARIATE

IF (minjx == maxjx)  THEN
  
!        ...DEGENERATE DISTRIBUTION...
  
  ix      = maxjx
  RETURN
ELSE IF (m-minjx < 10)  THEN
  
!        ...INVERSE TRANSFORMATION...
  
  IF (setup1 .OR. setup2)  THEN
    IF (k < n2) THEN
      w = EXP (con + afc(n2) + afc(n1+n2-k) - afc(n2-k) - afc(n1+n2))
    ELSE
      w = EXP (con + afc(n1) + afc(k) - afc(k-n2) - afc(n1+n2))
    END IF
  END IF
  
  10 p  = w
  ix = minjx
  u  = rand () * scale
  20 IF (u > p)  THEN
    u  = u - p
    p  = p * (n1-ix)*(k-ix)
    ix = ix + 1
    p  = p / ix / (n2-k+ix)
    IF (ix > maxjx)  GO TO 10
    GO TO 20
  END IF
ELSE
  
!        ...H2PE...
  
  IF (setup1 .OR. setup2)  THEN
    s     = SQRT ((tn-k) * k * n1 * n2 / (tn-1) / tn /tn)
    
!           ...REMARK:  D IS DEFINED IN REFERENCE WITHOUT INT.
!           THE TRUNCATION CENTERS THE CELL BOUNDARIES AT 0.5
    
    d     = INT (1.5*s) + .5
    xl    = m - d + .5
    xr    = m + d + .5
    a     = afc(m) + afc(n1-m) + afc(k-m) + afc(n2-k+m)
    kl    = EXP (a - afc(INT(xl)) - afc(INT(n1-xl))  &
            - afc(INT(k-xl)) - afc(INT(n2-k+xl)))
    kr    = EXP (a - afc(INT(xr-1)) - afc(INT(n1-xr+1))  &
            - afc(INT(k-xr+1)) - afc(INT(n2-k+xr-1)))
    lamdl = -LOG (xl * (n2-k+xl) / (n1-xl+1) / (k-xl+1))
    lamdr = -LOG ((n1-xr+1) * (k-xr+1) / xr / (n2-k+xr))
    p1    = d + d
    p2    = p1 + kl / lamdl
    p3    = p2 + kr / lamdr
  END IF
  
  30 u    = rand () * p3
  v       = rand ()
  IF (u < p1)  THEN
    
!           ...RECTANGULAR REGION...
    
    ix    = xl + u
  ELSE IF (u <= p2)  THEN
    
!           ...LEFT TAIL...
    
    ix    = xl + LOG(v)/lamdl
    IF (ix < minjx)  GO TO 30
    v     = v * (u-p1) * lamdl
  ELSE
    
!           ...RIGHT TAIL...
    
    ix    = xr - LOG(v)/lamdr
    IF (ix > maxjx)  GO TO 30
    v     = v * (u-p2) * lamdr
  END IF
  
!        ...ACCEPTANCE/REJECTION TEST...
  
  IF (m < 100 .OR. ix <= 50)  THEN
    
!           ...EXPLICIT EVALUATION...
    
    f     = 1.0
    IF (m < ix)  THEN
      DO  i = m+1,ix
        f   = f * (n1-i+1) * (k-i+1) / (n2-k+i) / i
      END DO
    ELSE IF (m > ix)  THEN
      DO  i = ix+1,m
        f   = f * i * (n2-k+i) / (n1-i) / (k-i)
      END DO
    END IF
    IF (v <= f)  THEN
      reject = .false.
    END IF
  ELSE
    
!        ...SQUEEZE USING UPPER AND LOWER BOUNDS...
    
    y   = ix
    y1  = y + 1.
    ym  = y - m
    yn  = n1 - y + 1.
    yk  = k - y + 1.
    nk  = n2 - k + y1
    r   = -ym / y1
    s   = ym / yn
    t   = ym / yk
    e   = -ym / nk
    g   = yn * yk / (y1*nk) - 1.
    dg  = 1.
    IF (g < 0.) dg = 1.+g
    gu  = g * (1.+g*(-.5+g/3.))
    gl  = gu - .25 * (g*g)**2 / dg
    xm  = m + .5
    xn  = n1 - m + .5
    xk  = k - m + .5
    nm  = n2 - k + xm
    ub  = y * gu - m * gl + deltau + xm * r * (1.+r*(-.5+r/3.))  &
          + xn * s * (1.+s*(-.5+s/3.)) + xk * t * (1.+t*(-.5+t/3.))  &
          + nm * e * (1.+e*(-.5+e/3.))
    
!           ...TEST AGAINST UPPER BOUND...
    
    alv = LOG(v)
    IF (alv > ub)  THEN
      reject = .true.
    ELSE
      
!              ...TEST AGAINST LOWER BOUND...
      
      dr = xm * (r*r)**2
      IF (r < 0.)  dr = dr / (1.+r)
      ds = xn * (s*s)**2
      IF (s < 0.)  ds = ds / (1.+s)
      dt = xk * (t*t)**2
      IF (t < 0.)  dt = dt / (1.+t)
      de = nm * (e*e)**2
      IF (e < 0.)  de = de / (1.+e)
      IF (alv < ub-.25*(dr+ds+dt+de) +(y+m)*(gl-gu)-deltal) THEN
        reject = .false.
      ELSE
        
!                 ...STIRLING'S FORMULA TO MACHINE ACCURACY...
        
        IF (alv <= a - afc(ix) - afc(n1-ix) - afc(k-ix) - afc(n2-k+ix) ) THEN
          reject = .false.
        ELSE
          reject = .true.
        END IF
      END IF
    END IF
  END IF
  IF (reject)  GO TO 30
END IF

!*****RETURN APPROPRIATE VARIATE

IF (kk + kk >= tn)  THEN
  IF (nn1 > nn2)  THEN
    ix = kk - nn2 + ix
  ELSE
    ix =  nn1 - ix
  END IF
ELSE
  IF (nn1 > nn2)  ix = kk - ix
END IF
jx = ix

RETURN
END SUBROUTINE h2pec



FUNCTION afc(i) RESULT(fn_val)

!     FUNCTION TO EVALUATE LOGARITHM OF THE FACTORIAL I
!        IF (I .GT. 7), USE STIRLING'S APPROXIMATION
!        (It actually uses the less accurate De Moivre approximation)
!           OTHERWISE,  USE TABLE LOOKUP

INTEGER, INTENT(IN)  :: i
REAL (dp)            :: fn_val

REAL (dp) :: di

REAL (dp), PARAMETER  :: al(8) = (/ 0._dp, 0._dp, 0.6931471806_dp,  &
                         1.791759469_dp, 3.178053830_dp, 4.787491743_dp, &
                         6.579251212_dp, 8.525161361_dp /)

IF (i <= 7)  THEN
  fn_val = al(i+1)
ELSE
  di  = i
  fn_val = (di+0.5_dp) * LOG(di) - di + 0.08333333333333_dp/di  &
           - 0.00277777777777_dp/di/di/di  + 0.9189385332_dp
END IF

RETURN
END FUNCTION afc



FUNCTION rand() RESULT(fn_val)

!  UNIFORM RANDOM NUMBER GENERATOR
!  REFERENCE:  L. SCHRAGE,
!     "A MORE PORTABLE FORTRAN RANDOM NUMBER GENERATOR,"
!     ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE, 5(1979), 132-138.

REAL                :: fn_val

INTEGER             :: xhi, xalo, leftlo, fhi, k
INTEGER, PARAMETER  :: a = 16807, b15 = 32768, b16 = 65536, p = 2147483647

xhi    = iseed / b16
xalo   = (iseed-xhi * b16) * a
leftlo = xalo/b16
fhi    = xhi * a + leftlo
k      = fhi / b15
iseed  = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k
IF (iseed < 0) iseed = iseed + p
fn_val = iseed * 4.656612875E-10

RETURN
END FUNCTION rand

END MODULE toms668




PROGRAM test_random_hyperg
! TEST 1
! ------
! Test the random generation of hypergeometric variables by generating
! the number of spades in a hand of 13 cards out of 52.
! 1 million replicates used.

USE random_hypergeometric
USE toms668
IMPLICIT NONE

REAL     :: t1, t2
INTEGER  :: freq(0:13), c, d, h, i, j, s

WRITE(*, *) '  TEST 1'
WRITE(*, *)

! First using r_hyperg
freq = 0
CALL CPU_TIME(t1)
DO i = 1, 1000000
  j = r_hyperg(52, 13, 13)
  freq(j) = freq(j) + 1
END DO
CALL CPU_TIME(t2)
WRITE(*, '(a, f8.2, a)') ' Time taken by r_hyperg = ', t2-t1, 'secs.'
WRITE(*, '(a/9i7/5i7)') ' Frequencies', freq

! Now using h_alias
freq = 0
CALL CPU_TIME(t1)
DO i = 1, 1000000
  j = h_alias(52, 13, 13)
  freq(j) = freq(j) + 1
END DO
CALL CPU_TIME(t2)
WRITE(*, '(a, f8.2, a)') ' Time taken by h_alias = ', t2-t1, 'secs.'
WRITE(*, '(a/9i7/5i7)') ' Frequencies', freq

! Now using TOMS algorithm 668

freq = 0
CALL CPU_TIME(t1)
DO i = 1, 1000000
  CALL h2pec(13, 13, 39, iseed, j)
  freq(j) = freq(j) + 1
END DO
CALL CPU_TIME(t2)
WRITE(*, '(a, f8.2, a)') ' Time taken by TOMS668 = ', t2-t1, 'secs.'
WRITE(*, '(a/9i7/5i7)') ' Frequencies', freq

! TEST 2
! ------
! This tests the set-up time as the parameters change on every call.
! The number of clubs, diamonds, hearts and spades are found in a hand
! of 13 cards out of 52 by calling the random number generator 3 times
! with different parameters.

WRITE(*, *)
WRITE(*, *) '  TEST 2'
WRITE(*, *)

! First using r_hyperg
freq = 0
CALL CPU_TIME(t1)
DO i = 1, 100000
  c = r_hyperg(52, 13, 13)
  d = r_hyperg(39, 13, 13-c)
  h = r_hyperg(26, 13, 13-c-d)
  s = 13 - c - d - h
  freq(c) = freq(c) + 1
  freq(d) = freq(d) + 1
  freq(h) = freq(h) + 1
  freq(s) = freq(s) + 1
END DO
CALL CPU_TIME(t2)
WRITE(*, '(a, f8.2, a)') ' Time taken by r_hyperg = ', t2-t1, 'secs.'
WRITE(*, '(a/9i7/5i7)') ' Frequencies', freq

! Now using h_alias

freq = 0
CALL CPU_TIME(t1)
DO i = 1, 100000
  c = h_alias(52, 13, 13)
  d = h_alias(39, 13, 13-c)
  h = h_alias(26, 13, 13-c-d)
  s = 13 - c - d - h
  freq(c) = freq(c) + 1
  freq(d) = freq(d) + 1
  freq(h) = freq(h) + 1
  freq(s) = freq(s) + 1
END DO
CALL CPU_TIME(t2)
WRITE(*, '(a, f8.2, a)') ' Time taken by h_alias = ', t2-t1, 'secs.'
WRITE(*, '(a/9i7/5i7)') ' Frequencies', freq

! Now using TOMS algorithm 668

freq = 0
CALL CPU_TIME(t1)
DO i = 1, 100000
  CALL h2pec(13, 13, 39, iseed, c)
  CALL h2pec(13-c, 13, 26, iseed, d)
  CALL h2pec(13-c-d, 13, 13, iseed, h)
  s = 13 - c - d - h
  freq(c) = freq(c) + 1
  freq(d) = freq(d) + 1
  freq(h) = freq(h) + 1
  freq(s) = freq(s) + 1
END DO
CALL CPU_TIME(t2)
WRITE(*, '(a, f8.2, a)') ' Time taken by TOMS668 = ', t2-t1, 'secs.'
WRITE(*, '(a/9i7/5i7)') ' Frequencies', freq

STOP
END PROGRAM test_random_hyperg
