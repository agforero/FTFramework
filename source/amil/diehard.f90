MODULE Diehard_Tests

! This is translated from George Marsaglia's DIEHARD package
! for testing random INTEGERS.

! Code converted using TO_F90 by Alan Miller
! Date: 2001-02-24  Time: 17:32:26

! Latest revision - 14 January 2004
! Corrections by Michael Bush incorporated 19 September 2001

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, SAVE  :: jkk = 1    ! Reset by function jkreset

CONTAINS


! First, some functions which were inline.

FUNCTION uni() RESULT(fn_val)
REAL  :: fn_val

! Put a one-line function here to provide the uni being tested:
fn_val = .5 + jtbl() * .5 ** 32

RETURN
END FUNCTION uni



FUNCTION ch(x) RESULT(fn_val)
REAL, INTENT(IN)  :: x
REAL              :: fn_val

!**** up and down runs test******************
fn_val = 1. - EXP(-.5*x) * (1. + .5*x + .125*x**2)

RETURN
END FUNCTION ch



FUNCTION kthrow() RESULT(fn_val)
INTEGER  :: fn_val

REAL, PARAMETER    :: cc = 6.*.5**32

fn_val = 2 + INT(cc*jtbl()+3.) + INT(cc*jtbl()+3.)

RETURN
END FUNCTION kthrow



FUNCTION ikbit(kr, mk) RESULT(fn_val)
INTEGER, INTENT(IN)  :: kr, mk
INTEGER              :: fn_val

!** ONE-LINE FUNCTION TO GENERATE 5-BIT LETTER IN CONVENIENT POSITION
fn_val = IAND(ISHFT(jtbl(),-kr), mk)

RETURN
END FUNCTION ikbit



SUBROUTINE sqsort(a, n, t)
 
!   NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!   BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

!   SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.

REAL, INTENT(IN OUT)  :: a(:)
INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: t(:)

INTEGER :: i, j, k, l, r, s, stackl(15), stackr(15)
REAL    :: w, ww, x

s = 1
stackl(1) = 1
stackr(1) = n

!     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10 l = stackl(s)
r = stackr(s)
s = s - 1

!     KEEP SPLITTING A(L), ... , A(R) UNTIL L>=R.

20 i = l
j = r
k = (l+r) / 2
x = a(k)

!     REPEAT UNTIL I > J.

30 DO
  IF (a(i) >= x) EXIT
  i = i + 1
END DO

50 IF (x < a(j)) THEN
  j = j - 1
  GO TO 50
END IF

IF (i <= j) THEN
  w = a(i)
  ww = t(i)
  a(i) = a(j)
  t(i) = t(j)
  a(j) = w
  t(j) = ww
  i = i + 1
  j = j - 1
  IF (i <= j) GO TO 30
END IF

IF (j-l >= r-i) THEN
  IF (l < j) THEN
    s = s + 1
    stackl(s) = l
    stackr(s) = j
  END IF
  l = i
ELSE
  IF (i < r) THEN
    s = s + 1
    stackl(s) = i
    stackr(s) = r
  END IF
  r = j
END IF
IF (l < r) GO TO 20
IF (s /= 0) GO TO 10

RETURN
END SUBROUTINE sqsort



! Now for the tests

SUBROUTINE cdbitst(filename)
!     This is test 5

!     THE OVERLAPPING 20-tuples TEST  BITSTREAM, 20 BITS PER WORD, N words
!     If n=2^22, should be 19205.3 missing 20-letter words, sigma 167.
!     If n=2^21, should be 141909  missing 20-letter words, sigma 428.
!     If n=2^20, should be 385750  missing 20-letter words, sigma 512

CHARACTER (LEN=25), INTENT(IN)  :: filename

INTEGER            :: w(0:32767), mbit(0:31)
CHARACTER (LEN=80) :: text(36), dum
REAL               :: mu, sigma, x
INTEGER            :: i, ib, ic, j, k, kount, kpow, l, nt, ntries, num

jkk = jkreset()
OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:16)
WRITE (*, 5000) text(1:16)
WRITE (3,5000) text(1:16)
CLOSE (4)

kpow = 21
ntries = 20
sigma = 428
!**** SET MASK BITS***************

      DO  i = 0, 31
        mbit(i) = ishft(1,i)
      END DO
!***********INITIALIZE*****************
WRITE (3,*) '      THE OVERLAPPING 20-tuples BITSTREAM TEST,'
WRITE (3,*) '           20 BITS PER WORD, 2^21 words.'
WRITE (*, *) '      THE OVERLAPPING 20-tuples BITSTREAM TEST,'
WRITE (*, *) '           20 BITS PER WORD, 2^21 words.'
WRITE (3,*) '   This test samples the bitstream 20 times.'
WRITE (*, *) '   This test samples the bitstream 20 times.'
mu = 2 ** 20 * EXP(-2.**(kpow-20))
WRITE (3,5100) mu, sigma
WRITE (*, 5100) mu, sigma

!*****MAIN LOOP*********
!**** GET INITIAL WORD
j = jtbl()
j = IAND(j, 2**20-1)
WRITE (*, 5200) filename

DO  nt = 1, ntries
!     ********SET W-TABLE TO ZEROS*******
  w(0:32767) = 0
!**** GENERATE 2**kpow OVERLAPPING WORDS**********
  DO  ic = 1, 2 ** (kpow-5)
    num = jtbl()
    DO  ib = 1, 32
!     *** GET NEW J *****
      j = ISHFT(IAND(j,2**19-1), 1) + IAND(num,1)
      num = ISHFT(num, -1)
!     *** GET BIT INDEX FROM LAST 5 BITS OF J  ***
      l = IAND(j,31)
!     *** GET TABLE INDEX FROM LEADING 15 BITS OF J***
      k = ISHFT(j, -5)
!     *** SET BIT L IN W(K) ***
      w(k) = IOR(w(k),mbit(l))
    END DO
  END DO
!     ********** COUNT NUMBER OF EMPTY CELLS *******
  kount = 0
  DO  k = 0, 32767
    DO  l = 0, 31
      IF (IAND(w(k),mbit(l)) == 0) kount = kount + 1
    END DO
  END DO
!     ****END OF MAIN LOOP****
  x = (kount - mu) / sigma
  WRITE (3,5300) nt, kount, x, phi(x)
  WRITE (*, 5300) nt, kount, x, phi(x)
END DO

jkk = jkreset()
WRITE (*, 5400)
WRITE (3,5400)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('  No. missing words should average', f9.0, ' with sigma= ',  &
             f4.0/ '--------------------------------------------------')
5200 FORMAT (' BITSTREAM test results for', a15)
5300 FORMAT (' tst no ',i2,': ',i7,' missing words, ', f7.2,  &
             ' sigmas from mean, p-value=', f7.5)
5400 FORMAT (/'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE cdbitst



SUBROUTINE d3sphere(filename)
!     This is test 11

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::              THE 3DSPHERES TEST                               ::
!     :: Choose  4000 random points in a cube of edge 1000.  At each   ::
!     :: point, center a sphere large enough to reach the next closest ::
!     :: point. Then the volume of the smallest such sphere is (very   ::
!     :: close to) exponentially distributed with mean 120pi/3.  Thus  ::
!     :: the radius cubed is exponential with mean 30. (The mean is    ::
!     :: obtained by extensive simulation).  The 3DSPHERES test gener- ::
!     :: ates 4000 such spheres 20 times.  Each min radius cubed leads ::
!     :: to a uniform variable by means of 1-exp(-r^3/30.), then a     ::
!     ::  KSTEST is done on the 20 p-values.                           ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


CHARACTER (LEN=25), INTENT(IN)  :: filename

REAL                 :: x(4000), y(4000), z(4000), p(20)
CHARACTER (LEN=80)   :: text(36)
CHARACTER (LEN=186)  :: dum
INTEGER              :: i, ij, j, n
REAL                 :: d, dmin, pv, r3, u, v, w

n = 4000
OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4, 5000) dum
READ (4, 5000) text(1:12)
WRITE (*, 5000) text(1:12)
WRITE (3, 5000) text(1:12)
CLOSE (4)

jkk = jkreset()
WRITE (3, 5100) filename
WRITE (*, 5100) filename
DO  ij = 1, 20
  dmin = 10000000.
  DO  i = 1, n
    x(i) = 500. + 0.2328306e-6 * jtbl()
  END DO
  CALL asort(x,n)
  DO  i = 1, n
    y(i) = 500. + 0.2328306e-6 * jtbl()
    z(i) = 500. + 0.2328306e-6 * jtbl()
  END DO
  loop40:  DO  i = 1, n
    u = x(i)
    v = y(i)
    w = z(i)
    DO  j = i + 1, n
      d = (u-x(j)) ** 2
      IF (d >= dmin) CYCLE loop40
      d = d + (v-y(j)) ** 2 + (w-z(j)) ** 2
      IF (d < dmin) dmin = d
    END DO
  END DO loop40
  r3 = dmin * SQRT(dmin)
  p(ij) = 1 - EXP(-r3/30.)
  WRITE (3, 5200) ij, r3, p(ij)
  WRITE (*, 5200) ij, r3, p(ij)
END DO
WRITE (3, 5300)
WRITE (*, 5300)
CALL kstest(p,20,pv)
WRITE (*, 5400) filename, pv
WRITE (3, 5400) filename, pv
jkk = jkreset()
RETURN

5000 FORMAT (a78)
5100 FORMAT ('               The 3DSPHERES test for file ', a15)
5200 FORMAT (' sample no: ', i2, '     r^3= ', f7.3, '     p-value= ', f7.5)
5300 FORMAT ('  A KS test is applied to those 20 p-values.'/   &
             '---------------------------------------------------------')
5400 FORMAT ('       3DSPHERES test for file ', a15, '      p-value= ', f8.6/  &
             '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE d3sphere



FUNCTION jtbl8() RESULT(fn_val)
INTEGER        :: fn_val

INTEGER, SAVE  :: b(4096)
INTEGER, SAVE  :: nleft = 0, j = 4097, jk = 1
INTEGER, SAVE  :: k

IF (j > 4096) THEN
  READ (1, REC=jk) b(1:4096)
  j = 1
  jk = jk + 1
END IF

IF (nleft == 0) THEN
  k = b(j)
  j = j + 1
  nleft = 4
END IF

fn_val = IAND(ISHFT(k, -24), 255)
k = ISHFT(k, 8)
nleft = nleft - 1

RETURN
END FUNCTION jtbl8



SUBROUTINE sqeez(filename)
!  This is test 12

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::      This is the SQEEZE test                                  ::
!     ::  Random integers are floated to get uniforms on [0,1). Start- ::
!     ::  ing with k=2^31=2147483647, the test finds j, the number of  ::
!     ::  iterations necessary to reduce k to 1, using the reduction   ::
!     ::  k=ceiling(k*U), with U provided by floating integers from    ::
!     ::  the file being tested.  Such j's are found 100,000 times,    ::
!     ::  then counts for the number of times j was <=6,7,...,47,>=48  ::
!     ::  are used to provide a chi-square test for cell frequencies.  ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!  SQUEEZE TEST.  How many iterations of k=k*uni()+1 are required
!  to squeeze k down to 1, starting with k=2147483647=2^31-1.
!  The exact distribution of the required j is used, with
!  a chi-square test based on 100,000 tries.
!  The mean of j is 23.064779, with variance 23.70971151.


CHARACTER (LEN=25), INTENT(IN)  :: filename

REAL                 :: chsq, sig, tbl(6:48)
CHARACTER (LEN=80)   :: text(10)
CHARACTER (LEN=198)  :: dum
INTEGER              :: i, j, k

REAL, PARAMETER  :: ex(6:48) = (/  &
    21.03, 57.79, 175.54, 467.32, 1107.83,  &
    2367.84, 4609.44, 8241.16, 13627.81, 20968.49, 30176.12,  &
    40801.97, 52042.03, 62838.28, 72056.37, 78694.51, 82067.55,  &
    81919.35, 78440.08, 72194.12, 63986.79, 54709.31, 45198.52,  &
    36136.61, 28000.28, 21055.67, 15386.52, 10940.20, 7577.96,  &
    5119.56, 3377.26, 2177.87, 1374.39, 849.70, 515.18, 306.66,  &
    179.39, 103.24, 58.51, 32.69, 18.03, 9.82, 11.21 /)

OPEN (4, FILE='tests.txt', STATUS='OLD')
jkk = jkreset()
READ (4,5000) dum
READ (4,5000) text(1:10)
WRITE (*, 5000) text(1:10)
WRITE (3,5000) text(1:10)
CLOSE (4)

tbl(6:48) = 0
DO  i = 1, 100000
  j = 0
  k = 2147483647
  20 k = k * uni() + 1
  j = j + 1
  IF (k > 1) GO TO 20
  j = MIN(MAX(j,6),48)
  tbl(j) = tbl(j) + 1.
END DO
chsq = 0
DO  i = 6, 48
  chsq = chsq + (tbl(i)-.1*ex(i)) ** 2 / (.1*ex(i))
END DO
sig = SQRT(84.)
WRITE (3,5100) filename
WRITE (*, 5100) filename
WRITE (3,5200)
WRITE (*, 5200)
WRITE (3,5300) ((tbl(i)-.1*ex(i))/SQRT(.1*ex(i)),i = 6,48)
WRITE (*, 5300) ((tbl(i)-.1*ex(i))/SQRT(.1*ex(i)),i=6,48)
WRITE (3,5400) chsq
WRITE (*, 5400) chsq
WRITE (*, 5500) (chsq-42.) / sig, chisq(chsq,42)
WRITE (3,5500) (chsq-42.) / sig, chisq(chsq,42)
jkk = jkreset()
WRITE (*, 5600)
WRITE (3,5600)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('            RESULTS OF SQUEEZE TEST FOR ', a15)
5200 FORMAT ('         Table of standardized frequency counts'/   &
             '     ( (obs-exp)/sqrt(exp) )^2'/   &
             '        for j taking values <=6,7,8,...,47,>=48:')
5300 FORMAT (6F8.1)
5400 FORMAT (t9, '   Chi-square with 42 degrees of freedom:', f7.3)
5500 FORMAT (t9, '      z-score=', f7.3, '  p-value=', f8.6/   &
             '______________________________________________________________')
5600 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE sqeez



SUBROUTINE cdpark(filename)
!     This is test 9

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::               THIS IS A PARKING LOT TEST                      ::
!     :: In a square of side 100, randomly "park" a car---a circle of  ::
!     :: radius 1.   Then try to park a 2nd, a 3rd, and so on, each    ::
!     :: time parking "by ear".  That is, if an attempt to park a car  ::
!     :: causes a crash with one already parked, try again at a new    ::
!     :: random location. (To avoid path problems, consider parking    ::
!     :: helicopters rather than cars.)   Each attempt leads to either ::
!     :: a crash or a success, the latter followed by an increment to  ::
!     :: the list of cars already parked. If we plot n:  the number of ::
!     :: attempts, versus k::  the number successfully parked, we get a::
!     :: curve that should be similar to those provided by a perfect   ::
!     :: random number generator.  Theory for the behavior of such a   ::
!     :: random curve seems beyond reach, and as graphics displays are ::
!     :: not available for this battery of tests, a simple characteriz ::
!     :: ation of the random experiment is used: k, the number of cars ::
!     :: successfully parked after n=12,000 attempts. Simulation shows ::
!     :: that k should average 3523 with sigma 21.9 and is very close  ::
!     :: to normally distributed.  Thus (k-3523)/21.9 should be a st-  ::
!     :: andard normal variable, which, converted to a uniform varia-  ::
!     :: ble, provides input to a KSTEST based on a sample of 10.      ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


CHARACTER (LEN=25), INTENT(IN)  :: filename

REAL                 :: x(4000), y(4000), g(10)
CHARACTER (LEN=80)   :: text(22)
CHARACTER (LEN=151)  :: dum
INTEGER              :: i, ij, k, n
REAL                 :: av, pp, s, sig, ss, w, z

INTEGER, PARAMETER  :: ntries = 12000
REAL, PARAMETER     :: sq = 100.
INTEGER, PARAMETER  :: nt = 10

jkk = jkreset()
s = 0.
ss = 0.
OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:22)
WRITE (*, 5000) text(1:22)
WRITE (3,5000) text(1:22)
CLOSE (4)

WRITE (3,5100) filename
WRITE (*, 5100) filename

DO  ij = 1, nt
  x(1) = sq * uni()
  y(1) = sq * uni()
  k = 1
  loop20:  DO  n = 1, ntries
    z = sq * uni()
    w = sq * uni()
    DO  i = 1, k
      IF (ABS(x(i)-z) <= 1. .AND. ABS(y(i)-w) <= 1.) CYCLE loop20
    END DO
    k = k + 1
    x(k) = z
    y(k) = w
  END DO loop20
  s = s + k
  ss = ss + k * k
  z = (k-3523.) / 21.9
  g(ij) = phi(z)
  WRITE (*, 5200) k, z, g(ij)
  WRITE (3,5200) k, z, g(ij)
END DO
av = s / nt
sig = ss / nt - av ** 2
WRITE (3,*)
WRITE (*, 5300) sq, av, SQRT(sig)
WRITE (3,5300) sq, av, SQRT(sig)
CALL kstest(g,10,pp)
WRITE (*, 5400) pp
WRITE (3,5400) pp
jkk = jkreset()
WRITE (*, 5500)
WRITE (3,5500)
RETURN

5000 FORMAT (a78)
5100 FORMAT (t11, ' CDPARK: result of ten tests on file ', a15/    &
             t11, '  Of 12,000 tries, the average no. of successes'/   &
             t16, '  should be 3523 with sigma=21.9')
5200 FORMAT (t11, '  Successes:',i5,'    z-score:', f7.3, ' p-value: ', f8.6 )
5300 FORMAT (t11, ' square size   avg. no.  parked   sample sigma'/   &
             t11, f7.0, f20.3, f13.3)
5400 FORMAT ('            KSTEST for the above 10: p= ', f8.6)
5500 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE cdpark



SUBROUTINE mindist(filename)
!  This is test 10

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::               THE MINIMUM DISTANCE TEST                       ::
!     :: It does this 100 times::   choose n=8000 random points in a   ::
!     :: square of side 10000.  Find d, the minimum distance between   ::
!     :: the (n^2-n)/2 pairs of points.  If the points are truly inde- ::
!     :: pendent uniform, then d^2, the square of the minimum distance ::
!     :: should be (very close to) exponentially distributed with mean ::
!     :: .995 .  Thus 1-exp(-d^2/.995) should be uniform on [0,1) and  ::
!     :: a KSTEST on the resulting 100 values serves as a test of uni- ::
!     :: formity for random points in the square. Test numbers=0 mod 5 ::
!     :: are printed but the KSTEST is based on the full set of 100    ::
!     :: random choices of 8000 points in the 10000x10000 square.      ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!  Minimum distance^2 between n  random points(x(i),y(i)).
!  Mean is about .64 for 4000 points in a square of side 1000.
!  and .995 for 8000 points in a square of side 10000.
!  Since distance^2 is approximately exponential with mean .04,
!  1.-exp(-d^2/.04) should be uniform on [0,1).  Thus a KS test.

CHARACTER (LEN=25), INTENT(IN)  :: filename

REAL       :: x(8000), g(100), y(8000)
INTEGER    :: i, ij, j, ns
REAL       :: d, dmin, p, sum, u, v

CHARACTER (LEN=80)  :: text(13)
CHARACTER (LEN=173) :: dum
! EQUIVALENCE (qq(1),xy(1))

OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:13)
WRITE (*, 5000) text(1:13)
WRITE (3, 5000) text(1:13)
CLOSE (4)

ns = 100
jkk = jkreset()
WRITE (*, 5100) filename
WRITE (3, 5100) filename

sum = 0.0
DO  ij = 1, ns
  dmin = 10000000.
  DO  i = 1, 8000
    x(i) = 5000. + jtbl() * 0.2328306E-5
    y(i) = 5000. + jtbl() * 0.2328306E-5
  END DO

! There was a crazy call to DSORT here in the F77 version.
! The double precision numbers in array QQ were obtained by
! equivalencing QQ to two consecutive reals from array XY!
  CALL sqsort(x, 8000, y)

  DO  i = 1, 7999
    u = x(i)
    v = y(i)
    DO  j = i+1, 8000

      d = (u-x(j)) ** 2 + (v-y(j)) ** 2
      dmin = MIN(d,dmin)

    END DO
  END DO

  d = dmin
  sum = sum + d

  g(ij) = 1. - EXP(-dmin/.995)
  IF (MOD(ij,5) == 0) THEN
    WRITE (3,5200) ij, d, sum / ij, g(ij)
    WRITE (*, 5200) ij, d, sum / ij, g(ij)
  END IF
END DO
WRITE (3,5300) filename
WRITE (*, 5300) filename
WRITE (3,5400)
WRITE (*, 5400)
CALL kstest(g,ns,p)
WRITE (3,5500) p
WRITE (*, 5500) p
WRITE (3,5600)
WRITE (*, 5600)
jkk = jkreset()
RETURN

5000 FORMAT (a78)
5100 FORMAT ('               This is the MINIMUM DISTANCE test'/   &
             '              for random integers in the file ', a15/   &
             t5, ' Sample no.    d^2     avg     equiv uni            ')
5200 FORMAT (i12, f10.4, f9.4, f12.6)
5300 FORMAT ('     MINIMUM DISTANCE TEST for ', a15)
5400 FORMAT (t11, 'Result of KS test on 20 transformed mindist^2''S:')
5500 FORMAT (t11, '                        p-value=', f8.6)
5600 FORMAT (/  '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/)
END SUBROUTINE mindist



FUNCTION jtbl() RESULT(fn_val)
INTEGER  :: fn_val

INTEGER, SAVE  :: b(4096)
INTEGER, SAVE  :: j = 4097

IF (j > 4096) THEN
  READ (1,REC=jkk) b(1:4096)
  j = 1
  jkk = jkk + 1
 
END IF
fn_val = b(j)
j = j + 1

RETURN
END FUNCTION jtbl



FUNCTION jkreset() RESULT(fn_val)
INTEGER  :: fn_val

jkk = 1
fn_val = 1

RETURN
END FUNCTION jkreset



SUBROUTINE runtest(filename)
!     This is test 14

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::     This is the RUNS test.  It counts runs up, and runs down, ::
!     :: in a sequence of uniform [0,1) variables, obtained by float-  ::
!     :: ing the 32-bit integers in the specified file. This example   ::
!     :: shows how runs are counted:  .123,.357,.789,.425,.224,.416,.95::
!     :: contains an up-run of length 3, a down-run of length 2 and an ::
!     :: up-run of (at least) 2, depending on the next values.  The    ::
!     :: covariance matrices for the runs-up and runs-down are well    ::
!     :: known, leading to chisquare tests for quadratic forms in the  ::
!     :: weak inverses of the covariance matrices.  Runs are counted   ::
!     :: for sequences of length 10,000.  This is done ten times. Then ::
!     :: repeated.                                                     ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


CHARACTER (LEN=25), INTENT(IN)  :: filename

REAL     :: x(10000), up(100), dn(100)
INTEGER  :: i, ifault, ij, ijkn, ns, nxs
REAL     :: dv, p, uv

CHARACTER (LEN=80)  :: text(13)
CHARACTER (LEN=219) :: dum

OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:13)
WRITE (*, 5000) text(1:13)
WRITE (3,5000) text(1:13)
CLOSE (4)

ns = 10
nxs = 10000
jkk = jkreset()
WRITE (3,5100) filename
WRITE (*, 5100) filename
DO  ijkn = 1, 2
  DO  ij = 1, ns
    DO  i = 1, nxs
      x(i) = jtbl() * 2.328306E-10
    END DO
    CALL udruns(x, nxs, uv, dv, ifault)
    up(ij) = ch(uv)
    dn(ij) = ch(dv)
  END DO
  CALL kstest(up,ns,p)
  WRITE (*, 5200) filename
  WRITE (3,5200) filename
  WRITE (3,5300) p
  WRITE (*, 5300) p
  CALL kstest(dn,ns,p)
  WRITE (3,5400) p
  WRITE (*, 5400) p
END DO
jkk = jkreset()
WRITE (*, 5500)
WRITE (3,5500)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('           The RUNS test for file ', a15/   &
             '     Up and down runs in a sample of 10000'/   &
             '_________________________________________________ ')
5200 FORMAT (t16, '  Run test for ', a15, ':')
5300 FORMAT (t5, '   runs up; ks test for 10 p''S:', f8.6)
5400 FORMAT (t5, ' runs down; ks test for 10 p''S:', f8.6)
5500 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE runtest



SUBROUTINE udruns(x, n, uv, dv, ifault)
!     Algorithm AS 157 Appl. Statist. (1981) vol. 30, No. 1
!     The Runs-up and Runs-down test.

REAL, INTENT(IN)      :: x(10000)
INTEGER, INTENT(IN)   :: n
REAL, INTENT(OUT)     :: uv
REAL, INTENT(OUT)     :: dv
INTEGER, INTENT(OUT)  :: ifault

INTEGER :: i, j, ucount(6), dcount(6), ru, rd
REAL    :: b(6), rn, uni

!     Set up the A and B matrices.
REAL    :: a(6,6) = RESHAPE( (/  &
           4529.4,    0.0,    0.0,    0.0,      0.0, 0.0,  &
           9044.9, 18097.,    0.0,    0.0,      0.0, 0.0,  &
           13568., 27139., 40721.,    0.0,      0.0, 0.0,  &
           18091., 36187., 54281., 72414.,      0.0, 0.0,  &
           22615., 45234., 67852., 90470.,  113262.,    0.0,  &
           27892., 55789., 83685., 111580., 139476., 172860. /), (/ 6,6 /) )

ifault = 0
IF (n >= 4000) THEN
  DO  j = 2, 6
    DO  i = 1, j-1
      a(j,i) = a(i,j)
    END DO
  END DO
  b(1) = 1. / 6.
  b(2) = 5. / 24.
  b(3) = 11. / 120.
  b(4) = 19. / 720.
  b(5) = 29. / 5040.
  b(6) = 1. / 840.
  DO  i = 1, 6
    ucount(i) = 0
    dcount(i) = 0
  END DO
!     The loop that ends at line 300 determines the number of
!     runs-up and runs-down of length i for i = 1(1)5 and the number
!     of runs-up and runs-down of length greater than or equal to 6.
  ru = 1
  rd = 1
  DO  j = 2, n
    IF (x(j)-x(j-1) < 0.0) THEN
      GO TO 50
    ELSE IF (x(j)-x(j-1) == 0.0) THEN
      GO TO 40
    ELSE
      GO TO 60
    END IF

    40 uni = .4
    IF (uni < .5) THEN
      GO TO 50
    ELSE
      GO TO 60
    END IF

    50 ucount(ru) = ucount(ru) + 1
    ru = 1
    IF (rd < 6) rd = rd + 1
    CYCLE
    60 dcount(rd) = dcount(rd) + 1
    rd = 1
    IF (ru < 6) ru = ru + 1
  END DO
  ucount(ru) = ucount(ru) + 1
  dcount(rd) = dcount(rd) + 1
!      print 21,ucount,dcount

!     Calculate the test statistics uv and dv.
  uv = 0.
  dv = 0.
  rn = n
  DO  i = 1, 6
    DO  j = 1, 6
      uv = uv + (ucount(i) - rn*b(i)) * (ucount(j) - rn * b(j)) * a(i,j)
      dv = dv + (dcount(i) - rn*b(i)) * (dcount(j) - rn * b(j)) * a(i,j)
    END DO
  END DO
  uv = uv / rn
  dv = dv / rn
ELSE
  ifault = n
END IF
RETURN

END SUBROUTINE udruns



SUBROUTINE craptest(filename)
!     This is test 15

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     :: This is the CRAPS TEST. It plays 200,000 games of craps, finds::
!     :: the number of wins and the number of throws necessary to end  ::
!     :: each game.  The number of wins should be (very close to) a    ::
!     :: normal with mean 200000p and variance 200000p(1-p), with      ::
!     :: p=244/495.  Throws necessary to complete the game can vary    ::
!     :: from 1 to infinity, but counts for all>21 are lumped with 21. ::
!     :: A chi-square test is made on the no.-of-throws cell counts.   ::
!     :: Each 32-bit integer from the test file provides the value for ::
!     :: the throw of a die, by floating to [0,1), multiplying by 6    ::
!     :: and taking 1 plus the integer part of the result.             ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


CHARACTER (LEN=25), INTENT(IN)  :: filename

REAL     :: av, e(21), ex, pthrows, pwins, sd, sum, t
INTEGER  :: i, iwin, k, lp, m, ng, nt(21), nthrows, nwins

CHARACTER (LEN=80)   :: text(12)
CHARACTER (LEN=232)  :: dum

jkk = jkreset()
OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:12)
WRITE (*, 5000) text(1:12)
WRITE (3,5000) text(1:12)
CLOSE (4)

e(1) = 1. / 3.
sum = e(1)
DO  k = 2, 20
  e(k) = (27.*(27./36.)**(k-2) + 40.*(26./36.)**(k-2) + 55.*(25./36.)**(k-2))  &
         / 648.
  sum = sum + e(k)
END DO
e(21) = 1. - sum
ng = 200000
nwins = 0
nt(1:21) = 0
DO  i = 1, ng
  lp = kthrow()
  nthrows = 1
  IF (lp == 7 .OR. lp == 11) THEN
    iwin = 1
    GO TO 40
  END IF
  IF (lp == 2 .OR. lp == 3 .OR. lp == 12) THEN
    iwin = 0
    GO TO 40
  END IF
  30 k = kthrow()
  nthrows = nthrows + 1
  IF (k == 7) THEN
    iwin = 0
    GO TO 40
  END IF
  IF (k == lp) THEN
    iwin = 1
    GO TO 40
  END IF
  GO TO 30
  40 m = MIN(21,nthrows)
  nt(m) = nt(m) + 1
  nwins = nwins + iwin
END DO
av = 244. * ng / 495.
sd = SQRT(av*251./495.)
t = (nwins-av) / sd
WRITE (*, 5100) filename
WRITE (3,5100) filename
WRITE (*, 5200) nwins, av
WRITE (3,5200) nwins, av
pwins = phi(t)
WRITE (*, 5300) nwins, t, pwins
WRITE (3,5300) nwins, t, pwins
sum = 0.
DO  i = 1, 21
  ex = ng * e(i)
  sum = sum + (nt(i)-ex) ** 2 / ex
END DO
pthrows = chisq(sum,20)
WRITE (3,5400) sum, pthrows
WRITE (*, 5400) sum, pthrows
sum = 0
DO  i = 1, 21
  ex = ng * e(i)
  sum = sum + (nt(i)-ex) ** 2 / ex
  WRITE (3,5500) i, nt(i), ex, (nt(i)-ex) ** 2 / ex, sum
  WRITE (*, 5500) i, nt(i), ex, (nt(i)-ex) ** 2 / ex, sum
END DO
WRITE (3,5600) filename
WRITE (*, 5600) filename
WRITE (3,5700) pwins
WRITE (*, 5700) pwins
WRITE (3,5800) pthrows
WRITE (*, 5800) pthrows
jkk = jkreset()
WRITE (*, 5900)
WRITE (3,5900)
RETURN

5000 FORMAT (a78)
5100 FORMAT (t16, ' Results of craps test for ', a15/   &
             '  No. of wins:  Observed Expected')
5200 FORMAT (t16, '          ',i12,f12.2)
5300 FORMAT (t16, i8, '= No. of wins, z-score=', f6.3, ' pvalue=', f7.5/   &
             '   Analysis of Throws-per-Game:')
5400 FORMAT (' Chisq=', f7.2, ' for 20 degrees of freedom, p= ', f8.5/   &
             t16, 'Throws Observed Expected  Chisq     Sum')
5500 FORMAT (i19, i9, f11.1, f8.3, f9.3)
5600 FORMAT ('            SUMMARY  FOR ', a15)
5700 FORMAT (t16, ' p-value for no. of wins:', f8.6)
5800 FORMAT (t16, ' p-value for throws/game:', f8.6)
5900 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE craptest



SUBROUTINE cdomso(filename)
!     This is test 6

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::                   THE BITSTREAM TEST                          ::
!     :: The file under test is viewed as a stream of bits. Call them  ::
!     :: b1,b2,... .  Consider an alphabet with two "letters", 0 and 1 ::
!     :: and think of the stream of bits as a succession of 20-letter  ::
!     :: "words", overlapping.  Thus the first word is b1b2...b20, the ::
!     :: second is b2b3...b21, and so on.  The bitstream test counts   ::
!     :: the number of missing 20-letter (20-bit) words in a string of ::
!     :: 2^21 overlapping 20-letter words.  There are 2^20 possible 20 ::
!     :: letter words.  For a truly random string of 2^21+19 bits, the ::
!     :: number of missing words j should be (very close to) normally  ::
!     :: distributed with mean 141,909 and sigma 428.  Thus            ::
!     ::  (j-141909)/428 should be a standard normal variate (z score) ::
!     :: that leads to a uniform [0,1) p value.  The test is repeated  ::
!     :: twenty times.                                                 ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!***** NUMBER OF MISSING WORDS IN A STRING OF 2**21 k-LETTER WORDS,****
!***** EACH LETTER 20/k BITS. THERE ARE 2**20 POSSIBLE WORDS**************
!***** EACH OF THE 32 BITS IN THE 2**15 W-TABLE IDENTIFIES A WORD*******
!******** mean should be 141,909 with sigma=290


CHARACTER (LEN=25), INTENT(IN)  :: filename

INTEGER            :: w(0:32767), mbit(0:31)
REAL, SAVE         :: sigs(3) = (/ 290., 295., 339. /)
CHARACTER (LEN= 4) :: ctest(3) = (/ 'OPSO', 'OQSO', ' DNA' /)
CHARACTER (LEN=80) :: text(36), dum
INTEGER            :: i, ic, INDEX, j, jk,  &
                      k, kij, kk, kount, kpow, krk, kr, l, lk, mk, mkk,  &
                      nt, ntries
REAL               :: TRUE, x


OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:36)
WRITE (*, 5000) text(1:36)
WRITE (3,5000) text(1:36)
CLOSE (4)

DO  jk = 1, 3
  k = 2 * jk
  IF (jk == 3) k = 10
  INDEX = (73-(k-9)**2) / 24
  WRITE (*, *) ctest(INDEX), ' test for generator ', filename
  WRITE (*, *)  &
         ' Output: No. missing words (mw), equiv normal variate (z), p-value (p)'
  WRITE (*, *) 'mw     z     p'
  WRITE (3,*) ctest(INDEX), ' test for generator ', filename
  WRITE (3,*)  &
         ' Output: No. missing words (mw), equiv normal variate (z), p-value (p)'
  WRITE (3,*) '   mw     z     p'
  ntries = 1
  kpow = 21
  DO  krk = 33 - 20 / k, 1, -1
    jkk = jkreset()
    kr = 33 - 20 / k - krk
    mk = 2 ** (20/k) - 1
    mkk = 2 ** (20-20/k) - 1
    lk = 20 / k
    DO  kij = 1, ntries
      kpow = 21
!                                  ****SET MASK BITS***************
      !mbit(0) = 1
      DO  i = 0, 31
        mbit(i) = ishft(1,i)
      END DO
!*********** INITIALIZE*****************
      true = 2 ** 20 * EXP(-2.**(kpow-20))
!*****MAIN LOOP*********
      DO  nt = 1, ntries
!                  ********SET W-TABLE TO ZEROS*******
        w(0:32767) = 0

!**** GET INITIAL WORD
        j = ikbit(kr, mk)
        DO  i = 1, k - 1
          j = 2 ** (20/k) * j + ikbit(kr, mk)
        END DO

!****  GENERATE 2**kpow OVERLAPPING WORDS********
        DO  ic = 1, 2 ** kpow
!         *** GET NEW J *****
          j = ISHFT(IAND(j,mkk),lk) + ikbit(kr, mk)
!         *** GET BIT INDEX FROM LAST 5 BITS OF J  ***
          l = IAND(j,31)
!         *** GET TABLE INDEX FROM LEADING 15 BITS OF J***
          kk = ISHFT(j, -5)
!         *** SET BIT L IN W(Kk) ***
          w(kk) = IOR(w(kk), mbit(l))
        END DO
!                    ********** COUNT NUMBER OF EMPTY CELLS *******
        kount = 0
        DO  kk = 0, 32767
          DO  l = 0, 31
            IF (IAND(w(kk),mbit(l)) == 0) kount = kount + 1
          END DO
        END DO
! ****END OF MAIN LOOP****

        x = (kount-true) / sigs(jk)
        WRITE (3,5100) ctest(INDEX), filename, 33 - 20 / k - kr,  &
            32 - kr, kount, x, phi(x)
        WRITE (*, 5100) ctest(INDEX), filename, 33 - 20 / k - kr, 32 - kr, &
                    kount, x, phi(x)
      END DO
      jkk = jkreset()
    END DO
  END DO
END DO
WRITE (*, 5200)
WRITE (3,5200)
RETURN

5000 FORMAT (a78)
5100 FORMAT (a8,' for ',a15,' using bits ',i2,' to ',i2,i14,f7.3,f7.4)
5200 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE cdomso



SUBROUTINE sknt1s(filename)
!     This is test 7

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::     This is the COUNT-THE-1's TEST on a stream of bytes.      ::
!     :: Consider the file under test as a stream of bytes (four per   ::
!     :: 32 bit integer).  Each byte can contain from 0 to 8 1's,      ::
!     :: with probabilities 1,8,28,56,70,56,28,8,1 over 256.  Now let  ::
!     :: the stream of bytes provide a string of overlapping  5-letter ::
!     :: words, each "letter" taking values A,B,C,D,E. The letters are ::
!     :: determined by the number of 1's in a byte::  0,1,or 2 yield A,::
!     :: 3 yields B, 4 yields C, 5 yields D and 6,7 or 8 yield E. Thus ::
!     :: we have a monkey at a typewriter hitting five keys with vari- ::
!     :: ous probabilities (37,56,70,56,37 over 256).  There are 5^5   ::
!     :: possible 5-letter words, and from a string of 256,000 (over-  ::
!     :: lapping) 5-letter words, counts are made on the frequencies   ::
!     :: for each word.   The quadratic form in the weak inverse of    ::
!     :: the covariance matrix of the cell counts provides a chisquare ::
!     :: test::  Q5-Q4, the difference of the naive Pearson sums of    ::
!     :: (OBS-EXP)^2/EXP on counts for 5- and 4-letter cell counts.    ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!     OBC: overlapping-bit-count from stream of bytes


CHARACTER (LEN=25), INTENT(IN)  :: filename

INTEGER  :: w, t(0:3124), s(0:624), kbits(0:255)
REAL     :: chsq, e, q4, q5, z
INTEGER  :: i, i1, i2, ii, j, jj, jk, ks, n

CHARACTER (LEN=80)  :: text(18)
CHARACTER (LEN=113) :: dum
INTEGER, PARAMETER  :: p(0:4) = (/ 37, 56, 70, 56, 37 /)

OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:18)
WRITE (*, 5000) text(1:18)
WRITE (3,5000) text(1:18)
CLOSE (4)

!******Create kbits table: kbits(j)=truncated no of bits in j, -128<=j<=127
!****Filename reads one byte at a time, as integer*1, so -128 to 127*****
DO  jj = 0, 255
  j = jj
  ks = 0
  DO  i = 1, 8
    ks = ks + IAND(j,1)
    j = ISHFT(j, -1)
  END DO
  IF (ks < 2) ks = 2
  IF (ks > 6) ks = 6
  kbits(jj) = ks - 2
END DO

n = 100
jkk = jkreset()
WRITE (3,5100) filename
WRITE (*, 5100) filename
DO  jk = 1, 2
  DO  i = 0, 5 ** 4 - 1
    s(i) = 0
  END DO
  DO  i = 0, 5 ** 5 - 1
    t(i) = 0
  END DO

!***** generate initial word with 5 random keystrokes:
  w = 5 ** 4 * kbits( jtbl8() ) + 5 ** 3 * kbits( jtbl8() ) +  &
      5 ** 2 * kbits( jtbl8() ) + 5 * kbits( jtbl8() ) + kbits( jtbl8() )
  DO  i2 = 1, n
    DO  i1 = 1, 25600
!******Erase leftmost letter of w:
      w = MOD(w,5**4)
!******Boost count for that 4-letter word:
      s(w) = s(w) + 1
!******Shift w left, add new letter, boost 5-letter word count:
      w = 5 * w + kbits( jtbl8() )
      t(w) = t(w) + 1
    END DO
  END DO

!****  Find q4: sum(obs-exp)^2/exp for 4-letter words
  q4 = 0
  DO  ii = 0, 5 ** 4 - 1
    i = ii
    e = 25600 * n
    DO  j = 0, 3
      e = e * p(MOD(i,5)) / 256.
      i = i / 5
    END DO
    q4 = q4 + (s(ii)-e) ** 2 / e
  END DO
!****  Find q5: sum(obs-exp)^2/exp for 5-letter words
  q5 = 0
  DO  ii = 0, 5 ** 5 - 1
    i = ii
    e = 25600 * n
    DO  j = 0, 4
      e = e * p(MOD(i,5)) / 256.
      i = i / 5
    END DO
    q5 = q5 + (t(ii) - e) ** 2 / e
  END DO
  chsq = q5 - q4
  z = (chsq-2500.) / SQRT(5000.)
  IF (jk == 1) THEN
    WRITE (*, 5200) 25600 * n
    WRITE (3,5200) 25600 * n
    WRITE (*, *) '  Results for COUNT-THE-1''S IN SUCCESSIVE BYTES:'
    WRITE (3,*) ' Results fo COUNT-THE-1''S IN SUCCESSIVE BYTES:'
  END IF
  WRITE (3,5300) filename, chsq, z, phi(z)
  WRITE (*, 5300) filename, chsq, z, phi(z)
END DO
jkk = jkreset()
WRITE (*, 5400)
WRITE (3,5400)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('   Test results for ',a15)
5200 FORMAT (' Chi-square with 5^5-5^4=2500 d.of f. for sample size:', i7/  &
     t32, ' chisquare  equiv normal  p-value')
5300 FORMAT (' byte stream for ', a15, f9.2, f11.3, f13.6)
5400 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE sknt1s



SUBROUTINE base5(w)
INTEGER, INTENT(IN)  :: w

INTEGER :: i, j, l(0:4)

j = w
DO  i = 4, 0, -1
  l(i) = MOD(j,5)
  j = j / 5
END DO
WRITE (*, 5000) w, l
RETURN

5000 FORMAT (i12,i3,4I1)
END SUBROUTINE base5



SUBROUTINE wknt1s(filename)
!     This is test 8

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::     This is the COUNT-THE-1's TEST for specific bytes.        ::
!     :: Consider the file under test as a stream of 32-bit integers.  ::
!     :: From each integer, a specific byte is chosen , say the left-  ::
!     :: most::  bits 1 to 8. Each byte can contain from 0 to 8 1's,   ::
!     :: with probabilitie 1,8,28,56,70,56,28,8,1 over 256.  Now let   ::
!     :: the specified bytes from successive integers provide a string ::
!     :: of (overlapping) 5-letter words, each "letter" taking values  ::
!     :: A,B,C,D,E. The letters are determined  by the number of 1's,  ::
!     :: in that byte::  0,1,or 2 ---> A, 3 ---> B, 4 ---> C, 5 ---> D,::
!     :: and  6,7 or 8 ---> E.  Thus we have a monkey at a typewriter  ::
!     :: hitting five keys with with various probabilities::  37,56,70,::
!     :: 56,37 over 256. There are 5^5 possible 5-letter words, and    ::
!     :: from a string of 256,000 (overlapping) 5-letter words, counts ::
!     :: are made on the frequencies for each word. The quadratic form ::
!     :: in the weak inverse of the covariance matrix of the cell      ::
!     :: counts provides a chisquare test::  Q5-Q4, the difference of  ::
!     :: the naive Pearson  sums of (OBS-EXP)^2/EXP on counts for 5-   ::
!     :: and 4-letter cell counts.                                     ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!     OBC: overlapping-bit-count in specified bytes

CHARACTER (LEN=25), INTENT(IN)  :: filename

INTEGER  :: w, t(0:3124), s(0:624)

CHARACTER (LEN=80)  :: text(20)
CHARACTER (LEN=131) :: dum
INTEGER, PARAMETER  :: p(0:4) = (/ 37, 56, 70, 56, 37 /)
INTEGER, PARAMETER  :: k(0:255) = (/  &
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 2, 0, 0, 0,  &
    1, 0, 1, 1, 2, 0, 1, 1, 2, 1, 2, 2, 3, 0, 0, 0, 1, 0, 1, 1, 2,  &
    0, 1, 1, 2, 1, 2, 2, 3, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2,  &
    3, 3, 4, 0, 0, 0, 1, 0, 1, 1, 2, 0, 1, 1, 2, 1, 2, 2, 3, 0, 1,  &
    1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2,  &
    3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4,  &
    3, 4, 4, 4, 0, 0, 0, 1, 0, 1, 1, 2, 0, 1, 1, 2, 1, 2, 2, 3, 0,  &
    1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2,  &
    2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3,  &
    4, 3, 4, 4, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,  &
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 4, 1, 2, 2, 3, 2,  &
    3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 4, 2, 3, 3, 4, 3, 4, 4, 4, 3, 4,  &
    4, 4, 4, 4, 4, 4 /)

INTEGER  :: i, i1, i2, ii, j, jk, l, m, n
REAL     :: chsq, e, q4, q5, z

n = 10
OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:20)
WRITE (*, 5000) text(1:20)
WRITE (3,5000) text(1:20)
CLOSE (4)

WRITE (*, 5100) filename
DO  jk = 1, 25
  jkk = jkreset()
  s(0:624) = 0
  t(0:3124) = 0
  i = k(IAND(ISHFT(jtbl(), -25+jk), 255))
  j = k(IAND(ISHFT(jtbl(), -25+jk), 255))
  l = k(IAND(ISHFT(jtbl(), -25+jk), 255))
  m = k(IAND(ISHFT(jtbl(), -25+jk), 255))
!***** generate initial word with 5 random keystrokes:
  w = 625 * k(IAND(ISHFT(jtbl(), -25+jk), 255)) +  &
      125 * i + 25 * j + 5 * l + m
  DO  i1 = 1, 25600
    DO  i2 = 1, n
!******Erase leftmost letter of w:
      w = MOD(w,5**4)
!******Boost count for that 4-letter word:
      s(w) = s(w) + 1
!******Shift w left, add new letter, boost 5-letter word count:
      i = k(IAND(ISHFT(jtbl(), -25+jk), 255))
      w = 5 * w + i
      t(w) = t(w) + 1
    END DO
  END DO

!****  Find q4: sum(obs-exp)**2/exp for 4-letter words
  q4 = 0
  DO  ii = 0, 5 ** 4 - 1
    i = ii
    e = 25600 * n
    DO  j = 1, 4
      e = e * p(MOD(i,5)) * 2. ** (-8)
      i = i / 5
    END DO
    q4 = q4 + (s(ii)-e) ** 2 / e
  END DO
!****  Find q5: sum(obs-exp)**2/exp for 5-letter words
  q5 = 0
  DO  ii = 0, 5 ** 5 - 1
    i = ii
    e = 25600 * n
    DO  j = 1, 5
      e = e * p(MOD(i,5)) * 2. ** (-8)
      i = i / 5
    END DO
    q5 = q5 + (t(ii)-e) ** 2 / e
  END DO
  chsq = q5 - q4
  z = (chsq-2500.) / SQRT(5000.)
  IF (jk == 1) THEN
    WRITE (3,5200) 25600 * n
    WRITE (*, 5200) 25600 * n
    WRITE (*, *) '  Results for COUNT-THE 1''S IN SPECIFIED BYTES:'
    WRITE (3,*) ' Results for COUNT-THE-1''S IN SPECIFIED BYTES:'
  END IF
  WRITE (3,5300) jk, jk + 7, chsq, z, phi(z)
  WRITE (*, 5300) jk, jk + 7, chsq, z, phi(z)
  jkk = jkreset()
END DO
WRITE (*, 5400)
WRITE (3,5400)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('   Test results for ',a15)
5200 FORMAT (' Chi-square with 5^5-5^4=2500 d.of f. for sample size:', i7/ &
             t10, '             chisquare  equiv normal  p value')
5300 FORMAT ('           bits ', i2, ' to ', i2, f9.2, f11.3, f13.6)
5400 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE wknt1s



SUBROUTINE cdbinrnk(filename)
!     This is test 4

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     :: This is the BINARY RANK TEST for 6x8 matrices.  From each of  ::
!     :: six random 32-bit integers from the generator under test, a   ::
!     :: specified byte is chosen, and the resulting six bytes form a  ::
!     :: 6x8 binary matrix whose rank is determined.  That rank can be ::
!     :: from 0 to 6, but ranks 0,1,2,3 are rare; their counts are     ::
!     :: pooled with those for rank 4. Ranks are found for 100,000     ::
!     :: random matrices, and a chi-square test is performed on        ::
!     :: counts for ranks 6,5 and <=4.                                 ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!*******Test ranks of 100,000 6x8 binary matrices**************
!*******Each row a byte from a RNG, overlapping rows*************


CHARACTER (LEN=25), INTENT(IN)  :: filename

REAL :: pp(0:25)

CHARACTER (LEN=80) :: text(10), dum
CHARACTER (LEN= 6) :: rk(4:6) = (/ ' r<=4 ', ' r =5 ', ' r =6 ' /)
INTEGER            :: r(31), k(2:6)
INTEGER            :: i, ij, ir, kk, kr, l, mr
REAL               :: e, pks, s, t

REAL, PARAMETER  :: p(2:6) =  &
           (/0.149858E-06, 0.808926E-04, 0.936197E-02, 0.217439, 0.773118 /)
!*** rank 2 to 6 with prob p(2),...,p(6); 2,3,4 pooled.

OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:10)
WRITE (*, 5000) text(1:10)
WRITE (3,5000) text(1:10)
CLOSE (4)

!       open(3,file='testout')
WRITE (3,5100) filename
WRITE (*, 5100) filename
DO  ij = 25, 1, -1
  jkk = jkreset()
  WRITE (*, 5200) filename
  WRITE (3, 5200) filename
  kr = ij - 1
  DO  kk = 2, 6
    k(kk) = 0
  END DO
  WRITE (*, 5300) 25 - kr, 32 - kr
  WRITE (3, 5300) 25 - kr, 32 - kr
  WRITE (*, 5400)
  WRITE (3, 5400)
  DO  l = 1, 100000
    DO  i = 1, 6
      r(i) = IAND(ISHFT(jtbl(), -kr), 255)
    END DO
    CALL rankb(r, 6, 8, ir)
    mr = MAX(4, ir)
    k(mr) = k(mr) + 1
  END DO
  s = 0
  DO  l = 4, 6
    IF (l > 4) THEN
      e = 100000 * p(l)
    ELSE
      e = 100000 * (p(2)+p(3)+p(4))
    END IF
    t = (k(l)-e) ** 2 / e
    s = s + t
    WRITE (3, 5500) rk(l), k(l), e, (k(l)-e) ** 2 / e, s
    WRITE (*, 5500) rk(l), k(l), e, (k(l)-e) ** 2 / e, s
  END DO
  pp(kr) = 1. - EXP(-s/2)
  WRITE (*, 5600) pp(kr)
  WRITE (3, 5600) pp(kr)
  jkk = jkreset()
END DO

WRITE (*, 5700)
WRITE (3,5700)
WRITE (*, 5800) (pp(i),i=24,0,-1)
WRITE (3,5800) (pp(i),i = 24,0,-1)
CALL asort(pp,25)
WRITE (*, 5900) filename
WRITE (3,5900) filename
CALL kstest(pp,25,pks)
WRITE (*, 6000) pks
WRITE (3,6000) pks
WRITE (*, 6100)
WRITE (3,6100)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('         Binary Rank Test for ', a15)
5200 FORMAT ('        Rank of a 6x8 binary matrix, '/   &
    '     rows formed from eight bits of the RNG ', a15)
5300 FORMAT ('     b-rank test for bits ', i2, ' to ', i2)
5400 FORMAT (t16, '      OBSERVED   EXPECTED     (O-E)^2/E      SUM')
5500 FORMAT (t7, a9, i12, f12.1, f12.3, f12.3)
5600 FORMAT (t13, '            p=1-exp(-SUM/2)=', f7.5)
5700 FORMAT ( '   TEST SUMMARY, 25 tests on 100,000 random 6x8 matrices'/   &
    ' These should be 25 uniform [0,1] random variables:')
5800 FORMAT (5F12.6)
5900 FORMAT ('   brank test summary for ', a15/   &
    '       The KS test for those 25 supposed UNI''S YIELDS')
6000 FORMAT ('                    KS p-value=', f8.6)
6100 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/)
END SUBROUTINE cdbinrnk



SUBROUTINE rankb(r, m, n, fn_val)

INTEGER, INTENT(IN OUT)  :: r(31)
INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(OUT)     :: fn_val

INTEGER  :: i, ii, j, k, x

INTEGER, SAVE  :: msk(31) = (/  &
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512,  &
    1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144,  &
    524288, 1048576, 2097152, 4194304, 8388608, 16777216,  &
    33554432, 67108864, 134217728, 268435456, 536870912, 0 /)

msk(31) = 32768 * 32768
fn_val = 0
j = n
i = 1
10 DO  ii = i, m
  IF (IAND(r(ii),msk(j)) == msk(j)) THEN
    x = r(ii)
    r(ii) = r(i)
    r(i) = x
    DO  k = i + 1, m
      IF (IAND(r(k),msk(j)) == msk(j)) r(k) = IEOR(r(k),x)
    END DO
    fn_val = fn_val + 1
    IF (i == m .OR. j == 1) RETURN
    j = j - 1
    i = i + 1
    GO TO 10
  END IF
END DO
j = j - 1
IF (j == 0) RETURN
GO TO 10
END SUBROUTINE rankb



SUBROUTINE rank3132(filename)
!     This is test 3

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     :: This is the BINARY RANK TEST for 31x31 matrices. The leftmost ::
!     :: 31 bits of 31 random integers from the test sequence are used ::
!     :: to form a 31x31 binary matrix over the field {0,1}. The rank  ::
!     :: is determined. That rank can be from 0 to 31, but ranks< 28   ::
!     :: are rare, and their counts are pooled with those for rank 28. ::
!     :: Ranks are found for 40,000 such random matrices and a chisqua-::
!     :: re test is performed on counts for ranks 31,30,29 and <=28.   ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     :: This is the BINARY RANK TEST for 32x32 matrices. A random 32x ::
!     :: 32 binary matrix is formed, each row a 32-bit random integer. ::
!     :: The rank is determined. That rank can be from 0 to 32, ranks  ::
!     :: less than 29 are rare, and their counts are pooled with those ::
!     :: for rank 29.  Ranks are found for 40,000 such random matrices ::
!     :: and a chisquare test is performed on counts for ranks  32,31, ::
!     :: 30 and <=29.                                                  ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! see original file \f\bprint.for that displays each step in the
!  rank reduction.
!  finds rank of 31x31 and 32x32 matrices.
! For the 31x31, uses 31 leftmost bits of a 32-bit integer
! to form a row of the binary matrix.
! For the 32x32, uses 32 full integer words for each of 32 rows
!      function mrank(r,m,n)
!   for nxn matrices, to at least 6 places,
!  the probability of rank n-2,n-1,n are all virtually the same.
!    r          p
!  <=29    .0052854502
!    30    .1283502644
!    31    .5775761902
!    32    .2887880952
!c**** Finds binary rank of m rows, n trailing bits each**********


CHARACTER (LEN=25), INTENT(IN)  :: filename

INTEGER            :: i, ij, ir, k, m, n, ntries, ntry, row(32), tbl(0:3)
REAL               :: d, e, s
CHARACTER (LEN=80) :: text(18), dum

REAL, PARAMETER    :: p(0:3) =  &
                      (/ .2887880952, .5775761902, .1283502644, .0052854502 /)

OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:18)
CLOSE (4)
DO  m = 31, 32
  jkk = jkreset()
  IF (m == 31) WRITE (*, 5000) text(1:9)
  IF (m == 31) WRITE (3,5000) text(1:9)
  IF (m == 32) WRITE (*, 5000) text(10:18)
  IF (m == 32) WRITE (3,5000) text(10:18)
  n = m
  WRITE (*, 5100) filename
  WRITE (3,5100) filename
  WRITE (*, 5200) m, n
  WRITE (3,5200) m, n
  WRITE (*, 5300) m
  WRITE (3,5300) m
  DO  i = 0, 3
    tbl(i) = 0
  END DO
  ntries = 40000
  DO  ij = 1, ntries
    DO  i = 1, m
      row(i) = ISHFT(jtbl(), -32+m)
    END DO
    ntry = ntry + 1
    CALL rank(row, m, n, ir)
    k = MIN(n - ir, 3)
    tbl(k) = tbl(k) + 1
  END DO
  s = 0
  WRITE (*, 5400)
  WRITE (3,5400)
  DO  i = 3, 0, -1
    e = p(i) * ntries
    d = (tbl(i)-e) ** 2 / e
    s = s + d
    WRITE (3,5500) n - i, tbl(i), e, d, s
    WRITE (*, 5500) n - i, tbl(i), e, d, s
  END DO
  WRITE (*, 5600) s, chisq(s,3)
  WRITE (3,5600) s, chisq(s,3)
  jkk = jkreset()
END DO

jkk = jkreset()
WRITE (*, 5700)
WRITE (3,5700)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('    Binary rank test for ', a15)
5200 FORMAT (t8, '  Rank test for ', i2, 'x', i2, ' binary matrices:')
5300 FORMAT (t8, ' rows from leftmost ', i2, ' bits of each 32-bit integer')
5400 FORMAT ('      rank   observed  expected (o-e)^2/e  sum')
5500 FORMAT (2I10, f10.1, f10.6, f9.3)
5600 FORMAT ('  chisquare= ', f6.3, ' for 3 d. of f.; p-value= ', f8.6/   &
             '--------------------------------------------------------------')
5700 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/)
END SUBROUTINE rank3132



SUBROUTINE rank(r, m, n, fn_val)

INTEGER, INTENT(IN OUT)  :: r(32)
INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(OUT)     :: fn_val

INTEGER  :: i, ii, j, k, x

INTEGER, SAVE  :: msk(32) = (/  &
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512,  &
    1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144,  &
    524288, 1048576, 2097152, 4194304, 8388608, 16777216,  &
    33554432, 67108864, 134217728, 268435456, 536870912, 0, 0 /)

msk(31) = ishft(1,30)
msk(32) = ishft(1,31)

fn_val = 0
j = n
i = 1
!33        call mprint(r,m,n)
!*****   find row that starts with a 1 in current column (33-j)
10 DO  ii = i, m
  IF (IAND(r(ii), msk(j)) == msk(j)) THEN
    x = r(ii)
    r(ii) = r(i)
    r(i) = x
    DO  k = i + 1, m
      IF (IAND(r(k), msk(j)) == msk(j)) r(k) = IEOR(r(k), x)
    END DO
    fn_val = fn_val + 1
!           print*,' good row',fn_val,i,x
    IF (i == m .OR. j == 1) RETURN
    j = j - 1
    i = i + 1
    GO TO 10
  END IF
END DO
j = j - 1
IF (j == 0) RETURN
GO TO 10
END SUBROUTINE rank



SUBROUTINE cdoperm5(filename)
!     This is test 2

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::            THE OVERLAPPING 5-PERMUTATION TEST                 ::
!     :: This is the OPERM5 test.  It looks at a sequence of one mill- ::
!     :: ion 32-bit random integers.  Each set of five consecutive     ::
!     :: integers can be in one of 120 states, for the 5! possible or- ::
!     :: derings of five numbers.  Thus the 5th, 6th, 7th,...numbers   ::
!     :: each provide a state. As many thousands of state transitions  ::
!     :: are observed,  cumulative counts are made of the number of    ::
!     :: occurences of each state.  Then the quadratic form in the     ::
!     :: weak inverse of the 120x120 covariance matrix yields a test   ::
!     :: equivalent to the likelihood ratio test that the 120 cell     ::
!     :: counts came from the specified (asymptotically) normal dis-   ::
!     :: tribution with the specified 120x120 covariance matrix (with  ::
!     :: rank 99).  This version uses 1,000,000 integers, twice.       ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!*** overlapping 5-permutations.  Uses 120x120 weak inverse *******
!*** of covariance matrix (in 60x60 blocks).
!***  69069 passes, Randu fails, Weyl fails, SR(15,17),SR(13,18) fail.
!***  F(2,1,*) and F(3,1,*) pass

CHARACTER (LEN=25), INTENT(IN)  :: filename

INTEGER            :: r(60,60), s(60,60), t(120)
INTEGER            :: u(1005)
CHARACTER (LEN=80) :: text(15), dum
INTEGER            :: i, ijkm, j, k, n
REAL               :: av, chsq, x, y

!**** divide r and s elements by (200000*n) for proper cov. inverse
!****    the rank is 99=50+49.
OPEN (4, FILE='operm5d.ata', STATUS='OLD')
READ (4,5000) ((r(i,j),j = i,60),i = 1,60)
READ (4,5000) ((s(i,j),j = i,60),i = 1,60)
CLOSE (4)

OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5100) dum
READ (4,5100) text(1:15)
WRITE (*, 5100) text(1:15)
WRITE (3, 5100) text(1:15)
CLOSE (4)

DO  ijkm = 1, 2
  DO  i = 1, 59
    DO  j = i + 1, 60
      r(j,i) = r(i,j)
      s(j,i) = s(i,j)
    END DO
  END DO
!**********************get counts t(1),...,t(120)******************
  jkk = jkreset()
  n = 1000
  t(1:120) = 0
  DO  i = 1001, 1005
    u(i) = jtbl()
  END DO
  DO  i = 1, n
    DO  j = 1, 5
      u(j) = u(1000+j)
    END DO
    DO  j = 1, 1000
      k = kp(u(j)) + 1
      t(k) = t(k) + 1
      u(j+5) = jtbl()
    END DO
  END DO
!*********************evalute quadratic form in weak inverse*******
  chsq = 0.
  av = n * 2000. / 120.
  DO  i = 1, 60
    x = t(i) + t(i+60) - av
    y = t(i) - t(i+60)
    DO  j = 1, 60
      chsq = chsq + x * r(i,j) * (t(j)+t(j+60)-av) + y * s(i,j) *  &
                    (t(j)-t(j+60))
    END DO
  END DO
  chsq = chsq / (2.e08*n)
  WRITE (3,5200) filename
  WRITE (*, 5200) filename
  WRITE (3,5300) chsq, chisq(chsq,99)
  WRITE (*, 5300) chsq, chisq(chsq,99)
END DO
jkk = jkreset()
RETURN

5000 FORMAT (8I10)
5100 FORMAT (a78)
5200 FORMAT ('           OPERM5 test for file ', a15/   &
             '     For a sample of 1,000,000 consecutive 5-tuples,')
5300 FORMAT (' chisquare for 99 degrees of freedom=', f7.3,'; p-value=', f8.6)
END SUBROUTINE cdoperm5



FUNCTION kp(c) RESULT(fn_val)

INTEGER, INTENT(IN)  :: c(5)
INTEGER              :: fn_val

!         real b(5),c(5)
INTEGER  :: b(5), i, j, l, s, t

INTEGER, PARAMETER  :: map(0:59) = (/  &
    39, 38, 37, 36, 41, 40, 54, 55, 56, 57, 58, 59, 49, 48,  &
    52, 53, 50, 51, 42, 43, 44, 45, 46, 47, 33, 32, 31, 30, 35,  &
    34, 12, 13, 14, 15, 16, 17, 29, 28, 24, 25, 27, 26, 21, 20,  &
    19, 18, 23, 22,  2,  3,  5,  4,  1,  0, 10, 11,  9,  8,  6, 7 /)

b(1:5) = c(1:5)
fn_val = 0
DO  i = 5, 2, -1
  t = b(1)
  l = 1
  DO  j = 2, i
    IF (b(j) >= t) THEN
      t = b(j)
      l = j
    END IF
  END DO
  fn_val = i * fn_val + l - 1
  s = b(i)
  b(i) = b(l)
  b(l) = s
END DO
IF (fn_val < 60) fn_val = map(fn_val)

RETURN
END FUNCTION kp



SUBROUTINE cdbday(filename)
!     This is test 1

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::            This is the BIRTHDAY SPACINGS TEST                 ::
!     :: Choose m birthdays in a year of n days.  List the spacings    ::
!     :: between the birthdays.  If j is the number of values that     ::
!     :: occur more than once in that list, then j is asymptotically   ::
!     :: Poisson distributed with mean m^3/(4n).  Experience shows n   ::
!     :: must be quite large, say n>=2^18, for comparing the results   ::
!     :: to the Poisson distribution with that mean.  This test uses   ::
!     :: n=2^24 and m=2^9,  so that the underlying distribution for j  ::
!     :: is taken to be Poisson with lambda=2^27/(2^26)=2.  A sample   ::
!     :: of 500 j's is taken, and a chi-square goodness of fit test    ::
!     :: provides a p value.  The first test uses bits 1-24 (counting  ::
!     :: from the left) from integers in the specified file.           ::
!     ::   Then the file is closed and reopened. Next, bits 2-25 are   ::
!     :: used to provide birthdays, then 3-26 and so on to bits 9-32.  ::
!     :: Each set of bits provides a p-value, and the nine p-values    ::
!     :: provide a sample for a KSTEST.                                ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!     PROGRAM BDAYTST
!     A PROGRAM TO DO  THE BIRTHDAY-SPACINGS TEST ON NBITS OF A UNI

CHARACTER (LEN=25), INTENT(IN)  :: filename

INTEGER  :: b(4096), c(4096), mspace(1000)
REAL     :: alam, pks(64), pp, s
INTEGER  :: i, j, jb, kr, l, lk, m, mask, nbits, nsampl

CHARACTER (LEN=80) :: text(18), dum

OPEN (4, FILE='tests.txt', STATUS='OLD')
!       read(4,766) (d(j),j=1,0)
READ (4,5000) text(1:18)
WRITE (*, 5000) text(1:18)
WRITE (3,5000) text(1:18)
CLOSE (4)

nbits = 24
m = 512
nsampl = 500
alam = REAL(m) ** 3 / 2.0_dp ** (nbits+2)
WRITE (*, 5200) m, nbits, alam
WRITE (3,5200) m, nbits, alam
WRITE (3,5100) filename
WRITE (*, 5100) filename
!      is=lw-nbits
mask = 2 ** (nbits-1) + 2 ** (nbits-1) - 1
DO  kr = 32 - nbits, 0, -1
  s = 0.
  jkk = jkreset()
  DO  j = 1, nsampl
    DO  i = 1, m
      jb = jtbl()
      b(i) = IAND(ISHFT(jb, -kr), mask)
    END DO
    CALL isort(b,m)
    c(1) = b(1)
    DO  i = 2, m
      c(i) = b(i) - b(i-1)
    END DO
    CALL isort(c,m)
    l = 0
    DO  i = 2, m
      lk = 0
      IF (c(i) == c(i-1)) THEN
        lk = lk + 1
        l = l + 1
      END IF
    END DO
    s = s + l
    mspace(j) = l
  END DO
  WRITE (3,5400) nsampl
  WRITE (*, 5400) nsampl
  WRITE (3,5300) filename, 33 - nbits - kr, 32 - kr, s / nsampl
  WRITE (*, 5300) filename, 33 - nbits - kr, 32 - kr, s / nsampl
  CALL chsqts(alam,mspace,nsampl,pp)
  pks(9-kr) = pp
  jkk = jkreset()
END DO

WRITE (*, 5500) pks(1:9)
WRITE (3,5500) pks(1:9)
CALL kstest(pks,9,pp)
WRITE (*, 5600) pp
WRITE (3,5600) pp
WRITE (*, 5700)
WRITE (3,5700)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('           Results for ', a15)
5200 FORMAT (' BIRTHDAY SPACINGS TEST, M=', i4, ' N=2**', i2, ' LAMBDA=', f8.4)
5300 FORMAT (t11, a16, ' using bits ', i2, ' to ', i2, f8.3, f10.6)
5400 FORMAT (t18, '  For a sample of size', i4, ':     mean   ')
5500 FORMAT ('   The 9 p-values were'/ f15.6, 4F10.6/ f15.6, 4F10.6)
5600 FORMAT ('  A KSTEST for the 9 p-values yields ',f8.6)
5700 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/)
END SUBROUTINE cdbday



!  A SUBROUTINE TO DO A CHISQUARE TEST ON N VALUES FROM
!  A DISCRETE DISTRIBUTION.  SET UP FOR POISSON.  CHANGE P'S FOR OTHERS.
!  REQUIRES ARRAY MSPACE(NSAMPL) THAT GIVES NO. OF DUPLICATE SPACINGS
!  IN EACH OF NSAMPL YEARS.

SUBROUTINE chsqts(lambda, mspace, nsampl, pp)

REAL, INTENT(IN)         :: lambda
INTEGER, INTENT(IN OUT)  :: mspace(1000)
INTEGER, INTENT(IN)      :: nsampl
REAL, INTENT(OUT)        :: pp

INTEGER  :: k(0:500)
REAL     :: ex(0:500), obs(0:500), ps(0:500)
INTEGER  :: i, j, l, lb, LT, m, nj
REAL     :: p, s

nj = INT(lambda + 4 * SQRT(lambda))
  ex(0:nj) = 0.0
  k(0:nj) = 0
  ps(0:nj) = 0.0
p = EXP(-lambda)
ps(0) = p * nsampl
k(0) = 0
j = 0
s = p * nsampl
IF (s > 5.) THEN
  j = 1
  ex(0) = s
  s = 0.
END IF

DO  i = 1, nj
  p = lambda * p / i
  ps(i) = ps(i-1) + p * nsampl
  s = s + p * nsampl
  k(i) = j
  IF (ps(i) > nsampl-5) THEN
    ex(j) = s + nsampl - ps(i)
    DO  l = i + 1, nsampl
      k(l) = j
    END DO
    EXIT
  END IF
  IF (s >= 5.) THEN
    ex(j) = s
    j = j + 1
    s = 0.
  END IF
END DO

obs(0:100) = 0.
DO  i = 1, nsampl
  l = k(mspace(i))
  obs(l) = obs(l) + 1
END DO
s = 0.
DO  m = 0, j
  s = s + (obs(m)-ex(m)) ** 2 / ex(m)
END DO
lb = 0
m = k(0)
WRITE (3,*) ' duplicate ', '      number       number '
WRITE (3,*) ' spacings  ', '     observed     expected'
WRITE (*, *) ' duplicate ', '      number       number '
WRITE (*, *) ' spacings  ', '     observed     expected'
DO  i = 1, 100
  IF (k(i) /= m) THEN
    LT = i - 1
    IF (lb /= LT) WRITE (*, 5000) lb, LT, obs(m), ex(m)
    IF (lb /= LT) WRITE (3, 5000) lb, LT, obs(m), ex(m)
    IF (lb == LT) WRITE (*, 5100) lb, obs(m), ex(m)
    IF (lb == LT) WRITE (3, 5100) lb, obs(m), ex(m)
    m = k(i)
    lb = i
    IF (m == j) EXIT
  END IF
END DO

WRITE (*, 5200) lb, obs(m), ex(m)
WRITE (3,5200) lb, obs(m), ex(m)
pp = chisq(s,j)
WRITE (*, 5300) j, s, pp
WRITE (3,5300) j, s, pp
WRITE (3,*) ' :::::::::::::::::::::::::::::::::::::::::'
WRITE (*, *) ' :::::::::::::::::::::::::::::::::::::::::'
!      print 32,pp
!      WRITE (3,32) pp
!  32  format(' UNI = ',f8.5,' (Probability of a better fit.)')
RETURN
5000 FORMAT (' ', i2, ' to ', i2, f13.0, f13.3)
5100 FORMAT (t4, i6, f13.0, f13.3)
5200 FORMAT (' ', i2, ' to INF', f12.0, f13.3)
5300 FORMAT (' Chisquare with ', i2, ' d.o.f. = ', f8.2, ' p-value= ', f8.6)
END SUBROUTINE chsqts



SUBROUTINE cdosum(filename)
!     This is test 13

!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ::             The  OVERLAPPING SUMS test                        ::
!     :: Integers are floated to get a sequence U(1),U(2),... of uni-  ::
!     :: form [0,1) variables.  Then overlapping sums,                 ::
!     ::   S(1)=U(1)+...+U(100), S2=U(2)+...+U(101),... are formed.    ::
!     :: The S's are virtually normal with a certain covariance mat-   ::
!     :: rix.  A linear transformation of the S's converts them to a   ::
!     :: sequence of independent standard normals, which are converted ::
!     :: to uniform variables for a KSTEST. The  p-values from ten     ::
!     :: KSTESTs are given still another KSTEST.                       ::
!     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


CHARACTER (LEN=25), INTENT(IN)  :: filename

REAL     :: x(100), y(100), t(199), u(100), w(10)
INTEGER  :: i, ij, ik, j, m
REAL     :: a, b, h, p, pk, qq, s

CHARACTER (LEN=80)   :: text(11)
CHARACTER (LEN=208)  :: dum
REAL, PARAMETER      :: f(0:100) =  &
 (/ 0.,    .0017, .0132, .0270, .0406, .0538, .0665, .0787,  &
    .0905, .1020, .1133, .1242, .1349, .1454, .1557, .1659, .1760,  &
    .1859, .1957, .2054, .2150, .2246, .2341, .2436, .2530, .2623,  &
    .2716, .2809, .2902, .2995, .3087, .3180, .3273, .3366, .3459,  &
    .3552, .3645, .3739, .3833, .3928, .4023, .4118, .4213, .4309,  &
    .4406, .4504, .4602, .4701, .4800, .4900, .5000, .5100, .5199,  &
    .5299, .5397, .5495, .5593, .5690, .5787, .5882, .5978, .6073,  &
    .6167, .6260, .6354, .6447, .6540, .6632, .6724, .6817, .6910,  &
    .7003, .7096, .7189, .7282, .7375, .7468, .7562, .7657, .7752,  &
    .7848, .7944, .8041, .8140, .8239, .8340, .8442, .8545, .8650,  &
    .8757, .8867, .8980, .9095, .9214, .9337, .9464, .9595, .9731,  &
    .9868, .9983, 1. /)

m = 100
jkk = jkreset()
OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:11)
WRITE (*, 5000) text(1:11)
WRITE (3,5000) text(1:11)
CLOSE (4)
DO  ik = 1, 10
  DO  ij = 1, 100
    s = 0.
    DO  i = 1, 199
      t(i) = jtbl() * .806549E-9
      IF (i <= m) s = s + t(i)
    END DO
    y(1) = s
    DO  j = 2, m
      y(j) = y(j-1) - t(j-1) + t(m+j-1)
    END DO
!***now y(j)=z(j)+...+z(99+j)
!*** They are virtually normal, mean 0, variance 100, but correlated.
!**** Now a matrix transformation of the y's: x=yM, will make the
!*** x's independent normal.
!***The y's covariance matrix T is Toeplitz with diagonals 100,99,...2,1
!***A Cholesky factorization of T: V'V=T provides M=V^(-1).  The
!***covariance of x=yM is M'TM=I.
!*** The columns of M have at most 3 non-zero elements.
    x(1) = y(1) / SQRT(m+0.)
    x(2) = -(m-1.) * y(1) / SQRT(m*(m+m-1.)) + y(2) * SQRT(m/(m+m-1.))
    qq = x(1) ** 2 + x(2) ** 2
    DO  i = 3, m
      a = m + m + 2 - i
      b = 4 * m + 2 - i - i
      x(i) = y(1)/SQRT(a*b) - SQRT((a-1.)/(b+2.))*y(i-1) + SQRT(a/b)*y(i)
      qq = qq + x(i) ** 2
    END DO
!****now the x's are a bunch of iid rnors with mean 0, variance 1.
!***Does sum(x(i)^2) behave as chisquare with m deg. freedom?
!****now convert  x(1),...,x(m) to uni's.
    DO  i = 1, m
      p = phi(x(i))
      h = 100. * p
      j = h
      x(i) = f(j) + (h-j) * (f(j+1)-f(j))
    END DO
!***test to see if the transformed x's are uniforms.
    CALL kstest(x,m,p)
    u(ij) = p
  END DO
!***Now do a KSTEST on the 100 p's from the tests for normality.
  CALL kstest(u,100,pk)
!***And a KSTEST on the 100 p's from the chisquare tests.
!       call kstest(uu,100,pq)
  w(ik) = pk
!       uq(ik)=pq
  WRITE (*,5200) ik, pk
  WRITE (3,5200) ik, pk
END DO
WRITE (3,5100) filename
WRITE (*,5100) filename
!     &'  Q p-value ',f8.6)
CALL kstest(w,10,p)
!      call kstest(uq,10,pq)
WRITE (3,5300) p
WRITE (*,5300) p
jkk = jkreset()
WRITE (*,5400)
WRITE (3,5400)
RETURN

5000 FORMAT (a78)
5100 FORMAT ('   Results of the OSUM test for ',a15)
5200 FORMAT (t16, ' Test no. ',i2,'      p-value ', 2F8.6)
5300 FORMAT (t8, ' KSTEST on the above 10 p-values: ', 2F8.6)
5400 FORMAT (/ '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/ )
END SUBROUTINE cdosum



SUBROUTINE asort(list, n)

REAL, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN)   :: n

INTEGER  :: iu(33), il(33)
INTEGER  :: i, ij, j, k, l, m
REAL     :: t, tt

m = 1
i = 1
j = n
10 IF (i >= j) GO TO 60
20 k = i
ij = (i+j) / 2
t = list(ij)
IF (list(i) > t) THEN
  list(ij) = list(i)
  list(i) = t
  t = list(ij)
END IF
l = j
IF (list(j) >= t) GO TO 40
list(ij) = list(j)
list(j) = t
t = list(ij)
IF (list(i) <= t) GO TO 40
list(ij) = list(i)
list(i) = t
t = list(ij)
GO TO 40

30 list(l) = list(k)
list(k) = tt
40 l = l - 1
IF (list(l) > t) GO TO 40
tt = list(l)
50 k = k + 1
IF (list(k) < t) GO TO 50
IF (k <= l) GO TO 30
IF (l-i > j-k) THEN
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  GO TO 70
END IF
il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 70

60 m = m - 1
IF (m <= 0) RETURN
i = il(m)
j = iu(m)
70 IF (j-i >= 11) GO TO 20
IF (i == 1) GO TO 10
i = i - 1
80 i = i + 1
IF (i == j) GO TO 60
t = list(i+1)
IF (list(i) <= t) GO TO 80
k = i
90 list(k+1) = list(k)
k = k - 1
IF (t < list(k)) GO TO 90
list(k+1) = t
GO TO 80
END SUBROUTINE asort



SUBROUTINE isort(list,n)

INTEGER, INTENT(IN OUT)  :: list(:)
INTEGER, INTENT(IN)      :: n

INTEGER  :: i, ij, j, k, l, m, t, tt
INTEGER  :: iu(33), il(33)

m = 1
i = 1
j = n
10 IF (i >= j) GO TO 60
20 k = i
ij = (i+j) / 2
t = list(ij)

IF (list(i) > t) THEN
  list(ij) = list(i)
  list(i) = t
  t = list(ij)
END IF
l = j
IF (list(j) >= t) GO TO 40
list(ij) = list(j)
list(j) = t
t = list(ij)
IF (list(i) <= t) GO TO 40
list(ij) = list(i)
list(i) = t
t = list(ij)
GO TO 40
30 list(l) = list(k)
list(k) = tt
40 l = l - 1
IF (list(l) > t) GO TO 40
tt = list(l)
50 k = k + 1
IF (list(k) < t) GO TO 50
IF (k <= l) GO TO 30
IF (l-i > j-k) THEN
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  GO TO 70
END IF
il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 70
60 m = m - 1
IF (m <= 0) RETURN
i = il(m)
j = iu(m)
70 IF (j-i >= 11) GO TO 20
IF (i == 1) GO TO 10
i = i - 1
80 i = i + 1
IF (i == j) GO TO 60
t = list(i+1)
IF (list(i) <= t) GO TO 80
k = i
90 list(k+1) = list(k)
k = k - 1
IF (t < list(k)) GO TO 90
list(k+1) = t
GO TO 80
END SUBROUTINE isort



SUBROUTINE kstest(y,n,p)
!      TO TEST WHETHER A SET OF N REAL NUMBERS IS DRAWN
!      FROM A UNIFORM DISTRIBUTION (KOLMOROGOV-SMIRNOV METHOD)
!      THE TEST IS BASED ON THE DISTANCE BETWEEN THE EMPIRICAL
!      AND THEORETICAL DISTRIBUTION FUNCTIONS
!       USAGE: CALL KSTEST(Y,N,P)
!      Y ...   ARRAY OF REAL NUMBERS HYPOTHETICALLY DRAWN
!              FROM A UNIFORM DISTRIBUTION ON (0,1)
!      N ...   NUMBER OF ELEMENTS IN 'Y'
!      P IS THE PROBABILITY ASSOCIATED WITH THE OBSERVED VALUE
!      OF THE ANDERSON-DARLING STATISTIC: N TIMES THE INTEGRAL
!      OF (FN(X)-X)**2/(X*(1-X))

REAL, INTENT(IN OUT)  :: y(:)
INTEGER, INTENT(IN)   :: n
REAL, INTENT(OUT)     :: p

INTEGER :: i, j, m
REAL    :: e, t, z

INTEGER, PARAMETER  :: l(8,10) = RESHAPE(  &
 (/ 40, 46, 37, 34, 27, 24, 20, 20, 88, 59, 43, 37, 29, 27, 20, 22,  &
    92, 63, 48, 41, 30, 30, 25, 24, 82, 59, 42, 37, 26, 28, 26, 22,  &
    62, 48, 33, 30, 23, 23, 22, 18, 49, 34, 22, 20, 16, 17, 17, 12,  &
    17, 17,  7,  8,  4,  7,  5,  1, 40, 18, 19, 14, 16, 13, 10,  9,  &
    59, 20, 10,  4,  1,  1,  0, -1, 41, 43, 36, 112, 15, 95, 32, 58 /), &
    (/ 8, 10 /) )

CALL asort(y, n)
z = -n * n
DO  i = 1, n
  t = y(i) * (1. - y(n+1-i))
  IF (t < 1.e-20) t = 1.e-20
  z = z - (i+i-1) * LOG(t)
END DO
z = z / n
p = 0.
IF (z >= .01) THEN
  IF (z <= 2.) THEN
    p = 2. * EXP(-1.2337/z) * (1.+z/8.-.04958*z*z/(1.325+z)) / SQRT(z)
  ELSE
    IF (z <= 4.) THEN
      p = 1. - .6621361 * EXP(-1.091638*z) - .95059 * EXP(-2.005138*z)
    ELSE
      p = 1. - .4938691 * EXP(-1.050321*z) - .5946335 * EXP(-1.527198*z)
    END IF
  END IF
END IF
m = MIN(n-2,8)
e = 0.
DO  j = 1, 10
  e = e + l(m,j) * sp(p,j) * .0001
END DO
IF (n > 10) e = 10. * e / n

RETURN
END SUBROUTINE kstest



FUNCTION sp(x, i) RESULT(fn_val)
REAL, INTENT(IN)     :: x
INTEGER, INTENT(IN)  :: i
REAL                 :: fn_val

REAL  :: t

fn_val = 0.
SELECT CASE ( i )
  CASE ( 1:7 )
    t = ABS(10.*x - 0.5 - i)
    IF (t > 1.5) RETURN
    IF (t <= .5) THEN
      fn_val = 1.5 - 2. * t * t
    ELSE
      fn_val = 2.25 - t * (3.-t)
    END IF

  CASE (    8)
    IF (x <= .8 .OR. x >= 1.) RETURN
    fn_val = 100. * (x-.9) ** 2 - 1.

  CASE (    9)
    IF (x <= 0. .OR. x >= .05) RETURN
    IF (x <= .01) THEN
      fn_val = -100. * x
    ELSE
      fn_val = 25. * (x-.05)
    END IF

  CASE (   10)
    IF (x <= .98 .OR. x >= 1.) RETURN
    fn_val = .1 - 10. * ABS(x-.99)

END SELECT

RETURN
END FUNCTION sp



FUNCTION phi(x) RESULT(fn_val)

REAL, INTENT(IN)  :: x
REAL              :: fn_val

REAL (dp), PARAMETER  :: v(0:14) = (/  &
   1.253314137315500_dp, .6556795424187985_dp,  &
    .4213692292880545_dp, .3045902987101033_dp, .2366523829135607_dp,  &
    .1928081047153158_dp, .1623776608968675_dp, .1401041834530502_dp,  &
    .1231319632579329_dp, .1097872825783083_dp,  &
    .9902859647173193D-1, .9017567550106468D-1,  &
    .8276628650136917D-1, .764757610162485D-1, .7106958053885211D-1 /)

REAL (dp)  :: a, b, cphi, h, pwr, sum, z
INTEGER    :: i, j

fn_val = .5 + SIGN(.5, x)
IF (ABS(x) > 7.) RETURN
cphi = .5_dp - SIGN(.5, x)
j = ABS(x) + .5_dp
j = MIN(j,14)
z = j
h = ABS(x) - z
a = v(j)
b = z * a - 1.0_dp
pwr = 1.0_dp
sum = a + h * b
DO  i = 2, 24 - j, 2
  a = (a + z*b) / i
  b = (b + z*a) / (i+1)
  pwr = pwr * h ** 2
  sum = sum + pwr * (a + h*b)
END DO
cphi = sum * EXP(-0.5_dp*x*x - 0.918938533204672_dp)
fn_val = 1.0_dp - cphi
IF (x < 0.0_dp) fn_val = cphi

RETURN
END FUNCTION phi



FUNCTION chisq(x, n) RESULT(fn_val)

REAL, INTENT(IN)     :: x
INTEGER, INTENT(IN)  :: n
REAL                 :: fn_val

INTEGER  :: i, l
REAL     :: d, s, t

fn_val = 0.
IF (x <= 0.) RETURN
IF (n > 20) THEN
  t = ((x/n)**.333333 - 1 + (.222222/n)) / SQRT(.222222/n)
  fn_val = phi(MIN(t,8.))
  RETURN
END IF
l = 4 - MOD(n,2)
d = MIN(1,n/3)
DO  i = l, n, 2
  d = d * x / (i-2)
  fn_val = fn_val + d
END DO
IF (l == 3) THEN
  s = SQRT(MIN(.5*x, 50.))
  fn_val = phi(s/.7071068) - EXP(-MIN(.5*x,50.)) * .564189 * fn_val / s
  RETURN
END IF
fn_val = 1. - EXP(-MIN(.5*x, 50.)) * (1. + fn_val)

RETURN
END FUNCTION chisq

END MODULE Diehard_Tests



PROGRAM diehard
USE Diehard_Tests
IMPLICIT NONE
CHARACTER (LEN=25)   :: filename, fileout
CHARACTER (LEN=80)   :: text(36)
CHARACTER (LEN=244)  :: dum
INTEGER              :: which(15)

OPEN (4, FILE='tests.txt', STATUS='OLD')
READ (4,5000) dum
READ (4,5000) text(1:15)
WRITE (*, 5000) text(1:15)
CLOSE (4)

WRITE (*, '(a)', ADVANCE='NO') ' Enter filename (<=25 characters): '
READ (*, 5100) filename
OPEN (1, FILE=filename, FORM='unformatted', ACCESS='direct', RECL=16384)

WRITE (*, *) ' Enter name of output file (<=25 characters): '
READ 5100, fileout
OPEN (3, FILE=fileout)

WRITE (*, 5000) text(1:14)
WRITE (3, 5000) text(1:14)
WRITE (*, *) ' Which tests do you want performed?'
WRITE (*, *) ' For all tests, enter 15 1`s:'
WRITE (*, *) ' 111111111111111'
WRITE (*, *) ' For, say, tests 1,3,7 and 14, enter'
WRITE (*, *) ' 101000100000010'
WRITE (*, *) '     HERE ARE YOUR CHOICES:'
WRITE (*, *) '     1  Birthday Spacings'
WRITE (*, *) '     2  Overlapping Permutations'
WRITE (*, *) '     3  Ranks of 31x31 and 32x32 matrices'
WRITE (*, *) '     4  Ranks of 6x8 Matrices'
WRITE (*, *) '     5  Monkey Tests on 20-bit Words'
WRITE (*, *) '     6  Monkey Tests OPSO,OQSO,DNA'
WRITE (*, *) '     7  Count the 1`s in a Stream of Bytes'
WRITE (*, *) '     8  Count the 1`s in Specific Bytes'
WRITE (*, *) '     9  Parking Lot Test'
WRITE (*, *) '    10  Minimum Distance Test'
WRITE (*, *) '    11  Random Spheres Test'
WRITE (*, *) '    12  The Sqeeze Test'
WRITE (*, *) '    13  Overlapping Sums Test'
WRITE (*, *) '    14  Runs Test'
WRITE (*, *) '    15  The Craps Test'
WRITE (*, *) ' Enter your choices, 1`s yes, 0`s no. using 15 columns:'
WRITE (*, 5200) ' 123456789012345'
WRITE(*, *)
READ (*, 5300) which

IF(which( 1) == 1) CALL cdbday(filename)
IF(which( 2) == 1) CALL cdoperm5(filename)
IF(which( 3) == 1) CALL rank3132(filename)
IF(which( 4) == 1) CALL cdbinrnk(filename)
IF(which( 5) == 1) CALL cdbitst(filename)
IF(which( 6) == 1) CALL cdomso(filename)
IF(which( 7) == 1) CALL sknt1s(filename)
IF(which( 8) == 1) CALL wknt1s(filename)
IF(which( 9) == 1) CALL cdpark(filename)
IF(which(10) == 1) CALL mindist(filename)
IF(which(11) == 1) CALL d3sphere(filename)
IF(which(12) == 1) CALL sqeez(filename)
IF(which(13) == 1) CALL cdosum(filename)
IF(which(14) == 1) CALL runtest(filename)
IF(which(15) == 1) CALL craptest(filename)

WRITE(*, 5400) fileout
STOP

5000 FORMAT (a78)
5100 FORMAT (a25)
5200 FORMAT (a)
5300 FORMAT (15I1)
5400 FORMAT (' Results of DIEHARD battery of tests sent to file ', a15)
END PROGRAM diehard
