!  Compliments of netlib   Sun Jul  6 04:10:36 EDT 1986
!  algorithm 615, collected algorithms from acm.
!  algorithm appeared in acm-trans. math. software, vol.10, no. 2, june, 1984,
!  pp. 202-206.  Authors: Armstrong, R.D., Beck, P.O. & Kung, M.T.

MODULE toms615
IMPLICIT NONE

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

! Replace the COMMON statements in TOMS615
REAL (dp), SAVE    :: acu, big, bsave(210), lu(20,20), pie(210), sad, sadt,  &
                      savtot(210), sig(6000), sigma(300), zlow, zsave(20)
INTEGER, SAVE      :: ibase(20), ik, ik1, incol(20), index(20), inprob(210), &
                      ipar(20), label(210), level
LOGICAL, SAVE      :: direct, intl

CONTAINS


SUBROUTINE kbest (x, y, m, n, iter, ifault, popt, minin, nmax, mmax, bval, &
                  idex, istat, zl)

!***********************************************************************

!  The purpose of this program is to determine the best subset of parameters
!  to fit a linear regression under an least absolute value criterion.
!  This program utilizes the simplex method of linear programming within a
!  branch-and-bound algorithm to solve the best subset problem.

!  The algorithm is based on the publication:
!     Armstrong,R.D. and M.T. Kung "An algorithm to select the best subset for
!     a least absolute value regression problem", Optimization in Statistics,
!     TIMS Studies of the Management Sciences, vol.33, 931-936 (1982).

!   formal parameters

!   x      real array        input:  values of independent variables such that
!          (nmax,mmax)               each row corresponds to an observation

!   y      real array        input:  values of the dependent variables
!            (nmax)

!   m      integer           input:  number of dependent variables

!   n      integer           input:  number of observations

!   iter   integer           output: number of iterations

!   ifault integer           output: failure indicator
!                               = 0  normal termination
!                               = 1  observation matrix does not have full
!                                    row rank (rank m)
!                               = 2  problem size out of range
!                               = 3  no pivot element found implying near
!                                    singular basis

!   popt    real             input:  percentage deviation from
!                                    optimality allowed

!   minin   integer          input:  minimum number of parameters in the
!                                    model.  best subset of size minin t
!                                    m is obtained.

!   nmax    integer          input:  dimension of rows in x (also y)

!   mmax    integer          input:  dimension of columns in x

!   bval    real array       output: array of optimal beta values for each
!                                    subset.  The beta values for the subset
!                                    of size m are stored in positions
!                                    bval(1), bval(2), ..., bval(m), for the
!                                    subset of size m-1 the values are stored
!                                    in positions bval(m+1), bval(m+2), ...,
!                                    bval(2m-1).  In general, the beta values
!                                    for the optimal subset of size k are
!                                    stored in positions bval(l), ...,
!                                    bval(l-k+1) where
!                                    l = (m*(m+1)-k*(k+1))/2 + 1

!  idex     integer array    output: beta index set for the optimal subset.
!        (((mmax+1)*mmax)/2)         This array is a  parallel array for bval;
!                                    i.e., if bval(j) = 2. and idex(j) = 7
!                                    then beta(7) = 2.7 in the associated
!                                    optimal subset.

!  istat    integer array    input:  parameter status array.
!            (mmax)

!                                        1 if beta(j) is required
!                                          in every model
!                     istat(j) =

!                                        0 otherwise

!  zl       real array       output: best objective value for each subse
!            (mmax)
!                                    zl(j) gives the best objective valu
!                                    for the subset with m-j+1 parameter
!***********************************************************************

! implementation notes:

!     1. the routine uses two machine dependent values acu and big.
!        (these are set at the beginning of this subroutine)
!        a) acu is used to test for "zero". acu should be set to approximately
!           100 * the relative machine accuracy of the system in use.
!        b) big is used to initialize the array zl (described above).
!           big should be set to the largest floating point value assignable.

!     2. both single and double precision versions are supplied.  This version
!        is in double precision.

!     3. array dimensions
!        the code is currently dimensioned to solve problems with up to 20
!        parameters and 300 observations. the dimension sizes are
!        determined as follows
!            20 mmax
!           300 nmax
!           210 (nmax+1)*mmax/2
!          6000 nmax*mmax
!            10 maximum value of ik .

!         if the dimension sizes are changed the common and dimension
!         statements will need to be modified in each subroutine.

!********************************************************************

!  subroutines:

!      calbet   - back-solves a sysytem of equations

!      calcpi   - forward-solves a sysytem of equations

!      kbest    - the driver

!      l1norm   - solves the initial regression problem with all
!                 parameters included in the model

!      phase2   - solves the current regression problem with a primal algorithm

!      setup    - determines the form of the current problem by
!                 choosing the parameter to leave based on a penalty

!      update   - updates lu decomposition matrix

!***********************************************************************

!     description of variables:

!     ik:    length of candidate list
!     beta:  estimates to the current subproblem
!     sad:   the minimum total absolute deviation
!     ibase: the index array of columns of the basis

!        n,m,nmax,mmax,x,y, are not changed in the subroutine;
!          sad is updated at each iteration

!     lu:    the lu decomposition of the current basis
!     index: the index array of rows of lu
!     tot:   the current rhs of the dual problem
!     sigma: indicator array to specify whether a nonbasic dual variable
!            is at upper or lower bound(+1 implies upper; -1 implies lower)
!     inext: a local array for the sort routine
!     rhs:   a local array for the calbet routine

!     ipar(i) = -k if the k-th parameter (beta) is out of model at level i
!             =  k if the k-th parameter is in

!     label(i)  saves the indices of the basic observations in the
!               predecessor path
!     zsave(i)  saves the objective values on the predecessor path

!**********************************************************************

! Conversion to be ELF90 compatible by Alan Miller
! e-mail:  amiller @ bigpond.net.au
! http://users.bigpond.net.au

! Latest revision - 8 December 1998

IMPLICIT NONE

REAL (dp), INTENT(IN)   :: x(300,20), y(300), popt
INTEGER, INTENT(IN)     :: m, n, nmax, mmax
INTEGER, INTENT(IN OUT) :: minin, istat(20)
INTEGER, INTENT(OUT)    :: iter, ifault, idex(210)
REAL (dp), INTENT(OUT)  :: bval(210), zl(20)

! Local variables
REAL (dp) :: beta(20), tot(20), pi(20), repp(20)
INTEGER   :: i, j, jinm, k2, k3, kkk1, l, nfree, numin
REAL (dp) :: popt1, zero = 0.0_dp

!     assign values to machine dependent constants

! acu = SQRT(d1mach(4))
! big = d1mach(2)

acu = SQRT( 2.0_dp * EPSILON(1.0_dp) )
big = HUGE(1.0_dp)

ifault = 0
ik = 5
ik1 = ik-1
popt1 = (100. - popt)/100.

!     test problem size

IF (n <= nmax .AND. m <= mmax) GO TO 10
ifault = 2
RETURN

10 DO i = 1,m
  index(i) = i
  tot(i) = zero
  incol(i) = i
END DO
CALL l1norm (x, y, n, m, ifault, iter, beta, pi, tot)
IF (ifault > 0) RETURN

!     initialization

!     save the initial solution

jinm = 0
DO i = 1,m
  bsave(i) = beta(i)
  bval(i) = beta(i)
  pie(i) = pi(i)
  label(i) = ibase(i)
  idex(i) = i
  inprob(i) = i
  savtot(i) = tot(i)
  zl(i) = big
  jinm = jinm+istat(i)
END DO
sig(1:n) = sigma(1:n)

level = 0
direct = .true.

IF (jinm > minin) minin = jinm

!     initialize number of free parameters in model

nfree = m-jinm

zl(m) = sad
zsave(m) = sad

IF (minin == m) RETURN

!       numin gives the number of parameters in the model

numin = m

!     start the main loop

50 numin = numin-1
level = level+1
IF (numin < minin .OR. nfree == 0) GO TO 120
k2 = (m*(m+1) - (numin+1)*(numin+2))/2
kkk1 = numin+1
zlow = zl(numin)*popt1
IF (direct) GO TO 80
DO i = 1,kkk1
  j = k2+i
  index(i) = i
  ibase(i) = label(j)
  incol(i) = inprob(j)
  pi(i) = pie(j)
  beta(i) = bsave(j)
  tot(i) = savtot(j)
END DO

!     load in the bounds from the immediate predecessor

k3 = (m-kkk1)*n
DO i = 1,n
  k3 = k3+1
  sigma(i) = sig(k3)
END DO
sad = zsave(kkk1)

!     set up a new problem

80 CALL setup (x, kkk1, ifault, istat, beta, tot, pi, repp)
IF (ifault > 0) RETURN
nfree = nfree - 1

!     call primal l.p. code

CALL phase2 (x, y, numin, n, iter, ifault, beta, tot, pi, repp)

!      save the solution data for later recall

k2 = k2 + kkk1 + 1
j = k2
DO i = 1,numin
  l = k2+i
  pie(l) = pi(i)
  label(l) = ibase(i)
  bsave(l) = beta(i)
  inprob(l) = incol(i)
  savtot(l) = tot(i)
END DO

!     save the objective value

zsave(numin) = sad

!     save the nonbasic bounds

k3 = (m-numin)*n
DO i = 1,n
  k3 = k3+1
  sig(k3) = sigma(i)
END DO

direct = .true.

!     check for a better solution

IF (sadt >=  zlow) GO TO 50

!     save the better solution

zl(numin) = sad
k2 = j
DO i = 1,numin
  l = k2 + i
  idex(l) = incol(i)
  bval(l) = beta(i)
END DO

!     go to top of loop

GO TO 50
120 numin = numin + 1

direct = .false.

!     must work back up the tree

level = level - 1
130 IF (ipar(level) > 0) GO TO 140
ipar(level) = -ipar(level)
j = ipar(level)
istat(j) = 1
numin = numin + 1
GO TO 50
140 j = ipar(level)
istat(j) = 0
nfree = nfree + 1
level = level - 1
IF (level > 0) GO TO 130
RETURN

END SUBROUTINE kbest


SUBROUTINE setup (x, m, ifault, istat, beta, tot, pi, repp)

!     this subroutine determines the form of the current subproblem

IMPLICIT NONE
REAL (dp), INTENT(IN)     :: x(300,20), beta(20)
REAL (dp), INTENT(IN OUT) :: tot(20), pi(20), repp(20)
INTEGER, INTENT(OUT)      :: ifault
INTEGER, INTENT(IN OUT)   :: m, istat(20)

! Local variables
INTEGER   :: i, ifixi, ii, iout, isave, kkk, l, leave, ll
REAL (dp) :: beside, one = 1.0_dp, penalt, ratio, rho, rhs(20), side, srepp,  &
             test, zero = 0.0_dp

!     if the immediate predecessor is in memory go on

IF (direct) GO TO 10

!     reconstruct the lu decomposition

kkk = 1
CALL update (kkk, ifault, x, 1, m)
IF (ifault > 0) RETURN

!     recalculate the basic pi values

kkk = 0

CALL calcpi (kkk, pi, tot, m)

!     determine parameter to leave based on penalty

10 penalt = big
kkk = 0
rhs(1:m) = zero

DO i = 1,m
  ii = index(i)
  ll = incol(ii)
  IF (istat(ll) == 1) CYCLE
  rho = SIGN(1._dp, beta(ii))
  rhs(ii) = one
  kkk = i-1
  CALL calcpi (kkk, repp, rhs, m)
  rhs(ii) = zero
  test = big
  DO l = 1,m
    IF (ABS(repp(l)) <= acu) CYCLE
    srepp = SIGN(1._dp,repp(l))
    ratio = ABS((1._dp - rho*srepp*pi(l))/repp(l))
    IF (ratio >= test) CYCLE
    side = rho*srepp
    test = ratio
    iout = l
  END DO
  
!     see if the penalty is less
  
  IF (test*ABS(beta(ii)) >= penalt) CYCLE
  leave = iout
  ifixi = ii
  beside = side
  penalt = test*ABS(beta(ii))
END DO

!     update the objective value

sad = sad+penalt

!     parameter incol(ifixi) will leave the model
!     pi(ibase(leave)) will leave the basis and go to bound

!     switch unwanted parameter and pi to the end of list

isave = incol(ifixi)
incol(ifixi) = incol(m)
incol(m) = isave
tot(ifixi) = tot(m)

!     label this node and flag the outgoing parameter

ipar(level) = -isave
istat(isave) = -1

isave = ibase(leave)
ibase(leave) = ibase(m)
ibase(m) = isave

!     place the outgoing pi at the correct bound

sigma(isave) = beside

!     form a new basis

m = m-1

DO i = 1,m
  index(i) = i
END DO

kkk = 1
CALL update (kkk, ifault, x, 1, m)
IF (ifault > 0) RETURN

!     update the tot array

DO i = 1,m
  ll = incol(i)
  tot(i) = tot(i) - beside*x(isave,ll)
END DO

RETURN

END SUBROUTINE setup


SUBROUTINE phase2 (x, y, m, n, iter, ifault, beta, tot, pi, repp)

!     solves lp with primal algorithm

IMPLICIT NONE

INTEGER, INTENT(IN)       :: m, n
INTEGER, INTENT(IN OUT)   :: iter
INTEGER, INTENT(OUT)      :: ifault
REAL (dp), INTENT(IN)     :: x(300,20), y(300)
REAL (dp), INTENT(IN OUT) :: tot(20)
REAL (dp), INTENT(OUT)    :: beta(20), pi(20), repp(20)

! Local variables
INTEGER   :: i, iin, iout, j, k, kan(10), kkk, krow, l, ll
REAL (dp) :: one = 1.0_dp, ratio, rho, rhs(20), srepp, test, xtwo,  &
             z(10), zero = 0.0_dp, zc, zz

!     z(i) stores the reduced cost values on the list
!     kan(i) stores the row number on the list
!     ik is the length of the list

!     calculate the simplex multipliers

DO i = 1,m
  k = ibase(i)
  pi(i) = y(k)
END DO
kkk = 0
CALL calbet (kkk, beta, pi, m)

!     recalculate the basic pi values

kkk = 0
CALL calcpi (kkk, pi, tot, m)

!     store the objective value used for termination check

sadt = sad
IF (sad >= zlow) RETURN

!     generate the candidate list

20 z(1:ik) = acu

!     calculate the reduced costs

DO j = 1,n
  IF (sigma(j) == zero) CYCLE
  zc = zero
  DO i = 1,m
    ll = incol(i)
    zc = zc+beta(i)*x(j,ll)
  END DO
  zz = (zc - y(j))*sigma(j)
  IF (zz <= z(ik)) CYCLE
  
!     rank the values in z(i), in descending order
  
  DO k = 1,ik1
    IF (zz <= z(k)) CYCLE
    i = ik
    50 z(i) = z(i-1)
    kan(i) = kan(i-1)
    i = i-1
    IF (i > k) GO TO 50
    z(k) = zz
    kan(k) = j
    CYCLE
  END DO
  z(ik) = zz
  kan(ik) = j
END DO
IF (z(1) <= acu) RETURN
krow = 1
zz = z(1)
iin = kan(1)
zc = zz*sigma(iin)
GO TO 100
80 krow = krow+1
IF (krow > ik) GO TO 20
IF (z(krow) <= acu) GO TO 20

!     determine the entering variable from the candidate list

iin = kan(krow)
zc = zero
DO i = 1,m
  ll = incol(i)
  zc = zc + beta(i)*x(iin,ll)
END DO
zc = zc - y(iin)
zz = zc*sigma(iin)
IF (zz <= acu) GO TO 80
100 rho = sigma(iin)

!     rho: sign of the incoming reduced cost
!     iin: incoming candidate

!     find the representation of the entering variable

DO i = 1,m
  ll = incol(i)
  rhs(i) = x(iin,ll)
END DO
kkk = 0
CALL calcpi (kkk, repp, rhs, m)

!     calculate the min ratio test to determine the leaving variable

test = 2.0
DO i = 1,m
  IF (ABS(repp(i)) <= acu) CYCLE
  srepp = SIGN(1.0_dp,repp(i))
  ratio = ABS((1._dp - rho*srepp*pi(i))/repp(i))
  IF (ratio >= test) CYCLE
  test = ratio
  iout = i
END DO

!     perform fathoming test before pivot

sadt = sad + test*zz
IF (sadt >= zlow) RETURN
sad = sadt

pi(1:m) = pi(1:m) + test*repp(1:m)*rho
IF (test < 2.0) GO TO 150

sigma(iin) = -sigma(iin)
xtwo = sigma(iin) + sigma(iin)
DO i = 1,m
  ll = incol(i)
  tot(i) = tot(i) - xtwo*x(iin,ll)
END DO
GO TO 80
150 k = ibase(iout)
sigma(k) = SIGN(1.0_dp, pi(iout))
pi(iout) = sigma(iin)*(1.0-test)
ibase(iout) = iin

DO i = 1,m
  ll = incol(i)
  tot(i) = tot(i) + sigma(iin)*x(iin,ll) - sigma(k)*x(k,ll)
END DO

kkk = iout
CALL update (iout, ifault, x, 1, m)

rhs(1:m) = zero
rhs(kkk) = one
kkk = kkk - 1
CALL calbet (kkk, repp, rhs, m)
DO l = 1,m
  beta(l) = beta(l) - zc*repp(l)
END DO
sigma(iin) = zero
iter = iter + 1
GO TO 80

END SUBROUTINE phase2


SUBROUTINE l1norm (x, y, n, m, ifault, iter, beta, pi, tot)

!     this routine solves the initial regression problem with all
!     the parameters included in the model

IMPLICIT NONE
INTEGER, INTENT(IN)       :: n, m
INTEGER, INTENT(OUT)      :: ifault, iter
REAL (dp), INTENT(IN)     :: x(300,20), y(300)
REAL (dp), INTENT(IN OUT) :: tot(20)
REAL (dp), INTENT(OUT)    :: beta(20), pi(20)

! Local variables
REAL (dp) :: d(300), delta(300), price(300), rhs(20)
INTEGER   :: i, iin, inext(300), iout, ipt, j, k, k1, kkk, kount, l
REAL (dp) :: ahalf, aone, one = 1.0_dp, ratio, rho, subt, t, test, total,  &
             val, zero = 0.0_dp, zzz

!        initial settings

iter = 0
aone = one + acu
ahalf = .5 + acu

!     set up initial lu decomposition

intl = .true.
k = 1
CALL update (k, ifault, x, n, m)
IF (ifault /= 0) RETURN
intl = .false.
DO i = 1,m
  k1 = ibase(i)
  rhs(i) = y(k1)
END DO
kkk = 0
CALL calbet (kkk, beta, rhs, m)

!        calculate initial d, tot and sigma vectors

DO j = 1,n
  val = DOT_PRODUCT( beta(1:m), x(j,1:m) )
  d(j) = y(j) - val
  sigma(j) = SIGN(1._dp, d(j))
END DO
DO j = 1,m
  rhs(j) = zero
  kkk = ibase(j)
  sigma(kkk) = zero
END DO
DO j = 1,n
  DO i = 1,m
    tot(i) = tot(i) - sigma(j)*x(j,i)
  END DO
END DO

!        main iterative loop begins

70 kkk = 0
CALL calcpi (kkk, pi, tot, m)

t = aone
k = 0
DO j = 1,m
  IF (ABS(pi(j)) < t) CYCLE
  k = j
  t = ABS(pi(j))
  rho = -SIGN(1._dp, pi(j))
END DO
IF (k == 0) GO TO 220
kkk = k - 1
rhs(k) = one
CALL calbet (kkk, beta, rhs, m)
rhs(k) = zero
DO i = 1,n
  delta(i) = zero
  IF (sigma(i) == zero) CYCLE
  DO j = 1,m
    delta(i) = delta(i) + beta(j)*x(i,j)
  END DO
  delta(i) = rho*delta(i)
END DO

!        perform partial sort of ratios

t = t*.5
kount = 0
ratio = big
total = ahalf
subt = zero
DO i = 1,n
  IF (delta(i)*sigma(i) <= acu) CYCLE
  test = d(i)/delta(i)
  IF (test >= ratio) CYCLE
  total = total + ABS(delta(i))
  IF (total - subt >= t) GO TO 110
  
!     insert i in list
  
  kount = kount+1
  price(kount) = test
  inext(kount) = i
  CYCLE
  
!     update sum and kick iin out of the list
  
  110 total = total - subt
  ratio = test
  ipt = 0
  kkk = 0
  
!     identify a new iin
  
  120 kkk = kkk + 1
  IF (kkk > kount) GO TO 130
  IF (price(kkk) <= ratio) GO TO 120
  ratio = price(kkk)
  ipt = kkk
  GO TO 120
  130 IF (ipt == 0) GO TO 150
  
!     switch values
  
  kkk = inext(ipt)
  subt = ABS(delta(kkk))
  IF (total - subt < t) GO TO 140
  price(ipt) = price(kount)
  inext(ipt) = inext(kount)
  kount = kount - 1
  GO TO 110
  140 iin = inext(ipt)
  inext(ipt) = i
  price(ipt) = test
  CYCLE
  150 iin = i
  subt = ABS(delta(i))
END DO

!        update basic indicators

DO j = 1,kount
  kkk = inext(j)
  zzz = sigma(kkk)
  DO l = 1,m
    subt = zzz*x(kkk,l)
    tot(l) = tot(l) + subt + subt
  END DO
  sigma(kkk) = -sigma(kkk)
END DO

iout = ibase(k)
delta(iout) = rho
DO l = 1,m
  tot(l) = tot(l) + rho*x(iout,l) + sigma(iin)*x(iin,l)
END DO
sigma(iout) = -rho
ibase(k) = iin
CALL update (k, ifault, x, 1, m)
d(1:n) = d(1:n) - ratio*delta(1:n)
sigma(iin) = zero
iter = iter + 1
GO TO 70

!        calculate optimal beta and sum of absolute deviations

220 DO i = 1,m
  k1 = ibase(i)
  rhs(i) = y(k1)
END DO
kkk = 0
CALL calbet (kkk, beta, rhs, m)
sad = SUM( ABS(d(1:n)) )
RETURN

END SUBROUTINE l1norm


SUBROUTINE calbet (kkk, beta, rhs, m)

!     subroutine calbet for back-solving system of equations

IMPLICIT NONE

INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN OUT) :: kkk
REAL (dp), INTENT(IN)   :: rhs(20)
REAL (dp), INTENT(OUT)  :: beta(20)

! Local variables
INTEGER   :: i, ii, ii1, k, k1, k2, kk, kkk1, m1
REAL (dp) :: zero = 0.0_dp

IF (m > 1) GO TO 10
k = index(1)
beta(k) = rhs(1)/lu(k,1)
RETURN

10 m1 = m-1
IF (kkk == 0) GO TO 30
DO i = 1,kkk
  k = index(i)
  beta(k) = zero
END DO
30 kkk = kkk + 1
k = index(kkk)
beta(k) = rhs(kkk)/lu(k,kkk)
IF (kkk == m) GO TO 60
kkk1 = kkk + 1
DO ii = kkk1,m
  k = index(ii)
  beta(k) = rhs(ii)
  ii1 = ii - 1
  DO i = kkk,ii1
    kk = index(i)
    beta(k) = beta(k) - lu(kk,ii)*beta(kk)
  END DO
  beta(k) = beta(k)/lu(k,ii)
END DO

60 DO ii = 1,m1
  k1 = m - ii
  k = index(k1)
  DO i = 1,ii
    kk = m - i + 1
    k2 = index(kk)
    beta(k) = beta(k) - lu(k2,k1)*beta(k2)
  END DO
END DO
RETURN

END SUBROUTINE calbet


SUBROUTINE update (kkk, ifault, x, n, m)

!     subroutine update for updating lu decomposition matrix

IMPLICIT NONE
INTEGER, INTENT(IN)      :: n, m
INTEGER, INTENT(IN OUT)  :: kkk
INTEGER, INTENT(OUT)     :: ifault
REAL (dp), INTENT(IN)    :: x(300,20)

INTEGER   :: i, icol, ii, ii1, irow, isave, j, k, kk, ll
REAL (dp) :: pivot, subt

irow = 0
DO ii = kkk,m
  IF (intl) GO TO 10
  irow = ibase(ii)
  GO TO 20
  10 irow = irow + 1
  ibase(ii) = irow
  IF (irow > n) ifault = 1
  IF (ifault /= 0) RETURN
  20 DO i = 1,m
    ll = incol(i)
    lu(i,ii) = x(irow,ll)
  END DO
  
!     set up representation of incoming row
  
  IF (ii == 1) GO TO 60
  ii1 = ii - 1
  DO icol = 1,ii1
    k = index(icol)
    subt = lu(k,ii)
    j = icol + 1
    DO i = j,m
      k = index(i)
      lu(k,ii) = lu(k,ii) - subt*lu(k,icol)
    END DO
  END DO
  
!     find maximum entry
  
  60 pivot = acu
  kk = 0
  DO i = ii,m
    k = index(i)
    IF (ABS(lu(k,ii)) <= pivot) CYCLE
    pivot = ABS(lu(k,ii))
    kk = i
  END DO
  IF (kk == 0) GO TO 10
  
!     switch order
  
  isave = index(kk)
  index(kk) = index(ii)
  index(ii) = isave
  
!     put in columns of lu one at a time
  
  IF (ii == m) CYCLE
  j = ii + 1
  DO i = j,m
    k = index(i)
    lu(k,ii) = lu(k,ii)/lu(isave,ii)
  END DO
END DO
kkk = irow
RETURN

END SUBROUTINE update


SUBROUTINE calcpi (kkk, pi, tot, m)

!     this routine forward solves a linear system

IMPLICIT NONE
INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN OUT) :: kkk
REAL (dp), INTENT(IN)   :: tot(20)
REAL (dp), INTENT(OUT)  :: pi(20)

! Local variables
INTEGER   :: i, ii, ii1, k, k1, k2, kkk1, m1
REAL (dp) :: zero = 0.0_dp

m1 = m-1

IF (m > 1) GO TO 10
k = index(m)
pi(m) = tot(k)/lu(k,m)
RETURN

10 pi(1:kkk) = zero

kkk = kkk + 1
k = index(kkk)
pi(kkk) = tot(k)
IF (kkk == m) GO TO 60
kkk1 = kkk + 1
DO ii = kkk1, m
  k = index(ii)
  pi(ii) = tot(k)
  ii1 = ii - 1
  pi(ii) = pi(ii) - DOT_PRODUCT( lu(k, kkk:ii1), pi(kkk:ii1) )
END DO

60 k = index(m)
pi(m) = pi(m)/lu(k,m)
DO ii = 1,m1
  k1 = m - ii
  k = index(k1)
  DO i = 1,ii
    k2 = m - i + 1
    pi(k1) = pi(k1) - lu(k,k2)*pi(k2)
  END DO
  pi(k1) = pi(k1)/lu(k,k1)
END DO

RETURN
END SUBROUTINE calcpi


END MODULE toms615
