MODULE global_minimum

! This is translated from the code of:
! Dr. Tibor Csendes
! Dept. of Applied Informatics, Jozsef Attila University
! H-6701 Szeged, Pf. 652, Hungary

! Phone: +36 62 544 305

! Fax: +36 62 420 292

! E-mail: csendes@inf.u-szeged.hu

! URL: http://www.inf.u-szeged.hu/~csendes/

IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS

!   ROUTINE NAME - GLOBAL
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-27  Time: 13:16:39

!-----------------------------------------------------------------------

!   LATEST REVISION - OKTOBER 23, 1986
!   Latest revision of Fortran 90 version - 30 January 2000

!   PURPOSE  - GLOBAL MINIMUM OF FUNCTION OF N VARIABLES
!              USING A LOCAL SEARCH METHOD

!   USAGE  - CALL GLOBAL (AMIN, AMAX, NPARM, M, N100, NG0, IPR,
!                         NSIG, X0, NC, F0)

!   ARGUMENTS
!  AMIN - VECTOR OF LENGTH NPARM CONTAINING THE LOWER BOUNDS OF THE
!         PARAMETERS, SO X(I) IS SEARCHED IN THE INTERVAL (AMIN(I), AMAX(I)).
!         (INPUT)
!  AMAX - VECTOR OF LENGTH NPARM CONTAINING THE UPPER BOUNDS OF THE
!         PARAMETERS. (INPUT)
! NPARM - NUMBER OF PARAMETERS, 1 <= NPARM <= 20. (INPUT)
!     M - NUMBER OF RESIDUAL FUNCTIONS, WHEN THE OBJECTIVE FUNCTION IS OF THE
!         FORM F1**2 + F2**2 +...+ FM**2, 1 <=M <= 100. (INPUT)
!         N.B. M is NOT used for this purpose!
!              It is passed to the user's routine FUNCT as a parameter.
!  N100 - NUMBER OF SAMPLE POINTS TO BE DRAWN UNIFORMLY IN ONE CYCLE,
!         20 <= N100 <= 10000. THE SUGGESTED VALUE IS 100*NPARM. (INPUT)
!   NG0 - NUMBER OF BEST POINTS SELECTED FROM THE ACTUAL SAMPLE, 1 <= NG0 <= 20.
!         THE SUGGESTED VALUE IS TWICE THE EXPECTED NUMBER OF LOCAL MINIMA.
!         (INPUT)
!   IPR - FORTRAN DATA SET REFERENCE NUMBER WHERE THE PRINTED OUTPUT BE SENT.
!         (INPUT)
!  NSIG - CONVERGENCE CRITERION, THE ACCURACY REQUIRED IN THE PARAMETER
!         ESTIMATES.  THIS CONVERGENCE CRITERION IS SATISFIED IF ON TWO
!         SUCCESSIVE ITERATIONS THE PARAMETER ESTIMATES AGREE, COMPONENT BY
!         COMPONENT, TO NSIG DIGITS.  THE SUGGESTED VALUE IS 6. (INPUT)
!    X0 - OUTPUT NPARM BY 20 MATRIX CONTAINING NC (UP TO 20) LOCAL MINIMIZERS
!         FOUND.
!    NC - NUMBER OF DIFFERENT LOCAL MINIMIZERS FOUND. (OUTPUT)
!    F0 - OUTPUT VECTOR OF NC (UP TO 20) OBJECTIVE FUNCTION VALUES, F0(I)
!         BELONGS TO THE PARAMETERS X0(1,I), X0(2,I),..., X0(NPARM,I).

!   REQUIRED ROUTINES - URDMN, FUN, LOCAL

!-----------------------------------------------------------------------

SUBROUTINE global(amin, amax, nparm, m, n100, ng0, ipr, nsig, x0, nc, f0)

REAL (dp), INTENT(IN)    :: amin(:)
REAL (dp), INTENT(IN)    :: amax(:)
INTEGER, INTENT(IN)      :: nparm
INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN OUT)  :: n100
INTEGER, INTENT(IN OUT)  :: ng0
INTEGER, INTENT(IN)      :: ipr
INTEGER, INTENT(IN)      :: nsig
REAL (dp), INTENT(OUT)   :: x0(:,:)
INTEGER, INTENT(OUT)     :: nc
REAL (dp), INTENT(OUT)   :: f0(:)

REAL (dp) :: x(nparm,100), x1(nparm,20), xcl(nparm,100), r(nparm), w(nparm)
INTEGER   :: ic(100), ic1(20)
INTEGER   :: i, i1, icc, icj, ig, ii, iii, im, in1, inum, inum1, inum2, it, &
             iv, j, jj, l1, maxfn, n, n0, n1, ncp, nfe, nfe1, ng, ng10, nm, &
             nn100, ns
REAL (dp) :: f(100), f1(20), fcl(100), y(nparm), mmin(nparm), mmax(nparm)
REAL (dp) :: a, alfa, b, b1, bb, fc, ff, fm, relcon
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, ten = 10._dp

DO  i=1,nparm
  mmin(i) = amin(i)
  mmax(i) = amax(i)
  IF (mmin(i) == mmax(i)) GO TO 460
END DO
b1 = one/nparm

IF (ng0 < 1) ng0 = 1
IF (ng0 > 20) ng0 = 20
IF (n100 < 20) n100 = 20
IF (n100 > 10000) n100 = 10000
IF (n100 >= 100) GO TO 10

nn100 = n100
n = 1
GO TO 15

10 nn100 = 100
n = n100/100
n100 = n*100

15 ng10 = 100
DO  i=1,ng10
  f(i) = HUGE(one)
  ic(i) = 0
END DO
DO  i=1,nparm
  mmax(i) = (mmax(i) - mmin(i))/two
  mmin(i) = mmin(i) + mmax(i)
END DO
alfa = 0.01_dp
nfe = 0
ng = 0
ns = 0
nc = 0
ncp = 1
n0 = 0
n1 = 0
im = 1
ig = 0
fm = HUGE(one)
maxfn = 500*nparm
relcon = ten**(-nsig)

!         SAMPLING
20 n0 = n0 + n100
nm = n0 - 1
ng = ng + ng0
ns = ns + 1
IF (ns*ng0 > 100) GO TO 465
b = (one - alfa**(one/REAL(nm)))**b1
bb = 0.1*b
DO  i1=1,n
  DO  j=1,nn100
    CALL RANDOM_NUMBER( r )
    DO  i=1,nparm
      y(i) = two*r(i) - one
    END DO
    CALL fun(y, fc, nparm, m, mmin, mmax)
    IF (fc >= fm) CYCLE
    f(im) = fc
    DO  i=1,nparm
      x(i,im) = y(i)
    END DO
    IF (im <= ng .AND. ic(im) > 0) ig = ig - 1
    ic(im) = 0
    im = 1
    fm = f(1)
    DO  i=2,ng10
      IF (f(i) < fm) CYCLE
      im = i
      fm = f(i)
    END DO
  END DO
END DO

nfe = nfe + n100
WRITE(ipr,901) n100
WRITE(*,901) n100
901 FORMAT(/' ', i5, ' FUNCTION EVALUATIONS USED FOR SAMPLING')

!        SORTING
inum = ng10 - 1
DO  i=1,inum
  im = i
  fm = f(i)
  inum1 = i + 1
  DO  j=inum1,ng10
    IF (f(j) >= fm) CYCLE
    im = j
    fm = f(j)
  END DO
  IF (im <= i) CYCLE
  a = fm
  DO  j=1,nparm
    y(j) = x(j,im)
  END DO
  IF (i > ng .OR. im <= ng) GO TO 55
  IF (ic(ng) == 0 .AND. ic(im) > 0) ig = ig+1
  IF (ic(ng) > 0 .AND. ic(im) == 0) ig = ig-1
  55 icc = ic(im)
  inum1 = im-i
  DO  j=1,inum1
    inum2 = im-j
    f(inum2+1) = f(inum2)
    ic(inum2+1) = ic(inum2)
    DO  jj=1,nparm
      x(jj,inum2+1) = x(jj,inum2)
    END DO
  END DO
  f(i) = a
  DO  j=1,nparm
    x(j,i) = y(j)
  END DO
  ic(i) = icc
END DO
IF (nc <= 0) GO TO 200

!     CLUSTERING TO X*
DO  iii=1,nc
  i = 1
  in1 = i
  fcl(i) = f0(iii)
  DO  j=1,nparm
    xcl(j,i) = x0(j,iii)
  END DO
  DO  j=1,ng
    IF (ic(j) /= iii) CYCLE
    in1 = in1+1
    DO  ii=1,nparm
      xcl(ii,in1) = x(ii,j)
    END DO
  END DO
  95 DO  j=1,ng
    IF (ic(j) /= 0) CYCLE
    IF (fcl(i) >= f(j)) CYCLE
    DO  l1=1,nparm
      w(l1) = ABS(xcl(l1,i)-x(l1,j))
    END DO
    a = zero
    DO  l1=1,nparm
      IF (w(l1) > a) a = w(l1)
    END DO
    IF (a >= b) CYCLE
    WRITE(ipr,902) iii
    WRITE(*,902) iii
    902 FORMAT(' SAMPLE POINT ADDED TO THE CLUSTER NO. ', i2)
    DO  ii=1,nparm
      w(ii) = x(ii,j)*mmax(ii) + mmin(ii)
    END DO
    WRITE(ipr,903) f(j), w(1:nparm)
    WRITE (*,903) f(j), w(1:nparm)
    903 FORMAT(' ', g14.8/ ('    ', 5(g14.8, ' ')))
    ig = ig+1
    IF (ig >= ng) GO TO 395
    in1 = in1+1
    fcl(in1) = f(j)
    DO  ii=1,nparm
      xcl(ii,in1) = x(ii,j)
    END DO
    ic(j) = iii
  END DO
  i = i+1
  IF (i <= in1) GO TO 95
END DO
IF (n1 <= 0) GO TO 200

!     CLUSTERING TO X1
DO  iii=1,n1
  i = 1
  in1 = i
  fcl(i) = f1(iii)
  DO  j=1,nparm
    xcl(j,i) = x1(j,iii)
  END DO
  155 DO  j=1,ng
    IF (ic(j) /= 0) CYCLE
    IF (fcl(i) >= f(j)) CYCLE
    DO  l1=1,nparm
      w(l1) = ABS(xcl(l1,i)-x(l1,j))
    END DO
    a = zero
    DO  l1=1,nparm
      IF (w(l1) > a) a=w(l1)
    END DO
    IF (a >= b) CYCLE
    WRITE(ipr,902) ic1(iii)
    WRITE(*,902) ic1(iii)
    DO  ii=1,nparm
      w(ii) = x(ii,j)*mmax(ii) + mmin(ii)
    END DO
    WRITE(ipr,903) f(j), w(1:nparm)
    WRITE(*,903) f(j), w(1:nparm)
    ig = ig+1
    IF (ig >= ng) GO TO 395
    in1 = in1+1
    fcl(in1) = f(j)
    DO  ii=1,nparm
      xcl(ii,in1) = x(ii,j)
    END DO
    ic(j) = ic1(iii)
  END DO
  i = i+1
  IF (i <= in1) GO TO 155
END DO

!     LOCAL SEARCH
200 it = 0
DO  i1=1,ng
  IF (ic(i1) /= 0) CYCLE
  DO  i=1,nparm
    y(i) = x(i,i1)
  END DO
  ff = f(i1)
  CALL local(m, nparm, relcon, maxfn, y, ff, nfe1, mmin, mmax)
  IF (nc <= 0) GO TO 290
  DO  iv=1,nc
    DO  l1=1,nparm
      w(l1) = ABS(x0(l1,iv) - y(l1))
    END DO
    a = zero
    DO  l1=1,nparm
      IF (w(l1) > a) a = w(l1)
    END DO
    IF (a < bb) GO TO 255
  END DO
  GO TO 290

!       NEW SEED-POINT
  255 n1 = n1 + 1
  WRITE(ipr,905) iv, nfe1
  WRITE(*,905) iv,nfe1
  905 FORMAT(' NEW SEED POINT ADDED TO THE CLUSTER NO. ', i2, ', NFEV=', i5)
  DO  ii=1,nparm
    w(ii) = x(ii,i1)*mmax(ii) + mmin(ii)
  END DO
  WRITE(ipr,903) ff, w(1:nparm)
  WRITE(*,903) ff, w(1:nparm)
  IF (ff >= f0(iv)) GO TO 280
  WRITE(ipr,906) iv, f0(iv), ff
  WRITE(*,906) iv,f0(iv),ff
  906 FORMAT(' *** IMPROVEMENT ON THE LOCAL MINIMUM NO. ',  &
             i2, ':', g14.8, ' FOR ', g14.8)
  w(1:nparm) = y(1:nparm)*mmax(1:nparm) + mmin(1:nparm)
  WRITE(ipr,903) ff, w(1:nparm)
  WRITE(*,903) ff, w(1:nparm)
  f0(iv) = ff
  DO  ii=1,nparm
    x0(ii,iv) = y(ii)
  END DO
  280 IF (n1 > 20) GO TO 470
  DO  ii=1,nparm
    x1(ii,n1) = x(ii,i1)
    xcl(ii,1) = x(ii,i1)
  END DO
  f1(n1) = f(i1)
  fcl(1) = f(i1)
  ic1(n1) = iv
  icj = iv
  GO TO 305

!     NEW LOCAL MINIMUM
  290 nc = nc+1
  ncp = ncp+1
  WRITE(ipr,907) nc, ff, nfe1
  WRITE(*,907) nc, ff, nfe1
  907 FORMAT(' *** THE LOCAL MINIMUM NO. ', i2, ': ', g14.8, ', NFEV=', i5)
  DO  ii=1,nparm
    w(ii) = y(ii)*mmax(ii) + mmin(ii)
  END DO
  WRITE(ipr,903) ff, w(1:nparm)
  WRITE(*,903) ff, w(1:nparm)
  DO  ii=1,nparm
    x0(ii,nc) = y(ii)
    xcl(ii,1) = y(ii)
  END DO
  fcl(1) = ff
  f0(nc) = ff
  IF (nc >= 20) GO TO 475
  it = 1
  icj = nc

!    CLUSTERING TO THE NEW POINT
  305 nfe = nfe + nfe1
  ic(i1) = icj
  ig = ig + 1
  IF (ig >= ng) EXIT
  i = 1
  in1 = i
  310 DO  j=1,ng
    IF (ic(j) /= 0) CYCLE
    IF (fcl(i) >= f(j)) CYCLE
    DO  l1=1,nparm
      w(l1) = ABS(xcl(l1,i) - x(l1,j))
    END DO
    a = zero
    DO  l1=1,nparm
      IF (w(l1) > a) a = w(l1)
    END DO
    IF (a >= b) CYCLE
    in1 = in1 + 1
    DO  ii=1,nparm
      xcl(ii,in1) = x(ii,j)
    END DO
    fcl(in1) = f(j)
    ic(j) = icj
    WRITE(ipr,902) icj
    WRITE(*,902) icj
    DO  ii=1,nparm
      w(ii) = x(ii,j)*mmax(ii) + mmin(ii)
    END DO
    WRITE(ipr,903) f(j), w(1:nparm)
    WRITE(*,903) f(j), w(1:nparm)
    ig = ig+1
    IF (ig >= ng) EXIT
  END DO
  i = i+1
  IF (i < in1) GO TO 310
END DO
IF (it /= 0) GO TO 20

!      PRINT RESULTS
395 WRITE(ipr,908)
WRITE(*,908)
908 FORMAT(///' LOCAL MINIMA FOUND:'//)
IF (nc <= 1) GO TO 430
inum = nc - 1
DO  i=1,inum
  im = i
  fm = f0(i)
  inum1 = i + 1
  DO  j=inum1,nc
    IF (f0(j) >= fm) CYCLE
    im = j
    fm = f0(j)
  END DO
  IF (im <= i) CYCLE
  a = fm
  y(1:nparm) = x0(1:nparm,im)
  inum1 = im-i
  DO  j=1,inum1
    inum2 = im - j
    f0(inum2+1) = f0(inum2)
    DO  jj=1,nparm
      x0(jj,inum2+1) = x0(jj,inum2)
    END DO
  END DO
  f0(i) = a
  x0(1:nparm,i) = y(1:nparm)
END DO
430 IF (nc <= 0) GO TO 445

DO  i=1,nc
  x0(1:nparm,i) = x0(1:nparm,i)*mmax(1:nparm) + mmin(1:nparm)
  WRITE(ipr,903) f0(i), x0(1:nparm,i)
  WRITE(*,903) f0(i), x0(1:nparm,i)
END DO

445 WRITE(ipr,911) nfe
WRITE(*,911) nfe
911 FORMAT(//' NORMAL TERMINATION AFTER ', i5, ' FUNCTION EVALUATIONS'//)
RETURN

460 WRITE(ipr,914)
WRITE(*,914)
914 FORMAT(' ***   DATA ERROR')
STOP

465 WRITE(ipr,915)
WRITE(*,915)
915 FORMAT(' ***   TOO MANY SAMPLES')
GO TO 395

470 WRITE(ipr,916)
WRITE(*,916)
916 FORMAT(' ***   TOO MANY NEW SEED POINTS')
GO TO 395

475 WRITE(ipr,917)
WRITE(*,917)
917 FORMAT(' ***   TOO MANY CLUSTERS')
GO TO 395
END SUBROUTINE global



SUBROUTINE fun(r, f, nparm, m, mmin, mmax)

REAL (dp), INTENT(IN)   :: r(:)
REAL (dp), INTENT(OUT)  :: f
INTEGER, INTENT(IN)     :: nparm
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: mmin(:)
REAL (dp), INTENT(IN)   :: mmax(:)

INTERFACE
  SUBROUTINE funct(x, f, nparm, m)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f
    INTEGER, INTENT(IN)     :: nparm, m
  END SUBROUTINE funct
END INTERFACE

REAL (dp) :: x(nparm)

! N.B. mmin = mid-point between lower & upper limits
!      mmax = (upper - lower limits) / 2

x(1:nparm) = mmax(1:nparm)*r(1:nparm) + mmin(1:nparm)
CALL funct(x, f, nparm, m)

RETURN
END SUBROUTINE fun



!   ROUTINE NAME - LOCAL
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-27  Time: 13:16:43

!-----------------------------------------------------------------------

!   LATEST REVISION - JULY 31, 1986

!   PURPOSE  - MINIMUM OF A FUNCTION OF N VARIABLES USING
!       A QUASI-NEWTON METHOD

!   USAGE  - CALL LOCAL (M, N, EPS, MAXFN, X, F, NFEV, mmin, mmax)

!   ARGUMENTS  M - THE NUMBER OF RESIDUAL FUNCTIONS (INPUT)
!       NOT USED IN THIS ROUTINE.
!   N - THE NUMBER OF PARAMETERS (I.E., THE LENGTH OF X) (INPUT)
!   EPS - CONVERGENCE CRITERION. (INPUT).  THE ACCURACY REQUIRED IN THE
!       PARAMETER ESTIMATES.
!       THIS CONVERGENCE CONDITION IS SATISFIED IF ON TWO SUCCESSIVE
!       ITERATIONS, THE PARAMETER ESTIMATES (I.E.,X(I), I=1,...,N) DIFFERS,
!       COMPONENT BY COMPONENT, BY AT MOST EPS.
!   MAXFN - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E.,
!       CALLS TO SUBROUTINE FUN) ALLOWED. (INPUT)
!   X - VECTOR OF LENGTH N CONTAINING PARAMETER VALUES.
!     ON INPUT, X MUST CONTAIN THE INITIAL PARAMETER ESTIMATES.
!     ON OUTPUT, X CONTAINS THE FINAL PARAMETER
!       ESTIMATES AS DETERMINED BY LOCAL.
!   F - A SCALAR CONTAINING THE VALUE OF THE FUNCTION
!       AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)
!   NFEV - THE NUMBER OF FUNCTION EVALUATIONS (OUTPUT)
!   W - A VECTOR OF LENGTH 3*N USED AS WORKING SPACE.
!   mmin    - A VECTOR OF LENGTH N CONTAINING THE LOWER
!              BOUNDS OF THE PARAMETERS, SO X(I) IS
!              SEARCHED IN THE INTERVAL (mmin(I),mmax(I)).
!              (INPUT)
!   mmax    - A VECTOR OF LENGTH N CONTAINING THE UPPER
!              BOUNDS OF THE PARAMETERS. (INPUT)

!   REQUIRED ROUTINES - UPDATE, FUN

!   FUN   - A USER SUPPLIED SUBROUTINE WHICH CALCULATES THE FUNCTION F FOR
!           GIVEN PARAMETER VALUES X(1),X(2),...,X(N).
!           THE CALLING SEQUENCE HAS THE FOLLOWING FORM
!              CALL FUN(X, F, N, M, mmin, mmax)
!           WHERE X IS A VEKTOR OF LENGTH N.
!           FUN MUST APPEAR IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM.
!           FUN MUST NOT ALTER THE VALUES OF X(I), I=1,...,N OR N.

!-----------------------------------------------------------------------

SUBROUTINE local (m, n, eps, maxfn, x, f, nfev, mmin, mmax)

! N.B. Argument W has been removed.

!       SPECIFICATIONS FOR ARGUMENTS

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: eps
INTEGER, INTENT(IN)        :: maxfn
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(OUT)     :: f
INTEGER, INTENT(OUT)       :: nfev
REAL (dp), INTENT(IN)      :: mmin(:)
REAL (dp), INTENT(IN)      :: mmax(:)


!       SPECIFICATIONS FOR LOCAL VARIABLES
REAL (dp) :: w(3*n)
INTEGER   :: ig, igg, im1, is, idiff, ir, ij, i, iopt, j, nm1, jj, jp1, l, &
             kj, k, link, itn, ii, jnt, np1, jb, nj, ier
REAL (dp) :: hh, hjj, v, df, relx, gs0, diff, aeps, alpha, ff, tot, f1, f2, &
             z, gys, dgs, sig, zz, hhh, ghh, g(n), h(n*n)
REAL (dp), PARAMETER :: reps = 1.1921E-07_dp, ax = 0.1_dp, zero = 0.0_dp,  &
                        one = 1.0_dp, half = 0.5_dp, seven = 7.0_dp,  &
                        five = 5.0_dp, twelve = 12.0_dp, p1 = 0.1_dp

!       INITIALIZATION
!       FIRST EXECUTABLE STATEMENT
iopt = 0

!   IOPT - OPTIONS SELECTOR. (INPUT)
!     IOPT = 0 CAUSES LOCAL TO INITIALIZE THE
!       HESSIAN MATRIX H TO THE IDENTITY MATRIX.
!     IOPT = 1 INDICATES THAT H HAS BEEN INITIALIZED
!       BY THE USER TO A POSITIVE DEFINITE MATRIX.
!     IOPT = 2 CAUSES LOCAL TO COMPUTE THE DIAGONAL
!       VALUES OF THE HESSIAN MATRIX AND SET H TO
!       A DIAGONAL MATRIX CONTAINING THESE VALUES.
!     IOPT = 3 CAUSES LOCAL TO COMPUTE AN ESTIMATE
!       OF THE HESSIAN IN H.
ier = 0
hh = SQRT(reps)
ig = n
igg = n+n
is = igg
idiff = 1
ir = n
w(1) = -one
w(2) = zero
w(3) = zero

!       EVALUATE FUNCTION AT STARTING POINT
g(1:n) = x(1:n)
CALL fun (g, f, n, m, mmin, mmax)
nfev = 1
IF (iopt == 1) GO TO 45

!       SET OFF-DIAGONAL ELEMENTS OF H TO 0.0
IF (n == 1) GO TO 25
ij = 2
DO  i=2,n
  DO  j=2,i
    h(ij) = zero
    ij = ij + 1
  END DO
  ij = ij + 1
END DO
IF (iopt /= 0) GO TO 25

!       SET DIAGONAL ELEMENTS OF H TO ONE
ij = 0
DO  i=1,n
  ij = ij+i
  h(ij) = one
END DO
GO TO 80

!       GET DIAGONAL ELEMENTS OF HESSIAN
25 im1 = 1
nm1 = 1
np1 = n+1
DO  i=2,np1
  hhh = hh*MAX(ABS(x(im1)), ax)
  g(im1) = x(im1) + hhh
  CALL fun (g, f2, n, m, mmin, mmax)
  g(im1) = g(im1) + hhh
  CALL fun (g, ff, n, m, mmin, mmax)
  h(nm1) = (ff-f2+f-f2)/(hhh*hhh)
  g(im1) = x(im1)
  im1 = i
  nm1 = i+nm1
END DO
nfev = nfev+n+n
IF (iopt /= 3 .OR. n == 1) GO TO 45

!       GET THE REST OF THE HESSIAN
jj = 1
ii = 2
DO  i=2,n
  ghh = hh*MAX(ABS(x(i)), ax)
  g(i) = x(i)+ghh
  CALL fun (g, f2, n, m, mmin, mmax)
  DO  j=1,jj
    hhh = hh*MAX(ABS(x(j)), ax)
    g(j) = x(j) + hhh
    CALL fun (g, ff, n, m, mmin, mmax)
    g(i) = x(i)
    CALL fun (g, f1, n, m, mmin, mmax)
!     H(II) = (FF-F1-F2+F)*SQREPS
    h(ii) = (ff-f1-f2+f)/(hhh*ghh)
    ii = ii+1
    g(j) = x(j)
  END DO
  jj = jj+1
  ii = ii+1
END DO
nfev = nfev + ((n*n-n)/2)

!       FACTOR H TO L*D*L-TRANSPOSE
45 ir = n
IF (n > 1) GO TO 50
IF (h(1) > zero) GO TO 80
h(1) = zero
ir = 0
GO TO 75
50 nm1 = n-1
jj = 0
DO  j=1,n
  jp1 = j+1
  jj = jj+j
  hjj = h(jj)
  IF (hjj > zero) GO TO 55
  h(jj) = zero
  ir = ir-1
  CYCLE
  55 IF (j == n) CYCLE
  ij = jj
  l = 0
  DO  i=jp1,n
    l = l+1
    ij = ij+i-1
    v = h(ij)/hjj
    kj = ij
    DO  k=i,n
      h(kj+l) = h(kj+l) - h(kj)*v
      kj = kj + k
    END DO
    h(ij) = v
  END DO
END DO
75 IF (ir /= n) THEN
  ier = 129
  GO TO 9000
END IF
80 itn = 0
df = -one

!       EVALUATE GRADIENT W(IG+I),I=1,...,N
85 link = 1
GO TO 260

!       BEGIN ITERATION LOOP
90 IF (nfev >= maxfn) GO TO 225
itn = itn+1
DO  i=1,n
  w(i) = -w(ig+i)
END DO

!       DETERMINE SEARCH DIRECTION W
!         BY SOLVING H*W = -G WHERE
!         H = L*D*L-TRANSPOSE
IF (ir < n) GO TO 125
!       N .EQ. 1
g(1) = w(1)
IF (n > 1) GO TO 100
w(1) = w(1)/h(1)
GO TO 125
!       N > 1
100 ii = 1

!       SOLVE L*W = -G
DO  i=2,n
  ij = ii
  ii = ii+i
  v = w(i)
  im1 = i-1
  DO  j=1,im1
    ij = ij+1
    v = v - h(ij)*w(j)
  END DO
  g(i) = v
  w(i) = v
END DO

!       SOLVE (D*LT)*Z = W WHERE
!                                    LT = L-TRANSPOSE
w(n) = w(n)/h(ii)
jj = ii
nm1 = n-1
DO  nj=1,nm1
!       J = N-1,N-2,...,1
  j = n-nj
  jp1 = j+1
  jj = jj-jp1
  v = w(j)/h(jj)
  ij = jj
  DO  i=jp1,n
    ij = ij+i-1
    v = v - h(ij)*w(i)
  END DO
  w(j) = v
END DO

!       DETERMINE STEP LENGTH ALPHA
125 relx = zero
gs0 = zero
DO  i=1,n
  w(is+i) = w(i)
  diff = ABS(w(i)) / MAX(ABS(x(i)),ax)
  relx = MAX(relx, diff)
  gs0 = gs0 + w(ig+i)*w(i)
END DO
IF (relx == zero) GO TO 230
aeps = eps/relx
ier = 130
IF (gs0 >= zero) GO TO 230
IF (df == zero) GO TO 230
ier = 0
alpha = (-df-df)/gs0
IF (alpha <= zero) alpha = one
alpha = MIN(alpha,one)
IF (idiff == 2) alpha = MAX(p1,alpha)
ff = f
tot = zero
jnt = 0

!       SEARCH ALONG X + ALPHA*W
135 IF (nfev >= maxfn) GO TO 225
DO  i=1,n
  w(i) = x(i) + alpha*w(is+i)
END DO
CALL fun (w, f1, n, m, mmin, mmax)
nfev = nfev+1
IF (f1 >= f) GO TO 165
f2 = f
tot = tot + alpha
145 ier = 0
f = f1
x(1:n) = w(1:n)
IF (jnt-1 < 0) THEN
  GO TO   155
ELSE IF (jnt-1 == 0) THEN
  GO TO 185
ELSE
  GO TO 190
END IF

155 IF (nfev >= maxfn) GO TO 225
DO  i=1,n
  w(i) = x(i) + alpha*w(is+i)
END DO
CALL fun (w, f1, n, m, mmin, mmax)
nfev = nfev+1
IF (f1 >= f) GO TO 190
IF (f1+f2 >= f+f .AND. seven*f1+five*f2 > twelve*f) jnt = 2
tot = tot + alpha
alpha = alpha + alpha
GO TO 145

165 IF (f == ff .AND. idiff == 2 .AND. relx > eps) ier = 130
IF (alpha < aeps) GO TO 230
IF (nfev >= maxfn) GO TO 225
alpha = half*alpha
DO  i=1,n
  w(i) = x(i) + alpha*w(is+i)
END DO
CALL fun (w, f2, n, m, mmin, mmax)
nfev = nfev + 1
IF (f2 >= f) GO TO 180
tot = tot + alpha
ier = 0
f = f2
x(1:n) = w(1:n)
GO TO 185

180 z = p1
IF (f1+f > f2+f2) z = one + half*(f-f1)/(f+f1-f2-f2)
z = MAX(p1,z)
alpha = z*alpha
jnt = 1
GO TO 135

185 IF (tot < aeps) GO TO 230
190 alpha = tot

!       SAVE OLD GRADIENT
DO  i=1,n
  w(i) = w(ig+i)
END DO

!       EVALUATE GRADIENT W(IG+I), I=1,...,N
link = 2
GO TO 260
200 IF (nfev >= maxfn) GO TO 225
gys = zero
DO  i=1,n
  gys = gys + w(ig+i)*w(is+i)
  w(igg+i) = w(i)
END DO
df = ff-f
dgs = gys-gs0
IF (dgs <= zero) GO TO 90
IF (dgs + alpha*gs0 > zero) GO TO 215

!       UPDATE HESSIAN H USING
!         COMPLEMENTARY DFP FORMULA
sig = one/gs0
ir = -ir
CALL update (h, n, w, sig, g, ir, 0, zero)
DO  i=1,n
  g(i) = w(ig+i) - w(igg+i)
END DO
sig = one/(alpha*dgs)
ir = -ir
CALL update (h, n, g, sig, w, ir, 0, zero)
GO TO 90

!       UPDATE HESSIAN USING DFP FORMULA
215 zz = alpha/(dgs - alpha*gs0)
sig = -zz
CALL update (h, n, w, sig, g, ir, 0, reps)
z = dgs*zz - one
DO  i=1,n
  g(i) = w(ig+i) + z*w(igg+i)
END DO
sig = one/(zz*dgs*dgs)
CALL update (h, n, g, sig, w, ir, 0, zero)
GO TO 90

!       MAXFN FUNCTION EVALUATIONS
225 GO TO 235
230 IF (idiff == 2) GO TO 235

!       CHANGE TO CENTRAL DIFFERENCES
idiff = 2
GO TO 85
235 IF (relx > eps .AND. ier == 0) GO TO 85

!       COMPUTE H = L*D*L-TRANSPOSE AND OUTPUT
IF (n == 1) GO TO 9000
np1 = n+1
nm1 = n-1
jj = (n*(np1))/2
DO  jb=1,nm1
  jp1 = np1-jb
  jj = jj-jp1
  hjj = h(jj)
  ij = jj
  l = 0
  DO  i=jp1,n
    l = l+1
    ij = ij+i-1
    v = h(ij)*hjj
    kj = ij
    DO  k=i,n
      h(kj+l) = h(kj+l) + h(kj)*v
      kj = kj + k
    END DO
    h(ij) = v
  END DO
  hjj = h(jj)
END DO
GO TO 9000

!        EVALUATE GRADIENT
260 IF (idiff == 2) GO TO 270

!       FORWARD DIFFERENCES
!         GRADIENT = W(IG+I), I=1,...,N
DO  i=1,n
  z = hh*MAX(ABS(x(i)),ax)
  zz = x(i)
  x(i) = zz + z
  CALL fun (x, f1, n, m, mmin, mmax)
  w(ig+i) = (f1-f)/z
  x(i) = zz
END DO
nfev = nfev+n
SELECT CASE ( link )
  CASE (    1)
    GO TO 90
  CASE (    2)
    GO TO 200
END SELECT

!       CENTRAL DIFFERENCES
!         GRADIENT = W(IG+I), I=1,...,N
270 DO  i=1,n
  z = hh*MAX(ABS(x(i)), ax)
  zz = x(i)
  x(i) = zz + z
  CALL fun (x, f1, n, m, mmin, mmax)
  x(i) = zz-z
  CALL fun (x, f2, n, m, mmin, mmax)
  w(ig+i) = (f1-f2)/(z+z)
  x(i) = zz
END DO
nfev = nfev+n+n
SELECT CASE ( link )
  CASE (    1)
    GO TO 90
  CASE (    2)
    GO TO 200
END SELECT

!       RETURN
9000 RETURN
END SUBROUTINE local



!   ROUTINE NAME - UPDATE

!-----------------------------------------------------------------------

!   LATEST REVISION     - JULY 31, 1986
!     (CHANGES IN COMMENTS)

!   PURPOSE  - NUCLEUS CALLED ONLY BY ROUTINE LOCAL

!   REQD. ROUTINES - NONE REQUIRED

!-----------------------------------------------------------------------

SUBROUTINE update (a, n, z, sig, w, ir, mk, eps)
!       SPECIFICATIONS FOR ARGUMENTS

REAL (dp), INTENT(OUT)     :: a(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: z(:)
REAL (dp), INTENT(IN)      :: sig
REAL (dp), INTENT(IN OUT)  :: w(:)
INTEGER, INTENT(OUT)       :: ir
INTEGER, INTENT(IN)        :: mk
REAL (dp), INTENT(IN)      :: eps


!       SPECIFICATIONS FOR LOCAL VARIABLES
INTEGER   :: j, jj, ij, jp1, i, ii, mm
REAL (dp) :: ti, v, tim, al, r, b, gm, y
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, four = 4.0_dp

!       UPDATE FACTORS GIVEN IN A
!         SIG*Z*Z-TRANSPOSE IS ADDED
!       FIRST EXECUTABLE STATEMENT
IF (n > 1) GO TO 5

!       N .EQ. 1
a(1) = a(1) + sig*z(1)*z(1)
ir = 1
IF (a(1) > zero) GO TO 9005
a(1) = zero
ir = 0
GO TO 9005

!       N > 1
5 IF (sig > zero) GO TO 65
IF (sig == zero .OR. ir == 0) GO TO 9005
ti = one/sig
jj = 0
IF (mk == 0) GO TO 15

!       L*W = Z ON INPUT
DO  j=1,n
  jj = jj+j
  IF (a(jj) /= zero) ti = ti + (w(j)*w(j))/a(jj)
END DO
GO TO 40

!       SOLVE L*W = Z
15 w(1:n) = z(1:n)
DO  j=1,n
  jj = jj+j
  v = w(j)
  IF (a(jj) > zero) GO TO 25
  w(j) = zero
  CYCLE
  25 ti = ti+(v*v)/a(jj)
  IF (j == n) CYCLE
  ij = jj
  jp1 = j+1
  DO  i=jp1,n
    ij = ij+i-1
    w(i) = w(i) - v*a(ij)
  END DO
END DO

!        SET TI, TIM AND W
40 IF (ir <= 0) GO TO 45
IF (ti > zero) GO TO 50
IF (mk-1 > 0) THEN
  GO TO  55
ELSE
  GO TO  65
END IF

45 ti = zero
ir = -ir-1
GO TO 55

50 ti = eps/sig
IF (eps == zero) ir = ir-1
55 tim = ti
ii = jj
i = n
DO  j=1,n
  IF (a(ii) /= zero) tim = ti - (w(i)*w(i))/a(ii)
  w(i) = ti
  ti = tim
  ii = ii-i
  i = i-1
END DO
mm = 1
GO TO 70

65 mm = 0
tim = one/sig
70 jj = 0

!       UPDATE A
DO  j=1,n
  jj = jj+j
  ij = jj
  jp1 = j+1

!       UPDATE A(J,J)
  v = z(j)
  IF (a(jj) > zero) GO TO 85

!       A(J,J) .EQ. ZERO
  IF (ir > 0 .OR. sig < zero .OR. v == zero) GO TO 80
  ir = 1-ir
  a(jj) = (v*v)/tim
  IF (j == n) GO TO 9005
  DO  i=jp1,n
    ij = ij+i-1
    a(ij) = z(i)/v
  END DO
  GO TO 9005
  80 ti = tim
  CYCLE

!       A(J,J) .GT. ZERO
  85 al = v/a(jj)
  ti = w(j)
  IF (mm == 0) ti = tim + v*al
  r = ti/tim
  a(jj) = r*a(jj)
  IF (r == zero) EXIT
  IF (j == n) EXIT

!       UPDATE REMAINDER OF COLUMN J
  b = al/ti
  IF (r > four) GO TO 95
  DO  i=jp1,n
    ij = ij+i-1
    z(i) = z(i) - v*a(ij)
    a(ij) = a(ij) + b*z(i)
  END DO
  GO TO 105
  95 gm = tim/ti
  DO  i=jp1,n
    ij = ij+i-1
    y = a(ij)
    a(ij) = b*z(i) + y*gm
    z(i) = z(i) - v*y
  END DO
  105  tim = ti
END DO
IF (ir < 0) ir = -ir

9005 RETURN
END SUBROUTINE update

END MODULE global_minimum
