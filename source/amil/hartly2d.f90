MODULE HartleyFFT_2D
IMPLICIT NONE
PUBLIC  :: fht2d

! From hamill@mailhub.nmg.sms.siemens.com Wed Apr 24 16:21:31 1996
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-02-04  Time: 13:01:18

! I see from the newsgroups that you are interested in multidimensional
! Hartley transforms.  The attached FORTRAN code works (I think ... but
! it's been a few years since I used it) for 2D and 1D.  For 3D see H.Hao
! and R.N.Bracewell, "A Three-Dimensional DFT Algorithm Using the Fast
! Hartley Transform," Proc IEEE vol. 75 no. 2 (Feb, 87) p.264-266.

! I do not vouch for the correctness of the information or the code, but
! perhaps it is a useful starting point.

! Good luck!

! Jim Hamill

!------------------------------------------------
! hamill@mailhub.nmg.sms.siemens.com
! James Hamill
! Siemens Medical Systems, Nuclear Medicine Group
! 2501 N. Barrington Rd.
! Hoffman Estates, IL  60195
! voice:(847)304-7760       FAX:(847)304-7706
!------------------------------------------------

CONTAINS


!------------------------------------------------------------------------------
!   FHT2D
!        A subroutine for the 2-dimensional Hartley transform of an
!        array  F2(N,N)  where N = 2**P
!        This routine call the 1-d Hartley transform routine FHTFR, for
!        rows and then for columns.
!-------------------------------------------------------------------------------

SUBROUTINE fht2d(nx, f2)
INTEGER, INTENT(IN)   :: nx
REAL, INTENT(IN OUT)  :: f2(0:nx-1,0:nx-1)

REAL     :: a, b, c, d, e, f(0:255)
INTEGER  :: i, j, n, ny, n1, p, u, v

p = INT( LOG(REAL(nx))/LOG(2.) + 0.5)
n  = 2**p
n1 = n - 1
DO j = 0,n1
  CALL fhtfr(p, f2(0,j))
END DO
DO i = 0,n1
  DO j = 0,n1
    f(j)=f2(i,j)
  END DO
  CALL fhtfr(p, f)
  DO j = 0,n1
    f2(i,j)=f(j)
  END DO
END DO

!  UP TO THIS POINT F2(U,V) NOW CONTAIN THE SEPARABLE
!                    CAS(2*pi*u*x)*CAS(2*pi*v*y)   transform

!  The following replacements form the nonseparable 2D-Hartley transform
!    Ref:   RN Bracewell, et al., "Fast Two-dimensional Hartley transform,"
!           Proc IEEE, 74:9,1282 (Sept 1986)

ny = nx
DO v = 1,ny/2-1
  DO u = 1,nx/2-1
    a = f2(u,v)
    b = f2(nx-u,v)
    c = f2(u,ny-v)
    d = f2(nx-u,ny-v)
    e = 0.5*( (a+d)-(b+c) )
    f2(u,v)       = a - e
    f2(nx-u,v)    = b + e
    f2(u,ny-v)    = c + e
    f2(nx-u,ny-v) = d - e
  END DO
END DO

RETURN
END SUBROUTINE fht2d


!-----------------------------------------------------------------------------
!   SUBROUTINE FHTFOR TAKES INPUT F() AND RETURNS THE DHT IN THE SAME F()
!   LENGTH OF F() IS 2**P
!-----------------------------------------------------------------------------
!     THIS IS THE CORRECTED VERSION OF FHTFOR.FOR
!     RN BRACEWELL,  THE HARTLEY TRANSFORM, Oxford Univ. Press (1986)
!     pp. 140-144
!-----------------------------------------------------------------------------

SUBROUTINE fhtfr(p, f)
INTEGER, INTENT(IN)   :: p
REAL, INTENT(IN OUT)  :: f(0:255)

REAL     :: h9, s9(256), t, t9(256)
INTEGER  :: d9, e9, i, j, k9, l9, n, n9, q9, sx9, u9, x9, y9,   &
    a9(0:64), m9(0:10), c9(0:10), v9(0:10)
INTEGER  :: pold = -99
REAL, PARAMETER     :: pi = 3.1415926536

IF (p == 1) THEN
  j    = f(0) + f(1)
  f(1) = f(0) - f(1)
  f(0) = j
  GO TO 9636   !   RETURN
END IF
IF(p == pold) GO TO 9400
n9 = 2**(p-2)
n = 4*n9
c9(5) = n-1
c9(6) = p-1

!   GET POWERS OF 2
i = 1
m9(0) = 1
m9(1) = 2
9204   m9(i+1) = m9(i) + m9(i)
i = i + 1
IF(i <= p) GO TO 9204

!  TRIGONOMETRIC FUNCTION TABLES
IF(n == 4) GO TO 9400
IF(n == 8) THEN
  s9(2) = 1.
  s9(1) = SIN(pi/4.)
  GO TO 9330
END IF
!   GET SINES
s9(n9) = 1.
!   COARSE SEED TABLE FOR SINES
DO i = 1,3
  s9(i*n9/4) = SIN(i*pi/8.)
END DO
!   INITIAL HALF SECANT
h9 = 1./2/COS(pi/16)
!   FILL SINE TABLE
c9(4) = p-4
DO i = 1,p-4
  c9(4) = c9(4)-1
  v9(0) = 0
  DO j = m9(c9(4)), n9-m9(c9(4)), m9(c9(4)+1)
    v9(1) = j + m9(c9(4))
    s9(j) = h9*(s9(v9(1))+v9(0))
    v9(0) = s9(v9(1))
  END DO
!   HALF SECANT RECURSION
  h9 = 1/SQRT(2+1/h9)
END DO

!   GET TANGENTS
9330 c9(0) = n9 - 1
DO i = 1,n9-1
  t9(i) = (1-s9(c9(0)))/s9(i)
  c9(0) = c9(0) - 1
END DO
t9(n9) = 1

!   FAST PERMUTE
!   FOR P = 2,3 PERMUTE DIRECTLY
9400 IF(p == 2) THEN
  v9(9) = f(1)
  f(1)  = f(2)
  f(2)  = v9(9)
  GO TO 9500
END IF
IF(p == 3) THEN
  v9(9) = f(1)
  f(1)  = f(4)
  f(4)  = v9(9)
  v9(9) = f(3)
  f(3)  = f(6)
  f(6)  = v9(9)
  GO TO 9500
END IF
IF(p == pold) GO TO 9420
!   FOR P=4, 5, 6 (Q9=2,3), SKIP STRUCTURE TABLE
q9 = p/2
c9(2) = m9(q9)
q9 = q9 + MOD(p,2)
IF(q9 == 2) THEN
  a9(0) = 0
  a9(1) = 2
  a9(2) = 1
  a9(3) = 3
  GO TO 9420
END IF
IF(q9 == 3) THEN
  a9(0) = 0
  a9(1) = 4
  a9(2) = 2
  a9(3) = 6
  a9(4) = 1
  a9(5) = 5
  a9(6) = 3
  a9(7) = 7
  GO TO 9420
END IF

!   SET UP STRUCTURE TABLE
a9(0) = 0
a9(1) = 1
DO i = 2,q9
  DO j = 0,m9(i-1)-1
    a9(j) = a9(j) + a9(j)
    a9(j+m9(i-1))=a9(j)+1
  END DO
END DO

!   PERMUTE
9420 DO i = 1,c9(2)-1
  v9(4)  = c9(2)*a9(i)
  v9(5)  = i
  v9(6)  = v9(4)
  v9(7)  = f(v9(5))
  f(v9(5)) = f(v9(6))
  f(v9(6)) = v9(7)
  DO j =1,a9(i)-1
    v9(5) = v9(5) + c9(2)
    v9(6) = v9(4) + a9(j)
    v9(7) = f(v9(5))
    f(v9(5)) = f(v9(6))
    f(v9(6)) = v9(7)
  END DO
END DO

!   STAGES 1 & 2
!   GET TWO-ELEMENT DHTs
9500 DO i = 0,n-2,2
  t      = f(i) - f(i+1)
  f(i)   = f(i) + f(i+1)
  f(i+1) = t
END DO

!   GET FOUR-ELEMENT DHTs
DO i = 0,n-4,4
  t      = f(i) - f(i+2)
  f(i)   = f(i) + f(i+2)
  f(i+2) = t
  t      = f(i+1)-f(i+3)
  f(i+1) = f(i+1)+f(i+3)
  f(i+3) = t
END DO

IF(p == 2) GO TO 9636  !  RETURN

!   STAGES 3, 4, ...
u9 = c9(6)
sx9 = 4
DO l9 = 2,c9(6)
  v9(2) = sx9 + sx9
  u9 = u9-1
  v9(3) = m9(u9-1)
  DO q9 = 0,c9(5),v9(2)
    i = q9
    d9 = i+sx9
    t     = f(i) - f(d9)
    f(i)  = f(i) + f(d9)
    f(d9) = t
    k9 = d9-1
    DO j=v9(3),n9,v9(3)
      i  = i  + 1
      d9 = i  + sx9
      e9 = k9 + sx9
      t  = f(d9) + f(e9)*t9(j)
      x9 = f(e9) - t*s9(j)
      y9 = x9*t9(j) + t
      t = f(i) + y9
      f(d9) = f(i)  - y9
      f(e9) = f(k9) + x9
      f(k9) = f(k9) - x9
      f(i)  = t
      k9 = k9 -1
    END DO
  END DO
  sx9= v9(2)
END DO

!   SAVE P SO TABLE GENERATION CAN BE SKIPPED ON NEXT CALL
9636 pold = p
RETURN
END SUBROUTINE fhtfr

END MODULE HartleyFFT_2D
