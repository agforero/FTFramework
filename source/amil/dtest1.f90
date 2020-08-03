PROGRAM dtest1
USE multidim_integrate
IMPLICIT NONE

INTEGER            :: key, n, mincls, maxcls, ifail, neval, nsub
INTEGER, PARAMETER :: ndim = 5, nf = ndim+1
REAL (dp)          :: a(ndim), b(ndim)
REAL (dp)          :: absest(nf), finest(nf), absreq, relreq

INTERFACE
  SUBROUTINE ftest(ndim, z, nfun, f)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
    INTEGER, INTENT(IN)    :: ndim, nfun
    REAL (dp), INTENT(IN)  :: z(:)
    REAL (dp), INTENT(OUT) :: f(:)
  END SUBROUTINE ftest
END INTERFACE

DO n = 1,ndim
  a(n) = 0
  b(n) = 1
END DO
mincls = 0
maxcls = 10000
key = 0
absreq = 0
relreq = 1.D-3
CALL dcuhre(ndim, nf, a, b, mincls, maxcls, ftest, absreq, relreq,  &
            key, 0, finest, absest, neval, ifail, nsub)
WRITE(*,9999) neval, ifail
9999 FORMAT ('        DCUHRE TEST RESULTS' // '     FTEST CALLS = ', i4,  &
	     ', IFAIL = ', i2 / '    N   ESTIMATED ERROR    INTEGRAL')
DO n = 1, nf
  WRITE(*, 9998) n, absest(n), finest(n)
  9998 FORMAT ("   ", i2, 2F15.6)
END DO
WRITE(*, '(a, i6)') ' No. of subregions generated = ', nsub

STOP
END PROGRAM dtest1


SUBROUTINE ftest(ndim, z, nfun, f)
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
INTEGER, INTENT(IN)    :: ndim, nfun
REAL (dp), INTENT(IN)  :: z(:)
REAL (dp), INTENT(OUT) :: f(:)

! Local variables
INTEGER   :: n
REAL (dp) :: sum

sum = 0
DO n = 1, ndim
  sum = sum + n*z(n)**2
END DO
f(1) = EXP(-sum/2)
f(2:nfun) = f(1) * z(1:nfun-1)

RETURN
END SUBROUTINE ftest
