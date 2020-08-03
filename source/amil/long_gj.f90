SUBROUTINE gaussj(a, n, b, m)

!  Gauss-Jordan solution of linear equations in quadruple precision.
!  Adapted from the routine GAUSSJ in the Numerical Recipes package.
!  A is an NxN quadruple-precision matrix, declared as quadruple-precision
!  A(:,:).
!  B contains M right-hand side vectors, where M = 1 in many cases,
!  and is declared as quadruple-precision B(:,:).
!  On output, A has been overwritten by its inverse, and B has been
!  replaced with the solutions.

! Programmer : Alan Miller, CSIRO Mathematical Information Sciences
! e-mail: alan @ mel.dms.csiro.au   WWW: http://www.mel.dms.csiro.au/~alan
! Fax: (+61) 3-9545-8080

! Latest revision - 6 December 1996

USE quadruple_precision
IMPLICIT NONE

INTEGER, INTENT(IN)                         :: n, m
TYPE (quad), INTENT(IN OUT), DIMENSION(:,:) :: a, b

!       Local variables

INTEGER, PARAMETER :: nmax = 50
INTEGER            :: ipiv(nmax), indxr(nmax), indxc(nmax), j, i, irow,  &
                      icol, k, l, ll
TYPE (quad)        :: pivinv, dum, one, zero
REAL (dp)          :: big

one%hi = 1._dp
one%lo = 0._dp
zero%hi = 0._dp
zero%lo = 0._dp

ipiv(1:n) = 0

!       Start of the main loop

DO i = 1, n
  big = 0.
  DO j = 1, n
    IF(ipiv(j) /= 1) THEN
      DO k = 1, n
        IF (ipiv(k) == 0) THEN
          IF (ABS(a(j,k)%hi) >= big) THEN
            big = ABS(a(j,k)%hi)
            irow = j
            icol = k
          END IF
        ELSE IF (ipiv(k) > 1) THEN
          WRITE(*, *) 'Singular matrix'
          RETURN
        END IF
      END DO
    END IF
  END DO
  ipiv(icol) = ipiv(icol)+1

!       Pivot has been selected.   Interchange rows, if necessary, to
!       put the element on the diagonal.

  IF (irow /= icol) THEN
    DO l = 1, n
      dum = a(irow,l)
      a(irow,l) = a(icol,l)
      a(icol,l) = dum
    END DO
    DO l = 1, m
      dum = b(irow,l)
      b(irow,l) = b(icol,l)
      b(icol,l) = dum
    END DO
  END IF
  indxr(i) = irow
  indxc(i) = icol

!       Now do the pivoting

  IF (a(icol,icol)%hi == 0.) THEN
    WRITE(*, *) 'Singular matrix.'
    RETURN
  END IF
  pivinv = one / a(icol,icol)
  a(icol,icol) = one
  DO l = 1, n
    a(icol, l) = a(icol, l) * pivinv
  END DO
  DO l = 1, m
    b(icol, l) = b(icol, l) * pivinv
  END DO
  DO ll = 1, n
    IF(ll /= icol)THEN
    dum = a(ll,icol)
    a(ll,icol) = zero
    DO l = 1, n
      a(ll, l) = a(ll, l) - a(icol, l) * dum
    END DO
    DO l = 1, m
      b(ll, l) = b(ll, l) - b(icol, l) * dum
    END DO
    END IF
  END DO
END DO

!       Finally unscramble A by re-arranging the columns

DO l = n, 1, -1
  IF(indxr(l) /= indxc(l))THEN
    DO k = 1, n
      dum = a(k,indxr(l))
      a(k,indxr(l)) = a(k,indxc(l))
      a(k,indxc(l)) = dum
    END DO
  END IF
END DO

RETURN
END SUBROUTINE gaussj
