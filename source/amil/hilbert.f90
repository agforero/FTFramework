PROGRAM hilbert

! Demonstration program 1 using the 'quad' package.
! Inversion of 10x10 Hilbert matrix by Cholesky factorization.

! Matrix A has elements a(i,j) = 1/(i+j-1).
! To obtain full accuracy, the elements are all multiplied by
! 9 x 5 x 7 x 11 x 13 x 17 x 19  to make them into integers so
! that the input matrix is exact.   The elements of the inverse
! will be multiplied by this integer at the end.

! The elements of the inverse should all be integers.

! Programmer: Alan Miller

! Latest Fortran 77 revision - 13 September 1986
! Fortran 90 version - 3 December 1996
! Latest revision - 3 October 1999

USE quadruple_precision
IMPLICIT NONE

TYPE (quad) :: a(55), total, exact_start, exact, relerr
REAL (dp)   :: det, scale_factor = 14549535.d0
INTEGER     :: i, j, nrows = 10, ier, ipos, ik, k, kj, kjstart

INTERFACE
  SUBROUTINE qchol(a, nrows, ier)
    USE quadruple_precision
    IMPLICIT NONE
    TYPE (quad), DIMENSION(:), INTENT(IN OUT) :: a
    INTEGER, INTENT(IN)                       :: nrows
    INTEGER, INTENT(OUT)                      :: ier
  END SUBROUTINE qchol
END INTERFACE

! Enter the elements of the Hilbert matrix, storing them in a
! 1-dimensional array up to the diagonal in each row.

ipos = 1
DO i = 1, nrows
  DO j = 1, i
    a(ipos) = quad(scale_factor/(i+j-1), 0._dp)
    ipos = ipos + 1
  END DO
END DO

! Get the Cholesky factorization.

CALL qchol(a, nrows, ier)
IF (ier /= 0) THEN
  WRITE(*, *) 'Error in qchol, number:', ier
  STOP
END IF

! Calculate the determinant by multiplying together the diagonal elements
! of the Cholesky factor, squaring the result, and remembering that each
! element was multiplied by scale_factor.

ipos = 1
det = 1.d0
DO i = 1, nrows
  det = det * a(ipos)%hi
  ipos = ipos + i + 1
END DO
det = det * det / scale_factor**nrows
WRITE(*, 900) det
900 FORMAT(' Determinant = ', g20.12/)

! Now invert the triangular matrix.

ipos = 1
DO i = 1, nrows
  kjstart = 1
  DO j = 1, i
    kj = kjstart
    ik = ipos
    IF (j < i) THEN
      total = quad(0._dp, 0._dp)
      DO k = j, i-1
        total = total - a(ik) * a(kj)
        ik = ik + 1
        kj = kj + k
      END DO
      a(ipos) = total / a(ik)
      ipos = ipos + 1
    ELSE
      a(ipos) = quad(1._dp, 0._dp) / a(ipos)
      ipos = ipos + 1
    END IF
    kjstart = kjstart + j + 1
  END DO
END DO

! We formed A = LL', so the inverse is L**(-T).L**(-1),
! where L**(-T) denotes the transpose of the inverse.

ipos = 1
kjstart = 1
DO i = 1, nrows
  DO j = 1, i
    ik = ipos
    total = quad(0._dp, 0._dp)
    kj = kjstart
    DO k = i, nrows
      total = total + a(ik) * a(kj)
      ik = ik + k
      kj = kj + k
    END DO
    a(ipos) = total * scale_factor
    ipos = ipos + 1
  END DO
  kjstart = kjstart + i + 1
END DO

! Output the result.
! The exact values for an n x n matrix are:-

!  (-1)**(i+j).(n+i-1)!(n+j-1)!
!      -------------------------------------
!      (i+j-1)[(i-1)!(j-1)!]**2.(n-i)!(n-j)!

WRITE(*, 910)
910 FORMAT(' row col         calculated           exact        Rel. error')
ipos = 1
exact_start = quad(100._dp, 0._dp)
DO i = 1, nrows
  exact = exact_start
  DO j = 1, i
    relerr = (a(ipos) - exact) / exact
    WRITE(*, 920) i, j, a(ipos)%hi, exact%hi, relerr%hi
    920 FORMAT(2I4, 2F20.4, g14.5)
    exact = -exact * (10+j)*(10-j)*(i+j-1) / REAL((i+j)*j*j)
    ipos = ipos + 1
  END DO
  exact_start = -exact_start * (10+i)*(10-i) / REAL((i+1)*i)
END DO

STOP
END PROGRAM hilbert



SUBROUTINE qchol(a, nrows, ier)

! Quadruple-precision Cholesky factorization of a +ve definite symmetric
! matrix, stored by rows up to the diagonal.   The input matrix is
! overwritten by the lower-triangular factor L, where A = LL'.

! Arguments:
! a(:)  INPUT.  Array containing the lower triangle of matrix A.
! nrows INPUT.  The number of rows (and columns) in A.
! ier   OUTPUT. Error indicator
!     = 0 if no error detected
!     = 1 if nrows < 1
!     = 2 if matrix not +ve definite

! Programmer: Alan Miller

! Latest revision - 12 September 1986
! Fortran 90 version - 3 December 1996

USE quadruple_precision
IMPLICIT NONE

TYPE (quad), DIMENSION(:), INTENT(IN OUT) :: a
INTEGER, INTENT(IN)                       :: nrows
INTEGER, INTENT(OUT)                      :: ier

! Local variables

INTEGER     :: i, j, ij, ikstart, jk, ik, k
TYPE (quad) :: total

IF (nrows <= 0) THEN
  ier = 1
  RETURN
END IF

ij = 1
ikstart = 1
DO i = 1, nrows
  jk = 1
  DO j = 1, i
    ik = ikstart
    total = a(ij)

! Form the sum  a(i,j) - l(i,k).l(j,k)

    DO k = 1, j-1
      total = total - a(ik) * a(jk)
      ik = ik + 1
      jk = jk + 1
    END DO

! On the diagonal  l(i,i) = sqrt(total)
! Off diagonal     l(i,j) = total/l(j,j)

    IF (i == j) THEN
      IF (total%hi <= 0._dp) THEN
        ier = 2
        RETURN
      END IF
      a(ij) = SQRT(total)
    ELSE
      a(ij) = total / a(jk)
      jk = jk + 1
    END IF
    ij = ij + 1
  END DO

! ikstart = location of first element in row i.

  ikstart = ij
END DO
ier = 0

RETURN
END SUBROUTINE qchol
