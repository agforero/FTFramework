PROGRAM symmetric_eigen

! Demonstration program 2 using the 'LONG' package.   Calculation of the
! eigenvalues & vectors of 2 large symmetric matrices using Kaiser's method.
! In this version, quadruple precision is used to orthogonalize the columns
! of the input matrix, but the final calculation of the eigenvalues and
! normalization of the eigenvectors is in double precision.

! Programmer: Alan Miller

! Latest revision - 28 September 1986
! Fortran 90 version - 6 December 1996
! Revised - 19 June 1997

USE quadruple_precision
IMPLICIT NONE

TYPE (quad)       :: a(44,44)
REAL (dp)         :: eigenv(44), trace, sume
INTEGER           :: i, j, nrows, ier
REAL              :: tstart, tend
CHARACTER (LEN=1) :: ch

INTERFACE
  SUBROUTINE qkaiser(a, dim_a, n, eigenv, trace, sume, ier)
    USE quadruple_precision
    IMPLICIT NONE
    INTEGER, INTENT(IN)                         :: dim_a, n
    INTEGER, INTENT(OUT)                        :: ier
    TYPE (quad), DIMENSION(:,:), INTENT(IN OUT) :: a
    REAL (dp), DIMENSION(:), INTENT(OUT)        :: eigenv
    REAL (dp), INTENT(OUT)                      :: trace, sume
  END SUBROUTINE qkaiser
END INTERFACE

! First example 30 x 30

!      ( 10  1  1  1                                  )
!      (  1 10  0  0  1                               )
!      (  1  0 10  0  0  1                            )
!      (  1  0  0  9  0  0  1                         )
!      (     1  0  0  9  0  0  1                      )
! A  = (                 . . . . . . . .              )
!      (                         1  0  0  2  0  0  1  )
!      (                            1  0  0  1  0  0  )
!      (                               1  0  0  1  0  )
!      (                                  1  0  0  1  )

! The published eigenvalues of this matrix are:-
! 11.80931 0238 10.74619 4183 10.17489 2583  9.47341 06874
!  9.21067 86475  8.68784 18564  8.11381 05886  8.03894 11158
!  7.76156 65166  7.01635 02620  7.00395 17985  6.96372 64186
!  6.00126 33582  6.00021 75223  5.99783 78494  4.99983 25857
!  4.99978 24777  4.99968 95662  3.99604 97565  3.99604 82014
!  3.99604 56418  2.96105 89161  2.96105 88842  2.96105 88358
!  1.78932 13529  1.78932 13527  1.78932 13522  0.25380 58171
!  0.25380 58171  0.25380 58170

DO i = 1, 30
  DO j = 1, 30
    a(i,j)%lo = 0._dp
    IF (i == j) THEN
      a(i,i)%hi = (33-i)/3
    ELSE IF (ABS(i-j) == 3) THEN
      a(i,j)%hi = 1._dp
    ELSE
      a(i,j)%hi = 0._dp
    END IF
  END DO
END DO
a(1,2)%hi = 1._dp
a(1,3)%hi = 1._dp
a(2,1)%hi = 1._dp
a(3,1)%hi = 1._dp

nrows = 30
CALL time_now(tstart)
CALL qkaiser(a, 44, nrows, eigenv, trace, sume, ier)
IF (ier == 0) THEN
  WRITE(*, 900) eigenv(1:nrows)
  900 FORMAT(' Eigenvalues:-', 11(/ 4F19.14))
  WRITE(*, *)
  WRITE(*, '(a, 2F19.14)') ' Trace & sum of eigenvalues = ', trace, sume
  WRITE(*, *) 'Trace = sum if all eigenvalues +ve'
ELSE
  WRITE(*, '(a, i4)') ' *** qkaiser error no. ', ier
END IF

CALL time_now(tend)
WRITE(*, '(a, f8.2, a)') ' Time taken = ', tend-tstart, 'secs.'

WRITE(*, *)
WRITE(*, *) 'Press RETURN to continue'
READ(*, '(a)') ch
WRITE(*, '(a)') ' ' // ch

! Second example 44 x 44
!
!      ( 5  2  1  1                               )
!      ( 2  6  3  1  1                            )
!      ( 1  3  6  3  1  1                         )
!      ( 1  1  3  6  3  1  1                      )
!      (    1  1  3  6  3  1  1                   )
! A  = (               . . . . . . . .            )
!      (                         1  1  3  6  3  1 )
!      (                            1  1  3  6  2 )
!      (                               1  1  2  5 )

! This is a particularly difficult matrix as it has 11 eigenvalues
! (nos. 15-25) very close together.   Their published values are:-
! 4.16250 38244 29765   4.14589 80337 50315
! 4.14077 17541 12427   4.12083 46534 01856
! 4.09474 53700 81830   4.07872 56924 26601
! 4.05284 15893 88415   4.03456 80076 76363
! 4.00521 19531 60050   4.00453 18458 00653
! 4.0 (exactly)

DO i = 1, 44
  DO j = 1, 44
    a(i,j)%lo = 0._dp
    IF (i == j) THEN
      a(i,i)%hi = 6._dp
    ELSE IF (ABS(i-j) == 1) THEN
      a(i,j)%hi = 3._dp
    ELSE IF (ABS(i-j) <= 3) THEN
      a(i,j)%hi = 1._dp
    ELSE
      a(i,j)%hi = 0._dp
    END IF
  END DO
END DO
a(1,1)%hi   = 5._dp
a(1,2)%hi   = 2._dp
a(2,1)%hi   = 2._dp
a(44,44)%hi = 5._dp
a(44,43)%hi = 2._dp
a(43,44)%hi = 2._dp

nrows = 44
CALL time_now(tstart)
CALL qkaiser(a, 44, nrows, eigenv, trace, sume, ier)
IF (ier == 0) THEN
  WRITE(*, 900) eigenv(1:nrows)
  WRITE(*, *)
  WRITE(*, '(a, 2F19.14)') ' Trace & sum of eigenvalues = ', trace, sume
  WRITE(*, *) 'Trace = sum if all eigenvalues +ve'
ELSE
  WRITE(*, '(a, i4)') ' *** qkaiser error no. ', ier
END IF

CALL time_now(tend)
WRITE(*, '(a, f8.2, a)') ' Time taken = ', tend-tstart, 'secs.'

STOP


CONTAINS


SUBROUTINE time_now(seconds)
REAL, INTENT(OUT) :: seconds

! Local variables
INTEGER :: t(8)

CALL DATE_AND_TIME(values=t)
seconds = 3600.*t(5) + 60.*t(6) + t(7) + 0.001*t(8)

RETURN
END SUBROUTINE time_now

END PROGRAM symmetric_eigen




SUBROUTINE qkaiser(a, dim_a, n, eigenv, trace, sume, ier)

! Calculate the eigenvalues & vectors of a +ve definite symmetric
! matrix in quadruple precision, using a method due to Kaiser.
! Reference:
! Kaiser, H.F. 'The JK method: a procedure for finding the eigenvalues of
!               a real symmetric matrix', Computer J., 15, 271-273, 1972.
! N.B. Strictly the method finds the squared eigenvalues.   It can be used
! with any symmetric matrix - it does not need to be +ve definite, but it
! does not know what signs to give to the eigenvalues.

! Arguments:
! A(DIM_A,*) INPUT.   Contains the matrix in quadruple-precision.
!            OUTPUT.  The columns contain the eigenvectors.
!                     N.B. The input matrix is overwritten.
! DIM_A      INPUT.   Second dimension of A in the calling program
! N          INPUT.   The order or no. of columns in A.
! EIGENV(N)  OUTPUT.  The eigenvalues (ordered).
! TRACE      OUTPUT.  The trace of A.
! SUME       OUTPUT.  Sum of the eigenvalues.   This should equal the trace
!                     if all the eigenvalues are >= 0.
!                     If not, the difference is equal to twice the sum
!                     of the negative eigenvalues.
! IER        OUTPUT.  Error indicator
!                     = 0 if no errors detected
!                     = 1 if N > DIM_A
!                     = 2 failed to converge in 12 iterations

! Programmer: Alan Miller

! Latest revision - 28 September 1986
! Fortran 90 version - 3 December 1996
! Revised - 5 August 1997

USE quadruple_precision
IMPLICIT NONE

INTEGER, INTENT(IN)                         :: dim_a, n
INTEGER, INTENT(OUT)                        :: ier
TYPE (quad), DIMENSION(:,:), INTENT(IN OUT) :: a
REAL (dp), DIMENSION(:), INTENT(OUT)        :: eigenv
REAL (dp), INTENT(OUT)                      :: trace, sume

! Local variables

TYPE (quad) :: halfp, q, temp, sine, cosine, qone
REAL (dp)   :: small = 1.E-08_dp, very_small = 1.E-25_dp, ss, eps1, eps2,  &
               xj, xk, p, absp, absq, tangent, cotan, zero = 0._dp, one = 1._dp
INTEGER     :: iter, nn, ncount, i, j, k

! Initial checks

IF (n < 1 .OR. n > dim_a) THEN
  ier = 1
  RETURN
END IF
ier = 0
qone%hi = one
qone%lo = zero

! Calculate trace & tolerances.

sume = zero
trace = zero
ss = zero
DO i = 1, n
  trace = trace + a(i,i)%hi
  DO j = 1, n
    ss = ss + a(i,j)%hi**2
  END DO
END DO
eps1 = small * ss/n
eps2 = very_small * ss/n

! Some initialization.

iter = 1
nn = n*(n-1)/2
ncount = nn

!----------------------------------------------------------------------

! Main loop to orthogonalize all pairs of columns starts here.

30 DO j = 1, n-1
  DO k = j+1, n

! Calculate planar rotation in double precision.
! This is fairly fast so it does not matter if we have to repeat
! it in quadruple precision.

    halfp%hi = zero
    q%hi = zero
    DO i = 1, n
      xj = a(i,j)%hi
      xk = a(i,k)%hi
      halfp%hi = halfp%hi + xj*xk
      q%hi = q%hi + (xj + xk) * (xj - xk)
    END DO
    p = SCALE(halfp%hi, 1)
    absp = ABS(p)
    absq = ABS(q%hi)

! If the cross-product in p or the difference of 2-norms in q%hi is
! small, repeat the calculation in quadruple precision.

    IF (absp < eps1 .OR. absq < eps1) THEN
      halfp%hi = zero
      halfp%lo = zero
      q%hi = zero
      q%lo = zero
      DO i = 1, n
        halfp = halfp + a(i,j) * a(i,k)
        q = q + (a(i,j) + a(i,k)) * (a(i,j) - a(i,k))
      END DO
      p = SCALE(halfp%hi, 1)
      absp = ABS(p)
      absq = ABS(q%hi)
    END IF

! Test if columns are almost orthogonal and in the correct order.

    IF (absp < eps2) THEN
      IF (q%hi >= zero) GO TO 80
    END IF

! The trigonometric bit.

    IF (absp <= absq) THEN
      tangent = absp / absq
      temp = exactmul2(tangent, tangent)
      cosine = qone / SQRT(temp + qone)
      temp%hi = tangent
      temp%lo = zero
      sine = temp * cosine
    ELSE
      cotan = absq / absp
      temp = exactmul2(cotan, cotan)
      sine = qone / SQRT(temp + qone)
      temp%hi = cotan
      temp%lo = zero
      cosine = temp * sine
    END IF
    temp = qone + cosine
    temp = SCALE(temp, -1)
    cosine = SQRT(temp)
    temp = SCALE(cosine, 1)
    sine = sine / temp
    IF (q%hi < zero) THEN
      temp = cosine
      cosine = sine
      sine = temp
    END IF
    IF (p < zero) THEN
      sine = -sine
    END IF

! Carry out the rotation in quadruple precision.

    WRITE(*, 900) iter, j, k, sine%hi
    900 FORMAT('+', 'iter:', i3, '   columns:', 2I4, '   sine:', g13.5)
    DO i = 1, n
      q = a(i,j) * cosine + a(i,k) * sine
      a(i,k) = a(i,k) * cosine - a(i,j) * sine
      a(i,j) = q
    END DO
    CYCLE

! No rotation was needed if this point is reached.
! Test for convergence.

    80 ncount = ncount - 1
    IF (ncount <= 0) GO TO 110

! End of loop over pairs of columns.

  END DO
END DO

! Increase count of iterations, then go back for more.

ncount = nn
iter = iter + 1
IF (iter <= 12) GO TO 30
ier = 2
WRITE(*, *)' *** Failed to converge ***'

!----------------------------------------------------------------------

! Converged, or gave up after 12 iterations.
! Eigenvalues squared are the squared lengths of columns.
! The columns of A contain the eigenvectors; scale to length 1.

110 DO j = 1, n
  ss = SUM( a(1:n,j)%hi**2 )
  eigenv(j) = SQRT(ss)
  sume = sume + eigenv(j)
  a(1:n,j)%hi = a(1:n,j)%hi / eigenv(j)
  a(1:n,j)%lo = zero
END DO

RETURN
END SUBROUTINE qkaiser
