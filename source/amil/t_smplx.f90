PROGRAM test_smplx
!     Test the F90 version of the NSWC routine SMPLX.

USE constants_NSWC
IMPLICIT NONE
INTEGER                :: n, pos, i, constr_type(9), m, ind, ibasis(9), iter, &
                          mxiter = 25, numle, numge, j
CHARACTER (LEN=79)     :: text, in_text
REAL (dp)              :: a(9,9), b(9), c(9), x(9), z, temp(9), rerr
REAL (dp), ALLOCATABLE :: bi(:,:)

INTERFACE
  SUBROUTINE smplx (a, b0, c, ka, m, n0, ind, ibasis, x, z, iter, mxiter,   &
                    numle, numge, bi, rerr)
    USE constants_NSWC
    IMPLICIT NONE
    INTEGER, INTENT(IN)                    :: ka, m, n0, mxiter, numle, numge
    INTEGER, INTENT(IN OUT)                :: ind
    INTEGER, DIMENSION(:), INTENT(IN OUT)  :: ibasis
    INTEGER, INTENT(OUT)                   :: iter
    REAL (dp), DIMENSION(:,:), INTENT(IN)  :: a
    REAL (dp), DIMENSION(:), INTENT(IN)    :: b0, c
    REAL (dp), INTENT(OUT)                 :: z, rerr
    REAL (dp), DIMENSION(:), INTENT(OUT)   :: x
    REAL (dp), DIMENSION(:,:), INTENT(OUT) :: bi
  END SUBROUTINE smplx
END INTERFACE

WRITE(*, *) 'How many variables? (<=9): '
READ(*, *) n
IF (n > 9) STOP 'Too many!'

WRITE(*, *) 'Enter constraint coefficients incl. constr. type (<=, >= or =)'
WRITE(*, *) 'Press ENTER again after the last constraint'
text = ' '
DO i = 1, n
  pos = 8*i-6
  WRITE(text(pos:pos+3), '("X(", i1, ")")') i
END DO
text(pos+4:pos+11) = 'Con.type'
text(pos+13:pos+15) = 'RHS'
WRITE(*, '(a)') text
WRITE(*, *)
numle = 0
numge = 0
DO i = 1, 9
  READ(*, '(a)') in_text
  IF (LEN_TRIM(in_text) == 0) EXIT
  READ(in_text, *) a(i,1:n)
  pos = INDEX(in_text, '=')
  IF (pos == 0) THEN
    STOP 'You must enter the constraint type (<=, >= or =)'
  END IF
  IF (in_text(pos-1:pos-1) == '<') THEN
    constr_type(i) = 1
    numle = numle + 1
  ELSE IF (in_text(pos-1:pos-1) == '>') THEN
    constr_type(i) = 2
    numge = numge + 1
  ELSE
    constr_type(i) = 3
  END IF
  READ(in_text(pos+1:), *) b(i)
END DO

m = i - 1
ALLOCATE( bi(m,m) )

!     Re-order the constraints if necessary
DO i = 1, numle
  IF (constr_type(i) /= 1) THEN
    DO j = numle+1, m
      IF (constr_type(j) == 1) THEN
        temp(1:n) = a(i, 1:n)
        a(i, 1:n) = a(j, 1:n)
        a(j, 1:n) = temp(1:n)
        z = b(i)
        b(i) = b(j)
        b(j) = z
        constr_type(j) = constr_type(i)
        constr_type(i) = 1
        EXIT
      END IF
    END DO
  END IF
END DO
DO i = numle+1, numle+numge
  IF (constr_type(i) /= 2) THEN
    DO j = numle+numge+1, m
      IF (constr_type(j) == 2) THEN
        temp(1:n) = a(i, 1:n)
        a(i, 1:n) = a(j, 1:n)
        a(j, 1:n) = temp(1:n)
        z = b(i)
        b(i) = b(j)
        b(j) = z
        constr_type(j) = constr_type(i)
        constr_type(i) = 2
        EXIT
      END IF
    END DO
  END IF
END DO

!     Enter the coefficients in the objective function.
WRITE(*, *) 'Enter the coefficients (costs) in the objective to be maximized'
pos = 8*n-2
text(pos:) = ' '
WRITE(*, '(a)') text
WRITE(*, *)
READ(*, '(a)') in_text
READ(in_text, *) c(1:n)

CALL smplx (a, b, c, 9, m, n, ind, ibasis, x, z, iter, mxiter,   &
            numle, numge, bi, rerr)

WRITE(*, '(a, i2, a, i3)') ' IND = ', ind, '  No. of iterations = ', iter
IF (ind == 0 .OR. i== 6) THEN
  WRITE(*, '(a, 9f7.2)') ' Solution: ', x(1:n)
  WRITE(*, '(a, f8.2)') ' Value of objective = ', z
  WRITE(*, *) 'The inverse matrix:'
  DO i = 1, m
    WRITE(*, '(9f8.3)') bi(i,1:m)
  END DO
END IF

WRITE(*, '(a, g12.4)') ' Relative error = ', rerr

STOP
END PROGRAM test_smplx
