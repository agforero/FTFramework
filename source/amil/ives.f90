SUBROUTINE ives(n)

! Generate all permutations of the N elements in array IX.
! This is algorithm 4a (Alternate Ives) from page 149 of:-
!   Sedgewick, R. Permutation generation methods, ACM Computing Surveys,
!                 9, 137-164, 1977.

! N.B. The user's subroutine for doing something with each new permutation
!      is called from this subroutine.

IMPLICIT NONE
INTEGER, INTENT(IN) :: n

! Local variables
INTEGER  :: ix(n), INDEX(15), i, itemp, j, np1mi

DO i = 1,n
  INDEX(i) = 1
END DO
i = 1

!     Insert initial call to user's process here.

DO
  j = INDEX(i)
  np1mi = n+1-i
  IF(j >= np1mi) GO TO 50
  IF(i == 2*(i/2)) GO TO 30
  itemp = ix(j)
  ix(j) = ix(j+1)
  ix(j+1) = itemp
  GO TO 40

  30 itemp = ix(i)
  ix(i) = ix(np1mi)
  ix(np1mi) = itemp

  40 INDEX(i) = j+1
  i = 1

  !     Insert main call to user's process here.

  CYCLE
  50 INDEX(i) = 1
  i = i+1
  IF(i > n) EXIT
END DO

RETURN
END SUBROUTINE ives
