Module m_inssor
Integer, Parameter :: kdp = selected_real_kind(15)
public :: inssor
private :: kdp
private :: R_inssor, I_inssor, D_inssor
interface inssor
  module procedure d_inssor, r_inssor, i_inssor
end interface inssor
contains

Subroutine D_inssor (XDONT)
!  Sorts XDONT into increasing order (Insertion sort)
! __________________________________________________________
!  This subroutine uses insertion sort. It does not use any
!  work array and is faster when XDONT is of very small size
!  (< 20), or already almost sorted, but worst case behavior
!  can happen fairly probably (initially inverse sorted).
!  In many cases, the quicksort or merge sort method is faster.
!  Michel Olagnon - Apr. 2000
! __________________________________________________________
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (InOut) :: XDONT
! __________________________________________________________
      Real (Kind=kdp) :: XWRK, XMIN
!
! __________________________________________________________
!
      Integer :: ICRS, IDCR, NDON
!
      NDON = Size (XDONT)
!
! We first bring the minimum to the first location in the array.
! That way, we will have a "guard", and when looking for the
! right place to insert a value, no loop test is necessary.
!
      If (XDONT (1) < XDONT (NDON)) Then
          XMIN = XDONT (1)
      Else
          XMIN = XDONT (NDON)
          XDONT (NDON) = XDONT (1)
      Endif
      Do IDCR = NDON-1, 2, -1
         XWRK = XDONT(IDCR)
         IF (XWRK < XMIN) Then
            XDONT (IDCR) = XMIN
            XMIN = XWRK
         End If
      End Do
      XDONT (1) = XMIN
!
! The first value is now the minimum
! Loop over the array, and when a value is smaller than
! the previous one, loop down to insert it at its right place.
!
      Do ICRS = 3, NDON
         XWRK = XDONT (ICRS)
         IDCR = ICRS - 1
         If (XWRK < XDONT(IDCR)) Then
            XDONT (ICRS) = XDONT (IDCR)
            IDCR = IDCR - 1
            Do
               If (XWRK >= XDONT(IDCR)) Exit
               XDONT (IDCR+1) = XDONT (IDCR)
               IDCR = IDCR - 1
            End Do
            XDONT (IDCR+1) = XWRK
         End If
      End Do
!
      Return
!
End Subroutine D_inssor

Subroutine R_inssor (XDONT)
!  Sorts XDONT into increasing order (Insertion sort)
! __________________________________________________________
!  This subroutine uses insertion sort. It does not use any
!  work array and is faster when XDONT is of very small size
!  (< 20), or already almost sorted, but worst case behavior
!  can happen fairly probably (initially inverse sorted).
!  In many cases, the quicksort or merge sort method is faster.
!  Michel Olagnon - Apr. 2000
! __________________________________________________________
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (InOut) :: XDONT
! __________________________________________________________
      Real :: XWRK, XMIN
!
! __________________________________________________________
!
      Integer :: ICRS, IDCR, NDON
!
      NDON = Size (XDONT)
!
! We first bring the minimum to the first location in the array.
! That way, we will have a "guard", and when looking for the
! right place to insert a value, no loop test is necessary.
!
      If (XDONT (1) < XDONT (NDON)) Then
          XMIN = XDONT (1)
      Else
          XMIN = XDONT (NDON)
          XDONT (NDON) = XDONT (1)
      Endif
      Do IDCR = NDON-1, 2, -1
         XWRK = XDONT(IDCR)
         IF (XWRK < XMIN) Then
            XDONT (IDCR) = XMIN
            XMIN = XWRK
         End If
      End Do
      XDONT (1) = XMIN
!
! The first value is now the minimum
! Loop over the array, and when a value is smaller than
! the previous one, loop down to insert it at its right place.
!
      Do ICRS = 3, NDON
         XWRK = XDONT (ICRS)
         IDCR = ICRS - 1
         If (XWRK < XDONT(IDCR)) Then
            XDONT (ICRS) = XDONT (IDCR)
            IDCR = IDCR - 1
            Do
               If (XWRK >= XDONT(IDCR)) Exit
               XDONT (IDCR+1) = XDONT (IDCR)
               IDCR = IDCR - 1
            End Do
            XDONT (IDCR+1) = XWRK
         End If
      End Do
!
      Return
!
End Subroutine R_inssor
Subroutine I_inssor (XDONT)
!  Sorts XDONT into increasing order (Insertion sort)
! __________________________________________________________
!  This subroutine uses insertion sort. It does not use any
!  work array and is faster when XDONT is of very small size
!  (< 20), or already almost sorted, but worst case behavior
!  can happen fairly probably (initially inverse sorted).
!  In many cases, the quicksort or merge sort method is faster.
!  Michel Olagnon - Apr. 2000
! __________________________________________________________
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (InOut)  :: XDONT
! __________________________________________________________
      Integer :: XWRK, XMIN
!
! __________________________________________________________
!
      Integer :: ICRS, IDCR, NDON
!
      NDON = Size (XDONT)
!
! We first bring the minimum to the first location in the array.
! That way, we will have a "guard", and when looking for the
! right place to insert a value, no loop test is necessary.
!
      If (XDONT (1) < XDONT (NDON)) Then
          XMIN = XDONT (1)
      Else
          XMIN = XDONT (NDON)
          XDONT (NDON) = XDONT (1)
      Endif
      Do IDCR = NDON-1, 2, -1
         XWRK = XDONT(IDCR)
         IF (XWRK < XMIN) Then
            XDONT (IDCR) = XMIN
            XMIN = XWRK
         End If
      End Do
      XDONT (1) = XMIN
!
! The first value is now the minimum
! Loop over the array, and when a value is smaller than
! the previous one, loop down to insert it at its right place.
!
      Do ICRS = 3, NDON
         XWRK = XDONT (ICRS)
         IDCR = ICRS - 1
         If (XWRK < XDONT(IDCR)) Then
            XDONT (ICRS) = XDONT (IDCR)
            IDCR = IDCR - 1
            Do
               If (XWRK >= XDONT(IDCR)) Exit
               XDONT (IDCR+1) = XDONT (IDCR)
               IDCR = IDCR - 1
            End Do
            XDONT (IDCR+1) = XWRK
         End If
      End Do
!
      Return
!
End Subroutine I_inssor
end module m_inssor
