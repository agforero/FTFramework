Module m_inspar
Integer, Parameter :: kdp = selected_real_kind(15)
public :: inspar
private :: kdp
private :: R_inspar, I_inspar, D_inspar
interface inspar
  module procedure d_inspar, r_inspar, i_inspar
end interface inspar
contains

Subroutine D_inspar (XDONT, NORD)
!  Sorts partially XDONT, bringing the NORD lowest values at the
!  begining of the array
! __________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It does not use any work array
!  and is faster when NORD is very small (2-5), but worst case
!  behavior can happen fairly probably (initially inverse sorted)
!  In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (InOut) :: XDONT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real (kind=kdp) :: XWRK, XWRK1
!
      Integer :: ICRS, IDCR
!
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XDONT(IDCR)) Exit
            XDONT (IDCR+1) = XDONT (IDCR)
         End Do
         XDONT (IDCR+1) = XWRK
      End Do
!
      XWRK1 = XDONT (NORD)
      Do ICRS = NORD + 1, SIZE (XDONT)
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            XDONT (ICRS) = XWRK1
            Do IDCR = NORD - 1, 1, - 1
               If (XWRK >= XDONT(IDCR)) Exit
               XDONT (IDCR+1) = XDONT (IDCR)
            End Do
            XDONT (IDCR+1) = XWRK
            XWRK1 = XDONT (NORD)
         End If
      End Do
!
!
End Subroutine D_inspar

Subroutine R_inspar (XDONT, NORD)
!  Sorts partially XDONT, bringing the NORD lowest values at the
!  begining of the array
! __________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It does not use any work array
!  and is faster when NORD is very small (2-5), but worst case
!  behavior can happen fairly probably (initially inverse sorted)
!  In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (InOut) :: XDONT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real    :: XWRK, XWRK1
!
      Integer :: ICRS, IDCR
!
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XDONT(IDCR)) Exit
            XDONT (IDCR+1) = XDONT (IDCR)
         End Do
         XDONT (IDCR+1) = XWRK
      End Do
!
      XWRK1 = XDONT (NORD)
      Do ICRS = NORD + 1, SIZE (XDONT)
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            XDONT (ICRS) = XWRK1
            Do IDCR = NORD - 1, 1, - 1
               If (XWRK >= XDONT(IDCR)) Exit
               XDONT (IDCR+1) = XDONT (IDCR)
            End Do
            XDONT (IDCR+1) = XWRK
            XWRK1 = XDONT (NORD)
         End If
      End Do
!
!
End Subroutine R_inspar
Subroutine I_inspar (XDONT, NORD)
!  Sorts partially XDONT, bringing the NORD lowest values at the
!  begining of the array
! __________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It does not use any work array
!  and is faster when NORD is very small (2-5), but worst case
!  behavior can happen fairly probably (initially inverse sorted)
!  In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (InOut)  :: XDONT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Integer :: XWRK, XWRK1
!
      Integer :: ICRS, IDCR
!
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XDONT(IDCR)) Exit
            XDONT (IDCR+1) = XDONT (IDCR)
         End Do
         XDONT (IDCR+1) = XWRK
      End Do
!
      XWRK1 = XDONT (NORD)
      Do ICRS = NORD + 1, SIZE (XDONT)
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            XDONT (ICRS) = XWRK1
            Do IDCR = NORD - 1, 1, - 1
               If (XWRK >= XDONT(IDCR)) Exit
               XDONT (IDCR+1) = XDONT (IDCR)
            End Do
            XDONT (IDCR+1) = XWRK
            XWRK1 = XDONT (NORD)
         End If
      End Do
!
!
End Subroutine I_inspar
end module m_inspar
