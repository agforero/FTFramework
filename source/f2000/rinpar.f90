Module m_rinpar
Integer, Parameter :: kdp = selected_real_kind(15)
public :: rinpar
private :: kdp
private :: R_rinpar, I_rinpar, D_rinpar
interface rinpar
  module procedure d_rinpar, r_rinpar, i_rinpar
end interface rinpar
contains

Subroutine D_rinpar (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD = size (IRNGT)
! __________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It does not use any work array
!  and is faster when NORD is very small (2-5), but worst case
!  behavior can happen fairly probably (initially inverse sorted)
!  In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real (kind=kdp) :: XWRK, XWRK1
!
      Integer :: ICRS, IDCR
!
      IRNGT (1) = 1
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XDONT(IRNGT(IDCR))) Exit
            IRNGT (IDCR+1) = IRNGT (IDCR)
         End Do
         IRNGT (IDCR+1) = ICRS
      End Do
!
      XWRK1 = XDONT (IRNGT(NORD))
      Do ICRS = NORD + 1, SIZE (XDONT)
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            Do IDCR = NORD - 1, 1, - 1
               If (XWRK >= XDONT(IRNGT(IDCR))) Exit
               IRNGT (IDCR+1) = IRNGT (IDCR)
            End Do
            IRNGT (IDCR+1) = ICRS
            XWRK1 = XDONT (IRNGT(NORD))
         End If
      End Do
!
!
End Subroutine D_rinpar

Subroutine R_rinpar (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD = size (IRNGT)
! __________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It does not use any work array
!  and is faster when NORD is very small (2-5), but worst case
!  behavior can happen fairly probably (initially inverse sorted)
!  In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real    :: XWRK, XWRK1
!
      Integer :: ICRS, IDCR
!
      IRNGT (1) = 1
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XDONT(IRNGT(IDCR))) Exit
            IRNGT (IDCR+1) = IRNGT (IDCR)
         End Do
         IRNGT (IDCR+1) = ICRS
      End Do
!
      XWRK1 = XDONT (IRNGT(NORD))
      Do ICRS = NORD + 1, SIZE (XDONT)
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            Do IDCR = NORD - 1, 1, - 1
               If (XWRK >= XDONT(IRNGT(IDCR))) Exit
               IRNGT (IDCR+1) = IRNGT (IDCR)
            End Do
            IRNGT (IDCR+1) = ICRS
            XWRK1 = XDONT (IRNGT(NORD))
         End If
      End Do
!
!
End Subroutine R_rinpar
Subroutine I_rinpar (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD = size (IRNGT)
! __________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It does not use any work array
!  and is faster when NORD is very small (2-5), but worst case
!  behavior can happen fairly probably (initially inverse sorted)
!  In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Integer :: XWRK, XWRK1
!
      Integer :: ICRS, IDCR
!
      IRNGT (1) = 1
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XDONT(IRNGT(IDCR))) Exit
            IRNGT (IDCR+1) = IRNGT (IDCR)
         End Do
         IRNGT (IDCR+1) = ICRS
      End Do
!
      XWRK1 = XDONT (IRNGT(NORD))
      Do ICRS = NORD + 1, SIZE (XDONT)
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            Do IDCR = NORD - 1, 1, - 1
               If (XWRK >= XDONT(IRNGT(IDCR))) Exit
               IRNGT (IDCR+1) = IRNGT (IDCR)
            End Do
            IRNGT (IDCR+1) = ICRS
            XWRK1 = XDONT (IRNGT(NORD))
         End If
      End Do
!
!
End Subroutine I_rinpar
end module m_rinpar
