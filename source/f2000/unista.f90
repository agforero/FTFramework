Module m_unista
Use m_uniinv
Private
Integer, Parameter :: kdp = selected_real_kind(15)
public :: unista
private :: kdp
private :: R_unista, I_unista, D_unista
interface unista
  module procedure d_unista, r_unista, i_unista
end interface unista
contains

Subroutine D_unista (XDONT, NUNI)
!   UNISTA = (Stable unique) Removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (InOut) :: XDONT
      Integer, Intent (Out) :: NUNI
! __________________________________________________________
!
      Integer, Dimension (Size(XDONT)) :: IWRKT
      Logical, Dimension (Size(XDONT)) :: IFMPTYT
      Integer :: ICRS
! __________________________________________________________
      Call UNIINV (XDONT, IWRKT)
      IFMPTYT = .True.
      NUNI = 0
      Do ICRS = 1, Size(XDONT)
         If (IFMPTYT(IWRKT(ICRS))) Then
            IFMPTYT(IWRKT(ICRS)) = .False.
            NUNI = NUNI + 1
            XDONT (NUNI) = XDONT (ICRS)
         End If
      End Do
      Return
!
End Subroutine D_unista

Subroutine R_unista (XDONT, NUNI)
!   UNISTA = (Stable unique) Removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (InOut) :: XDONT
      Integer, Intent (Out) :: NUNI
! __________________________________________________________
!
      Integer, Dimension (Size(XDONT)) :: IWRKT
      Logical, Dimension (Size(XDONT)) :: IFMPTYT
      Integer :: ICRS
! __________________________________________________________
      Call UNIINV (XDONT, IWRKT)
      IFMPTYT = .True.
      NUNI = 0
      Do ICRS = 1, Size(XDONT)
         If (IFMPTYT(IWRKT(ICRS))) Then
            IFMPTYT(IWRKT(ICRS)) = .False.
            NUNI = NUNI + 1
            XDONT (NUNI) = XDONT (ICRS)
         End If
      End Do
      Return
!
End Subroutine R_unista
Subroutine I_unista (XDONT, NUNI)
!   UNISTA = (Stable unique) Removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (InOut)  :: XDONT
      Integer, Intent (Out) :: NUNI
! __________________________________________________________
!
      Integer, Dimension (Size(XDONT)) :: IWRKT
      Logical, Dimension (Size(XDONT)) :: IFMPTYT
      Integer :: ICRS
! __________________________________________________________
      Call UNIINV (XDONT, IWRKT)
      IFMPTYT = .True.
      NUNI = 0
      Do ICRS = 1, Size(XDONT)
         If (IFMPTYT(IWRKT(ICRS))) Then
            IFMPTYT(IWRKT(ICRS)) = .False.
            NUNI = NUNI + 1
            XDONT (NUNI) = XDONT (ICRS)
         End If
      End Do
      Return
!
End Subroutine I_unista
end module m_unista
