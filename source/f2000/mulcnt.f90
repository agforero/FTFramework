Module m_mulcnt
Use m_uniinv
Private
Integer, Parameter :: kdp = selected_real_kind(15)
public :: mulcnt
private :: kdp
private :: R_mulcnt, I_mulcnt, D_mulcnt
interface mulcnt
  module procedure d_mulcnt, r_mulcnt, i_mulcnt
end interface mulcnt
contains

Subroutine D_mulcnt (XDONT, IMULT)
!   MULCNT = Give for each array value its multiplicity
!            (number of times that it appears in the array)
! __________________________________________________________
!  Michel Olagnon - Mar. 2000
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IMULT
! __________________________________________________________
!
      Integer, Dimension (Size(XDONT)) :: IWRKT
      Integer, Dimension (Size(XDONT)) :: ICNTT
      Integer :: ICRS
! __________________________________________________________
      Call UNIINV (XDONT, IWRKT)
      ICNTT = 0
      Do ICRS = 1, Size(XDONT)
            ICNTT(IWRKT(ICRS)) = ICNTT(IWRKT(ICRS)) + 1
      End Do
      Do ICRS = 1, Size(XDONT)
            IMULT(ICRS) = ICNTT(IWRKT(ICRS))
      End Do

!
End Subroutine D_mulcnt

Subroutine R_mulcnt (XDONT, IMULT)
!   MULCNT = Give for each array value its multiplicity
!            (number of times that it appears in the array)
! __________________________________________________________
!  Michel Olagnon - Mar. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IMULT
! __________________________________________________________
!
      Integer, Dimension (Size(XDONT)) :: IWRKT
      Integer, Dimension (Size(XDONT)) :: ICNTT
      Integer :: ICRS
! __________________________________________________________
      Call UNIINV (XDONT, IWRKT)
      ICNTT = 0
      Do ICRS = 1, Size(XDONT)
            ICNTT(IWRKT(ICRS)) = ICNTT(IWRKT(ICRS)) + 1
      End Do
      Do ICRS = 1, Size(XDONT)
            IMULT(ICRS) = ICNTT(IWRKT(ICRS))
      End Do

!
End Subroutine R_mulcnt
Subroutine I_mulcnt (XDONT, IMULT)
!   MULCNT = Give for each array value its multiplicity
!            (number of times that it appears in the array)
! __________________________________________________________
!  Michel Olagnon - Mar. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IMULT
! __________________________________________________________
!
      Integer, Dimension (Size(XDONT)) :: IWRKT
      Integer, Dimension (Size(XDONT)) :: ICNTT
      Integer :: ICRS
! __________________________________________________________
      Call UNIINV (XDONT, IWRKT)
      ICNTT = 0
      Do ICRS = 1, Size(XDONT)
            ICNTT(IWRKT(ICRS)) = ICNTT(IWRKT(ICRS)) + 1
      End Do
      Do ICRS = 1, Size(XDONT)
            IMULT(ICRS) = ICNTT(IWRKT(ICRS))
      End Do

!
End Subroutine I_mulcnt
end module m_mulcnt
