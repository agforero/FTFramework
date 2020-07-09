Module m_fndnth
Integer, Parameter :: kdp = selected_real_kind(15)
public :: fndnth
private :: kdp
private :: R_fndnth, I_fndnth, D_fndnth
interface fndnth
  module procedure d_fndnth, r_fndnth, i_fndnth
end interface fndnth
contains

Function D_fndnth (XDONT, NORD) Result (FNDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! ______________________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It is faster when NORD is very small (2-5),
!  and it requires only a workarray of size NORD and type of XDONT,
!  but worst case behavior can happen fairly probably (initially inverse
!  sorted). In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! __________________________________________________________
      Real (Kind=kdp), Dimension (:), Intent (In) :: XDONT
      Real (Kind=kdp) :: FNDNTH
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real (Kind=kdp), Dimension (NORD) :: XWRKT
      Real (Kind=kdp) :: XWRK, XWRK1
!
!
      Integer :: ICRS, IDCR, ILOW, NDON
!
      XWRKT (1) = XDONT (1)
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XWRKT(IDCR)) Exit
            XWRKT (IDCR+1) = XWRKT (IDCR)
         End Do
         XWRKT (IDCR+1) = XWRK
      End Do
!
      NDON = SIZE (XDONT)
      XWRK1 = XWRKT (NORD)
      ILOW = 2*NORD - NDON
      Do ICRS = NORD + 1, NDON
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            Do IDCR = NORD - 1, MAX (1, ILOW) , - 1
               If (XWRK >= XWRKT(IDCR)) Exit
               XWRKT (IDCR+1) = XWRKT (IDCR)
            End Do
            XWRKT (IDCR+1) = XWRK
            XWRK1 = XWRKT(NORD)
         End If
         ILOW = ILOW + 1
      End Do
      FNDNTH = XWRK1

!
End Function D_fndnth

Function R_fndnth (XDONT, NORD) Result (FNDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! ______________________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It is faster when NORD is very small (2-5),
!  and it requires only a workarray of size NORD and type of XDONT,
!  but worst case behavior can happen fairly probably (initially inverse
!  sorted). In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Real :: FNDNTH
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real, Dimension (NORD) :: XWRKT
      Real :: XWRK, XWRK1
!
!
      Integer :: ICRS, IDCR, ILOW, NDON
!
      XWRKT (1) = XDONT (1)
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XWRKT(IDCR)) Exit
            XWRKT (IDCR+1) = XWRKT (IDCR)
         End Do
         XWRKT (IDCR+1) = XWRK
      End Do
!
      NDON = SIZE (XDONT)
      XWRK1 = XWRKT (NORD)
      ILOW = 2*NORD - NDON
      Do ICRS = NORD + 1, NDON
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            Do IDCR = NORD - 1, MAX (1, ILOW) , - 1
               If (XWRK >= XWRKT(IDCR)) Exit
               XWRKT (IDCR+1) = XWRKT (IDCR)
            End Do
            XWRKT (IDCR+1) = XWRK
            XWRK1 = XWRKT(NORD)
         End If
         ILOW = ILOW + 1
      End Do
      FNDNTH = XWRK1

!
End Function R_fndnth
Function I_fndnth (XDONT, NORD) Result (FNDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! ______________________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It is faster when NORD is very small (2-5),
!  and it requires only a workarray of size NORD and type of XDONT,
!  but worst case behavior can happen fairly probably (initially inverse
!  sorted). In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In) :: XDONT
      Integer :: fndnth
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Integer, Dimension (NORD) :: XWRKT
      Integer :: XWRK, XWRK1
!
!
      Integer :: ICRS, IDCR, ILOW, NDON
!
      XWRKT (1) = XDONT (1)
      Do ICRS = 2, NORD
         XWRK = XDONT (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XWRKT(IDCR)) Exit
            XWRKT (IDCR+1) = XWRKT (IDCR)
         End Do
         XWRKT (IDCR+1) = XWRK
      End Do
!
      NDON = SIZE (XDONT)
      XWRK1 = XWRKT (NORD)
      ILOW = 2*NORD - NDON
      Do ICRS = NORD + 1, NDON
         If (XDONT(ICRS) < XWRK1) Then
            XWRK = XDONT (ICRS)
            Do IDCR = NORD - 1, MAX (1, ILOW) , - 1
               If (XWRK >= XWRKT(IDCR)) Exit
               XWRKT (IDCR+1) = XWRKT (IDCR)
            End Do
            XWRKT (IDCR+1) = XWRK
            XWRK1 = XWRKT(NORD)
         End If
         ILOW = ILOW + 1
      End Do
      FNDNTH = XWRK1

!
End Function I_fndnth
end module m_fndnth
