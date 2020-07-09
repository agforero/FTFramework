Module m_refpar
Integer, Parameter :: kdp = selected_real_kind(15)
public :: refpar
private :: kdp
private :: R_refpar, I_refpar, D_refpar
interface refpar
  module procedure d_refpar, r_refpar, i_refpar
end interface refpar
contains

Subroutine D_refpar (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm. It uses
!  a temporary array, where it stores the partially ranked indices
!  of the values. It iterates until it can bring the number of
!  values lower than the pivot to exactly NORD, and then uses an
!  insertion sort to rank this set, since it is supposedly small.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real (kind=kdp) :: XPIV, XWRK
! __________________________________________________________
!
      Integer, Dimension (SIZE(XDONT)) :: IWRKT
      Integer :: NDON, ICRS, IDEB, IDCR, IFIN, IMIL, IWRK
!
      NDON = SIZE (XDONT)
!
      Do ICRS = 1, NDON
         IWRKT (ICRS) = ICRS
      End Do
      IDEB = 1
      IFIN = NDON
      Do
         If (IDEB >= IFIN) Exit
         IMIL = (IDEB+IFIN) / 2
!
!  One chooses a pivot, median of 1st, last, and middle values
!
         If (XDONT(IWRKT(IMIL)) < XDONT(IWRKT(IDEB))) Then
            IWRK = IWRKT (IDEB)
            IWRKT (IDEB) = IWRKT (IMIL)
            IWRKT (IMIL) = IWRK
         End If
         If (XDONT(IWRKT(IMIL)) > XDONT(IWRKT(IFIN))) Then
            IWRK = IWRKT (IFIN)
            IWRKT (IFIN) = IWRKT (IMIL)
            IWRKT (IMIL) = IWRK
            If (XDONT(IWRKT(IMIL)) < XDONT(IWRKT(IDEB))) Then
               IWRK = IWRKT (IDEB)
               IWRKT (IDEB) = IWRKT (IMIL)
               IWRKT (IMIL) = IWRK
            End If
         End If
         If ((IFIN-IDEB) < 3) Exit
         XPIV = XDONT (IWRKT(IMIL))
!
!  One exchanges values to put those > pivot in the end and
!  those <= pivot at the beginning
!
         ICRS = IDEB
         IDCR = IFIN
         ECH2: Do
            Do
               ICRS = ICRS + 1
               If (ICRS >= IDCR) Then
!
!  the first  >  pivot is IWRKT(IDCR)
!  the last   <= pivot is IWRKT(ICRS-1)
!  Note: If one arrives here on the first iteration, then
!        the pivot is the maximum of the set, the last value is equal
!        to it, and one can reduce by one the size of the set to process,
!        as if XDONT (IWRKT(IFIN)) > XPIV
!
                  Exit ECH2
!
               End If
               If (XDONT(IWRKT(ICRS)) > XPIV) Exit
            End Do
            Do
               If (XDONT(IWRKT(IDCR)) <= XPIV) Exit
               IDCR = IDCR - 1
               If (ICRS >= IDCR) Then
!
!  The last value < pivot is always IWRKT(ICRS-1)
!
                  Exit ECH2
               End If
            End Do
!
            IWRK = IWRKT (IDCR)
            IWRKT (IDCR) = IWRKT (ICRS)
            IWRKT (ICRS) = IWRK
         End Do ECH2
!
!  One restricts further processing to find the fractile value
!
         If (ICRS <= NORD) IDEB = ICRS
         If (ICRS > NORD) IFIN = ICRS - 1
      End Do
!
!  Now, we only need to complete ranking of the 1:NORD set
!  Assuming NORD is small, we use a simple insertion sort
!
      Do ICRS = 2, NORD
         IWRK = IWRKT (ICRS)
         XWRK = XDONT (IWRK)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK <= XDONT(IWRKT(IDCR))) Then
               IWRKT (IDCR+1) = IWRKT (IDCR)
            Else
               Exit
            End If
         End Do
         IWRKT (IDCR+1) = IWRK
      End Do
      IRNGT (1:NORD) = IWRKT (1:NORD)
      Return
!
End Subroutine D_refpar

Subroutine R_refpar (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm. It uses
!  a temporary array, where it stores the partially ranked indices
!  of the values. It iterates until it can bring the number of
!  values lower than the pivot to exactly NORD, and then uses an
!  insertion sort to rank this set, since it is supposedly small.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real :: XPIV, XWRK
! __________________________________________________________
!
      Integer, Dimension (SIZE(XDONT)) :: IWRKT
      Integer :: NDON, ICRS, IDEB, IDCR, IFIN, IMIL, IWRK
!
      NDON = SIZE (XDONT)
!
      Do ICRS = 1, NDON
         IWRKT (ICRS) = ICRS
      End Do
      IDEB = 1
      IFIN = NDON
      Do
         If (IDEB >= IFIN) Exit
         IMIL = (IDEB+IFIN) / 2
!
!  One chooses a pivot, median of 1st, last, and middle values
!
         If (XDONT(IWRKT(IMIL)) < XDONT(IWRKT(IDEB))) Then
            IWRK = IWRKT (IDEB)
            IWRKT (IDEB) = IWRKT (IMIL)
            IWRKT (IMIL) = IWRK
         End If
         If (XDONT(IWRKT(IMIL)) > XDONT(IWRKT(IFIN))) Then
            IWRK = IWRKT (IFIN)
            IWRKT (IFIN) = IWRKT (IMIL)
            IWRKT (IMIL) = IWRK
            If (XDONT(IWRKT(IMIL)) < XDONT(IWRKT(IDEB))) Then
               IWRK = IWRKT (IDEB)
               IWRKT (IDEB) = IWRKT (IMIL)
               IWRKT (IMIL) = IWRK
            End If
         End If
         If ((IFIN-IDEB) < 3) Exit
         XPIV = XDONT (IWRKT(IMIL))
!
!  One exchanges values to put those > pivot in the end and
!  those <= pivot at the beginning
!
         ICRS = IDEB
         IDCR = IFIN
         ECH2: Do
            Do
               ICRS = ICRS + 1
               If (ICRS >= IDCR) Then
!
!  the first  >  pivot is IWRKT(IDCR)
!  the last   <= pivot is IWRKT(ICRS-1)
!  Note: If one arrives here on the first iteration, then
!        the pivot is the maximum of the set, the last value is equal
!        to it, and one can reduce by one the size of the set to process,
!        as if XDONT (IWRKT(IFIN)) > XPIV
!
                  Exit ECH2
!
               End If
               If (XDONT(IWRKT(ICRS)) > XPIV) Exit
            End Do
            Do
               If (XDONT(IWRKT(IDCR)) <= XPIV) Exit
               IDCR = IDCR - 1
               If (ICRS >= IDCR) Then
!
!  The last value < pivot is always IWRKT(ICRS-1)
!
                  Exit ECH2
               End If
            End Do
!
            IWRK = IWRKT (IDCR)
            IWRKT (IDCR) = IWRKT (ICRS)
            IWRKT (ICRS) = IWRK
         End Do ECH2
!
!  One restricts further processing to find the fractile value
!
         If (ICRS <= NORD) IDEB = ICRS
         If (ICRS > NORD) IFIN = ICRS - 1
      End Do
!
!  Now, we only need to complete ranking of the 1:NORD set
!  Assuming NORD is small, we use a simple insertion sort
!
      Do ICRS = 2, NORD
         IWRK = IWRKT (ICRS)
         XWRK = XDONT (IWRK)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK <= XDONT(IWRKT(IDCR))) Then
               IWRKT (IDCR+1) = IWRKT (IDCR)
            Else
               Exit
            End If
         End Do
         IWRKT (IDCR+1) = IWRK
      End Do
      IRNGT (1:NORD) = IWRKT (1:NORD)
      Return
!
End Subroutine R_refpar
Subroutine I_refpar (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm. It uses
!  a temporary array, where it stores the partially ranked indices
!  of the values. It iterates until it can bring the number of
!  values lower than the pivot to exactly NORD, and then uses an
!  insertion sort to rank this set, since it is supposedly small.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Integer :: XPIV, XWRK
!
      Integer, Dimension (SIZE(XDONT)) :: IWRKT
      Integer :: NDON, ICRS, IDEB, IDCR, IFIN, IMIL, IWRK
!
      NDON = SIZE (XDONT)
!
      Do ICRS = 1, NDON
         IWRKT (ICRS) = ICRS
      End Do
      IDEB = 1
      IFIN = NDON
      Do
         If (IDEB >= IFIN) Exit
         IMIL = (IDEB+IFIN) / 2
!
!  One chooses a pivot, median of 1st, last, and middle values
!
         If (XDONT(IWRKT(IMIL)) < XDONT(IWRKT(IDEB))) Then
            IWRK = IWRKT (IDEB)
            IWRKT (IDEB) = IWRKT (IMIL)
            IWRKT (IMIL) = IWRK
         End If
         If (XDONT(IWRKT(IMIL)) > XDONT(IWRKT(IFIN))) Then
            IWRK = IWRKT (IFIN)
            IWRKT (IFIN) = IWRKT (IMIL)
            IWRKT (IMIL) = IWRK
            If (XDONT(IWRKT(IMIL)) < XDONT(IWRKT(IDEB))) Then
               IWRK = IWRKT (IDEB)
               IWRKT (IDEB) = IWRKT (IMIL)
               IWRKT (IMIL) = IWRK
            End If
         End If
         If ((IFIN-IDEB) < 3) Exit
         XPIV = XDONT (IWRKT(IMIL))
!
!  One exchanges values to put those > pivot in the end and
!  those <= pivot at the beginning
!
         ICRS = IDEB
         IDCR = IFIN
         ECH2: Do
            Do
               ICRS = ICRS + 1
               If (ICRS >= IDCR) Then
!
!  the first  >  pivot is IWRKT(IDCR)
!  the last   <= pivot is IWRKT(ICRS-1)
!  Note: If one arrives here on the first iteration, then
!        the pivot is the maximum of the set, the last value is equal
!        to it, and one can reduce by one the size of the set to process,
!        as if XDONT (IWRKT(IFIN)) > XPIV
!
                  Exit ECH2
!
               End If
               If (XDONT(IWRKT(ICRS)) > XPIV) Exit
            End Do
            Do
               If (XDONT(IWRKT(IDCR)) <= XPIV) Exit
               IDCR = IDCR - 1
               If (ICRS >= IDCR) Then
!
!  The last value < pivot is always IWRKT(ICRS-1)
!
                  Exit ECH2
               End If
            End Do
!
            IWRK = IWRKT (IDCR)
            IWRKT (IDCR) = IWRKT (ICRS)
            IWRKT (ICRS) = IWRK
         End Do ECH2
!
!  One restricts further processing to find the fractile value
!
         If (ICRS <= NORD) IDEB = ICRS
         If (ICRS > NORD) IFIN = ICRS - 1
      End Do
!
!  Now, we only need to complete ranking of the 1:NORD set
!  Assuming NORD is small, we use a simple insertion sort
!
      Do ICRS = 2, NORD
         IWRK = IWRKT (ICRS)
         XWRK = XDONT (IWRK)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK <= XDONT(IWRKT(IDCR))) Then
               IWRKT (IDCR+1) = IWRKT (IDCR)
            Else
               Exit
            End If
         End Do
         IWRKT (IDCR+1) = IWRK
      End Do
      IRNGT (1:NORD) = IWRKT (1:NORD)
      Return
!
End Subroutine I_refpar
end module m_refpar
