Module m_rapknr
Integer, Parameter :: kdp = selected_real_kind(15)
public :: rapknr
private :: kdp
private :: R_rapknr, I_rapknr, D_rapknr
interface rapknr
  module procedure d_rapknr, r_rapknr, i_rapknr
end interface rapknr
contains

Subroutine D_rapknr (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD, in decreasing order.
!  rapknr = (rnkpar backwards)
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values larger than the pivot
!  (IHIGT), and the indices of values smaller than the pivot
!  that we might still need later on (ILOWT). It iterates
!  until it can bring the number of values in IHIGT to
!  exactly NORD, and then uses an insertion sort to rank
!  this set, since it is supposedly small.
!  Michel Olagnon - Feb. 2011
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real (kind=kdp) :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      Integer, Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      Integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      Integer :: IDEB, JDEB, IMIL, IFIN, NWRK, ICRS, IDCR, ILOW
      Integer :: JLM2, JLM1, JHM2, JHM1
!
      NDON = SIZE (XDONT)
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      If (NDON < 2) Then
         If (NORD >= 1) IRNGT (1) = 1
         Return
      End If
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of high values.
!
      If (XDONT(2) < XDONT(1)) Then
         ILOWT (1) = 2
         IHIGT (1) = 1
      Else
         ILOWT (1) = 1
         IHIGT (1) = 2
      End If
!
      If (NDON < 3) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         Return
      End If
! ---
      If (XDONT(3) > XDONT(ILOWT(1))) Then
         ILOWT (2) = ILOWT (1)
         If (XDONT(3) > XDONT(IHIGT(1))) Then
            ILOWT (1) = IHIGT (1)
            IHIGT (1) = 3
         Else
            ILOWT (1) = 3
         End If
      Else
         ILOWT (2) = 3
      End If
! ---
      If (NDON < 4) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         If (NORD >= 3) IRNGT (3) = ILOWT (2)
         Return
      End If
!
      If (XDONT(NDON) > XDONT(ILOWT(1))) Then
         ILOWT (3) = ILOWT (2)
         ILOWT (2) = ILOWT (1)
         If (XDONT(NDON) > XDONT(IHIGT(1))) Then
            ILOWT (1) = IHIGT (1)
            IHIGT (1) = NDON
         Else
            ILOWT (1) = NDON
         End If
      Else
         if (XDONT (NDON) > XDONT (ILOWT(2))) Then
            ILOWT (3) = ILOWT (2)
            ILOWT (2) = NDON
         else
            ILOWT (3) = NDON
         endif 
      End If
!
      If (NDON < 5) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         If (NORD >= 3) IRNGT (3) = ILOWT (2)
         If (NORD >= 4) IRNGT (4) = ILOWT (3)
         Return
      End If
! ---
      JDEB = 0
      IDEB = JDEB + 1
      JHIG = IDEB
      JLOW = 3
      XPIV = XDONT (IHIGT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                   (XDONT(ILOWT(3))-XDONT(IHIGT(IDEB)))
      If (XPIV >= XDONT(ILOWT(1))) Then
         XPIV = XDONT (IHIGT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                      (XDONT(ILOWT(2))-XDONT(IHIGT(IDEB)))
         If (XPIV >= XDONT(ILOWT(1))) &
             XPIV = XDONT (IHIGT(IDEB)) + REAL (2*NORD) / REAL (NDON+NORD) * &
                                          (XDONT(ILOWT(1))-XDONT(IHIGT(IDEB)))
      End If
      XPIV0 = XPIV
! ---
!  One puts values < pivot in the end and those >= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the ILOWT array as soon as we have more
!  than enough values in IHIGT.
!
!
      If (XDONT(NDON) < XPIV) Then
         ICRS = 3
         Do
            ICRS = ICRS + 1
            If (XDONT(ICRS) < XPIV) Then
               If (ICRS >= NDON) Exit
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
            Else
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
               If (JHIG >= NORD) Exit
            End If
         End Do
!
!  One restricts further processing because it is no use
!  to store more low values
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) >= XPIV) Then
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ICRS
               Else If (ICRS >= NDON) Then
                  Exit
               End If
            End Do
         End If
!
! ---
      Else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         Do ICRS = 4, NDON - 1
            If (XDONT(ICRS) < XPIV) Then
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
            Else
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
               If (JHIG >= NORD) Exit
            End If
         End Do
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) >= XPIV) Then
                  If (ICRS >= NDON) Exit
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ICRS
               End If
            End Do
         End If
      End If
! ---
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      Do
         if (JHIG == NORD) Exit
         If (JHM2 == JHIG .And. JLM2 == JLOW) Then
!
!   We are oscillating. Perturbate by bringing JHIG closer by one
!   to NORD
!
           If (NORD > JHIG) Then
                XMAX = XDONT (ILOWT(1))
                ILOW = 1
                Do ICRS = 2, JLOW
                   If (XDONT(ILOWT(ICRS)) > XMAX) Then
                      XMAX = XDONT (ILOWT(ICRS))
                      ILOW = ICRS
                   End If
                End Do
!
                JHIG = JHIG + 1
                IHIGT (JHIG) = ILOWT (ILOW)
                ILOWT (ILOW) = ILOWT (JLOW)
                JLOW = JLOW - 1
             Else
                IHIG = IHIGT (JHIG)
                XMIN = XDONT (IHIG)
                Do ICRS = 1, JHIG
                   If (XDONT(IHIGT(ICRS)) < XMIN) Then
                      IWRK = IHIGT (ICRS)
                      XMIN = XDONT (IWRK)
                      IHIGT (ICRS) = IHIG
                      IHIG = IWRK
                   End If
                End Do
                JHIG = JHIG - 1
             End If
         End If
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
! ---
!   We try to bring the number of values in the high values set
!   closer to NORD.
!
        Select Case (NORD-JHIG)
         Case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            Select Case (JLOW)
!!!!!           CASE DEFAULT
!!!!!              write (*,*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the low values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            Case (2)
               If (XDONT(ILOWT(1)) >= XDONT(ILOWT(2))) Then
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (1)
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (2)
               Else
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (2)
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (1)
               End If
               Exit
! ---
            Case (3)
!
!
               IWRK1 = ILOWT (1)
               IWRK2 = ILOWT (2)
               IWRK3 = ILOWT (3)
               If (XDONT(IWRK2) > XDONT(IWRK1)) Then
                  ILOWT (1) = IWRK2
                  ILOWT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) < XDONT(IWRK3)) Then
                  ILOWT (3) = IWRK2
                  ILOWT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) > XDONT(ILOWT(1))) Then
                     ILOWT (2) = ILOWT (1)
                     ILOWT (1) = IWRK2
                  End If
               End If
               JLOW = 0
               Do ICRS = JHIG + 1, NORD
                  JLOW = JLOW + 1
                  IHIGT (ICRS) = ILOWT (JLOW)
               End Do
               JHIG = NORD
               Exit
! ---
            Case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JLOW
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = ILOWT (1)
               IWRK2 = ILOWT (2)
               IWRK3 = ILOWT (IFIN)
               If (XDONT(IWRK2) > XDONT(IWRK1)) Then
                  ILOWT (1) = IWRK2
                  ILOWT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) < XDONT(IWRK3)) Then
                  ILOWT (IFIN) = IWRK2
                  ILOWT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) > XDONT(IHIGT(1))) Then
                     ILOWT (2) = ILOWT (1)
                     ILOWT (1) = IWRK2
                  End If
               End If
!
               JDEB = JHIG
               NWRK = NORD - JHIG
               IWRK1 = ILOWT (1)
               JHIG = JHIG + 1
               IHIGT (JHIG) = IWRK1
               XPIV = XDONT (IWRK1) + REAL (NWRK) / REAL (NORD+NWRK) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(IWRK1))
!
!  One takes values >= pivot to IHIGT
!  Again, 2 parts, one where we take care of the remaining
!  low values because we might still need them, and the
!  other when we know that we will have more than enough
!  high values in the end.
! ---
               JLOW = 0
               Do ICRS = 2, IFIN
                  If (XDONT(ILOWT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     If (JHIG >= NORD) Exit
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(ILOWT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  End If
               End Do
           End Select
! ---
!
         Case (1)
!
!  Only 1 value is missing in high part
!
            XMAX = XDONT (ILOWT(1))
            ILOW = 1
            Do ICRS = 2, JLOW
               If (XDONT(ILOWT(ICRS)) > XMAX) Then
                  XMAX = XDONT (ILOWT(ICRS))
                  ILOW = ICRS
               End If
            End Do
!
            JHIG = JHIG + 1
            IHIGT (JHIG) = ILOWT (ILOW)
            Exit
!
!
         Case (0)
!
!  Low part is exactly what we want
!
            Exit
! ---
!
         Case (-5:-1)
!
!  Only few values too many in high part
!
            IRNGT (1) = IHIGT (1)
            Do ICRS = 2, NORD
               IWRK = IHIGT (ICRS)
               XWRK = XDONT (IWRK)
               Do IDCR = ICRS - 1, 1, - 1
                  If (XWRK > XDONT(IRNGT(IDCR))) Then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                     Exit
                  End If
               End Do
               IRNGT (IDCR+1) = IWRK
            End Do
!
            XWRK1 = XDONT (IRNGT(NORD))
            Do ICRS = NORD + 1, JHIG
               If (XDONT(IHIGT (ICRS)) > XWRK1) Then
                  XWRK = XDONT (IHIGT (ICRS))
                  Do IDCR = NORD - 1, 1, - 1
                     If (XWRK <= XDONT(IRNGT(IDCR))) Exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  End Do
                  IRNGT (IDCR+1) = IHIGT (ICRS)
                  XWRK1 = XDONT (IRNGT(NORD))
               End If
            End Do
!
            Return
!
!
         Case (:-6)
!
! last case: too many values in high part
! ---
            IDEB = JDEB + 1
            IMIL = (JHIG+IDEB) / 2
            IFIN = JHIG
! ---
!  One chooses a pivot from 1st, last, and middle values
!
            If (XDONT(IHIGT(IMIL)) > XDONT(IHIGT(IDEB))) Then
               IWRK = IHIGT (IDEB)
               IHIGT (IDEB) = IHIGT (IMIL)
               IHIGT (IMIL) = IWRK
            End If
            If (XDONT(IHIGT(IMIL)) < XDONT(IHIGT(IFIN))) Then
               IWRK = IHIGT (IFIN)
               IHIGT (IFIN) = IHIGT (IMIL)
               IHIGT (IMIL) = IWRK
               If (XDONT(IHIGT(IMIL)) > XDONT(IHIGT(IDEB))) Then
                  IWRK = IHIGT (IDEB)
                  IHIGT (IDEB) = IHIGT (IMIL)
                  IHIGT (IMIL) = IWRK
               End If
            End If
            If (IFIN <= 3) Exit
! ---
            XPIV = XDONT (IHIGT(1)) + REAL(NORD)/REAL(JHIG+NORD) * &
                                      (XDONT(IHIGT(IFIN))-XDONT(IHIGT(1)))
            If (JDEB > 0) Then
               If (XPIV <= XPIV0) &
                   XPIV = XPIV0 + REAL(2*NORD-JDEB)/REAL (JHIG+NORD) * &
                                  (XDONT(IHIGT(IFIN))-XPIV0)
            Else
               IDEB = 1
            End If
!
!  One takes values < XPIV to ILOWT
!  However, we do not process the first values if we have been
!  through the case when we did not have enough high values
! ---
            JLOW = 0
            JHIG = JDEB
! ---
            If (XDONT(IHIGT(IFIN)) < XPIV) Then
               ICRS = JDEB
               Do
                 ICRS = ICRS + 1
                  If (XDONT(IHIGT(ICRS)) < XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     If (ICRS >= IFIN) Exit
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                     If (JHIG >= NORD) Exit
                  End If
               End Do
! ---
               If (ICRS < IFIN) Then
                  Do
                     ICRS = ICRS + 1
                     If (XDONT(IHIGT(ICRS)) >= XPIV) Then
                        JHIG = JHIG + 1
                        IHIGT (JHIG) = IHIGT (ICRS)
                     Else
                        If (ICRS >= IFIN) Exit
                     End If
                  End Do
               End If
           Else
               Do ICRS = IDEB, IFIN
                  If (XDONT(IHIGT(ICRS)) < XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                     If (JHIG >= NORD) Exit
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(IHIGT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  End If
               End Do
            End If
!
         End Select
!
      End Do
! ---
!  Now, we only need to complete ranking of the 1:NORD set
!  Assuming NORD is small, we use a simple insertion sort
!
      IRNGT (1) = IHIGT (1)
      Do ICRS = 2, NORD
         IWRK = IHIGT (ICRS)
         XWRK = XDONT (IWRK)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK > XDONT(IRNGT(IDCR))) Then
               IRNGT (IDCR+1) = IRNGT (IDCR)
            Else
               Exit
            End If
         End Do
         IRNGT (IDCR+1) = IWRK
      End Do
     Return
!
!
End Subroutine D_rapknr

Subroutine R_rapknr (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD, in decreasing order.
!  rapknr = (rnkpar backwards)
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values larger than the pivot
!  (IHIGT), and the indices of values smaller than the pivot
!  that we might still need later on (ILOWT). It iterates
!  until it can bring the number of values in IHIGT to
!  exactly NORD, and then uses an insertion sort to rank
!  this set, since it is supposedly small.
!  Michel Olagnon - Feb. 2011
! __________________________________________________________
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real    :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      Integer, Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      Integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      Integer :: IDEB, JDEB, IMIL, IFIN, NWRK, ICRS, IDCR, ILOW
      Integer :: JLM2, JLM1, JHM2, JHM1
!
      NDON = SIZE (XDONT)
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      If (NDON < 2) Then
         If (NORD >= 1) IRNGT (1) = 1
         Return
      End If
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of high values.
!
      If (XDONT(2) < XDONT(1)) Then
         ILOWT (1) = 2
         IHIGT (1) = 1
      Else
         ILOWT (1) = 1
         IHIGT (1) = 2
      End If
!
      If (NDON < 3) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         Return
      End If
! ---
      If (XDONT(3) > XDONT(ILOWT(1))) Then
         ILOWT (2) = ILOWT (1)
         If (XDONT(3) > XDONT(IHIGT(1))) Then
            ILOWT (1) = IHIGT (1)
            IHIGT (1) = 3
         Else
            ILOWT (1) = 3
         End If
      Else
         ILOWT (2) = 3
      End If
! ---
      If (NDON < 4) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         If (NORD >= 3) IRNGT (3) = ILOWT (2)
         Return
      End If
!
      If (XDONT(NDON) > XDONT(ILOWT(1))) Then
         ILOWT (3) = ILOWT (2)
         ILOWT (2) = ILOWT (1)
         If (XDONT(NDON) > XDONT(IHIGT(1))) Then
            ILOWT (1) = IHIGT (1)
            IHIGT (1) = NDON
         Else
            ILOWT (1) = NDON
         End If
      Else
         if (XDONT (NDON) > XDONT (ILOWT(2))) Then
            ILOWT (3) = ILOWT (2)
            ILOWT (2) = NDON
         else
            ILOWT (3) = NDON
         endif 
      End If
!
      If (NDON < 5) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         If (NORD >= 3) IRNGT (3) = ILOWT (2)
         If (NORD >= 4) IRNGT (4) = ILOWT (3)
         Return
      End If
! ---
      JDEB = 0
      IDEB = JDEB + 1
      JHIG = IDEB
      JLOW = 3
      XPIV = XDONT (IHIGT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                   (XDONT(ILOWT(3))-XDONT(IHIGT(IDEB)))
      If (XPIV >= XDONT(ILOWT(1))) Then
         XPIV = XDONT (IHIGT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                      (XDONT(ILOWT(2))-XDONT(IHIGT(IDEB)))
         If (XPIV >= XDONT(ILOWT(1))) &
             XPIV = XDONT (IHIGT(IDEB)) + REAL (2*NORD) / REAL (NDON+NORD) * &
                                          (XDONT(ILOWT(1))-XDONT(IHIGT(IDEB)))
      End If
      XPIV0 = XPIV
! ---
!  One puts values < pivot in the end and those >= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the ILOWT array as soon as we have more
!  than enough values in IHIGT.
!
!
      If (XDONT(NDON) < XPIV) Then
         ICRS = 3
         Do
            ICRS = ICRS + 1
            If (XDONT(ICRS) < XPIV) Then
               If (ICRS >= NDON) Exit
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
            Else
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
               If (JHIG >= NORD) Exit
            End If
         End Do
!
!  One restricts further processing because it is no use
!  to store more low values
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) >= XPIV) Then
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ICRS
               Else If (ICRS >= NDON) Then
                  Exit
               End If
            End Do
         End If
!
! ---
      Else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         Do ICRS = 4, NDON - 1
            If (XDONT(ICRS) < XPIV) Then
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
            Else
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
               If (JHIG >= NORD) Exit
            End If
         End Do
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) >= XPIV) Then
                  If (ICRS >= NDON) Exit
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ICRS
               End If
            End Do
         End If
      End If
! ---
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      Do
         if (JHIG == NORD) Exit
         If (JHM2 == JHIG .And. JLM2 == JLOW) Then
!
!   We are oscillating. Perturbate by bringing JHIG closer by one
!   to NORD
!
           If (NORD > JHIG) Then
                XMAX = XDONT (ILOWT(1))
                ILOW = 1
                Do ICRS = 2, JLOW
                   If (XDONT(ILOWT(ICRS)) > XMAX) Then
                      XMAX = XDONT (ILOWT(ICRS))
                      ILOW = ICRS
                   End If
                End Do
!
                JHIG = JHIG + 1
                IHIGT (JHIG) = ILOWT (ILOW)
                ILOWT (ILOW) = ILOWT (JLOW)
                JLOW = JLOW - 1
             Else
                IHIG = IHIGT (JHIG)
                XMIN = XDONT (IHIG)
                Do ICRS = 1, JHIG
                   If (XDONT(IHIGT(ICRS)) < XMIN) Then
                      IWRK = IHIGT (ICRS)
                      XMIN = XDONT (IWRK)
                      IHIGT (ICRS) = IHIG
                      IHIG = IWRK
                   End If
                End Do
                JHIG = JHIG - 1
             End If
         End If
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
! ---
!   We try to bring the number of values in the high values set
!   closer to NORD.
!
        Select Case (NORD-JHIG)
         Case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            Select Case (JLOW)
!!!!!           CASE DEFAULT
!!!!!              write (*,*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the low values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            Case (2)
               If (XDONT(ILOWT(1)) >= XDONT(ILOWT(2))) Then
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (1)
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (2)
               Else
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (2)
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (1)
               End If
               Exit
! ---
            Case (3)
!
!
               IWRK1 = ILOWT (1)
               IWRK2 = ILOWT (2)
               IWRK3 = ILOWT (3)
               If (XDONT(IWRK2) > XDONT(IWRK1)) Then
                  ILOWT (1) = IWRK2
                  ILOWT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) < XDONT(IWRK3)) Then
                  ILOWT (3) = IWRK2
                  ILOWT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) > XDONT(ILOWT(1))) Then
                     ILOWT (2) = ILOWT (1)
                     ILOWT (1) = IWRK2
                  End If
               End If
               JLOW = 0
               Do ICRS = JHIG + 1, NORD
                  JLOW = JLOW + 1
                  IHIGT (ICRS) = ILOWT (JLOW)
               End Do
               JHIG = NORD
               Exit
! ---
            Case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JLOW
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = ILOWT (1)
               IWRK2 = ILOWT (2)
               IWRK3 = ILOWT (IFIN)
               If (XDONT(IWRK2) > XDONT(IWRK1)) Then
                  ILOWT (1) = IWRK2
                  ILOWT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) < XDONT(IWRK3)) Then
                  ILOWT (IFIN) = IWRK2
                  ILOWT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) > XDONT(IHIGT(1))) Then
                     ILOWT (2) = ILOWT (1)
                     ILOWT (1) = IWRK2
                  End If
               End If
!
               JDEB = JHIG
               NWRK = NORD - JHIG
               IWRK1 = ILOWT (1)
               JHIG = JHIG + 1
               IHIGT (JHIG) = IWRK1
               XPIV = XDONT (IWRK1) + REAL (NWRK) / REAL (NORD+NWRK) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(IWRK1))
!
!  One takes values >= pivot to IHIGT
!  Again, 2 parts, one where we take care of the remaining
!  low values because we might still need them, and the
!  other when we know that we will have more than enough
!  high values in the end.
! ---
               JLOW = 0
               Do ICRS = 2, IFIN
                  If (XDONT(ILOWT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     If (JHIG >= NORD) Exit
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(ILOWT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  End If
               End Do
           End Select
! ---
!
         Case (1)
!
!  Only 1 value is missing in high part
!
            XMAX = XDONT (ILOWT(1))
            ILOW = 1
            Do ICRS = 2, JLOW
               If (XDONT(ILOWT(ICRS)) > XMAX) Then
                  XMAX = XDONT (ILOWT(ICRS))
                  ILOW = ICRS
               End If
            End Do
!
            JHIG = JHIG + 1
            IHIGT (JHIG) = ILOWT (ILOW)
            Exit
!
!
         Case (0)
!
!  Low part is exactly what we want
!
            Exit
! ---
!
         Case (-5:-1)
!
!  Only few values too many in high part
!
            IRNGT (1) = IHIGT (1)
            Do ICRS = 2, NORD
               IWRK = IHIGT (ICRS)
               XWRK = XDONT (IWRK)
               Do IDCR = ICRS - 1, 1, - 1
                  If (XWRK > XDONT(IRNGT(IDCR))) Then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                     Exit
                  End If
               End Do
               IRNGT (IDCR+1) = IWRK
            End Do
!
            XWRK1 = XDONT (IRNGT(NORD))
            Do ICRS = NORD + 1, JHIG
               If (XDONT(IHIGT (ICRS)) > XWRK1) Then
                  XWRK = XDONT (IHIGT (ICRS))
                  Do IDCR = NORD - 1, 1, - 1
                     If (XWRK <= XDONT(IRNGT(IDCR))) Exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  End Do
                  IRNGT (IDCR+1) = IHIGT (ICRS)
                  XWRK1 = XDONT (IRNGT(NORD))
               End If
            End Do
!
            Return
!
!
         Case (:-6)
!
! last case: too many values in high part
! ---
            IDEB = JDEB + 1
            IMIL = (JHIG+IDEB) / 2
            IFIN = JHIG
! ---
!  One chooses a pivot from 1st, last, and middle values
!
            If (XDONT(IHIGT(IMIL)) > XDONT(IHIGT(IDEB))) Then
               IWRK = IHIGT (IDEB)
               IHIGT (IDEB) = IHIGT (IMIL)
               IHIGT (IMIL) = IWRK
            End If
            If (XDONT(IHIGT(IMIL)) < XDONT(IHIGT(IFIN))) Then
               IWRK = IHIGT (IFIN)
               IHIGT (IFIN) = IHIGT (IMIL)
               IHIGT (IMIL) = IWRK
               If (XDONT(IHIGT(IMIL)) > XDONT(IHIGT(IDEB))) Then
                  IWRK = IHIGT (IDEB)
                  IHIGT (IDEB) = IHIGT (IMIL)
                  IHIGT (IMIL) = IWRK
               End If
            End If
            If (IFIN <= 3) Exit
! ---
            XPIV = XDONT (IHIGT(1)) + REAL(NORD)/REAL(JHIG+NORD) * &
                                      (XDONT(IHIGT(IFIN))-XDONT(IHIGT(1)))
            If (JDEB > 0) Then
               If (XPIV <= XPIV0) &
                   XPIV = XPIV0 + REAL(2*NORD-JDEB)/REAL (JHIG+NORD) * &
                                  (XDONT(IHIGT(IFIN))-XPIV0)
            Else
               IDEB = 1
            End If
!
!  One takes values < XPIV to ILOWT
!  However, we do not process the first values if we have been
!  through the case when we did not have enough high values
! ---
            JLOW = 0
            JHIG = JDEB
! ---
            If (XDONT(IHIGT(IFIN)) < XPIV) Then
               ICRS = JDEB
               Do
                 ICRS = ICRS + 1
                  If (XDONT(IHIGT(ICRS)) < XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     If (ICRS >= IFIN) Exit
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                     If (JHIG >= NORD) Exit
                  End If
               End Do
! ---
               If (ICRS < IFIN) Then
                  Do
                     ICRS = ICRS + 1
                     If (XDONT(IHIGT(ICRS)) >= XPIV) Then
                        JHIG = JHIG + 1
                        IHIGT (JHIG) = IHIGT (ICRS)
                     Else
                        If (ICRS >= IFIN) Exit
                     End If
                  End Do
               End If
           Else
               Do ICRS = IDEB, IFIN
                  If (XDONT(IHIGT(ICRS)) < XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                     If (JHIG >= NORD) Exit
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(IHIGT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  End If
               End Do
            End If
!
         End Select
!
      End Do
! ---
!  Now, we only need to complete ranking of the 1:NORD set
!  Assuming NORD is small, we use a simple insertion sort
!
      IRNGT (1) = IHIGT (1)
      Do ICRS = 2, NORD
         IWRK = IHIGT (ICRS)
         XWRK = XDONT (IWRK)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK > XDONT(IRNGT(IDCR))) Then
               IRNGT (IDCR+1) = IRNGT (IDCR)
            Else
               Exit
            End If
         End Do
         IRNGT (IDCR+1) = IWRK
      End Do
     Return
!
!
End Subroutine R_rapknr
Subroutine I_rapknr (XDONT, IRNGT, NORD)
!  Ranks partially XDONT by IRNGT, up to order NORD, in decreasing order.
!  rapknr = (rnkpar backwards)
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values larger than the pivot
!  (IHIGT), and the indices of values smaller than the pivot
!  that we might still need later on (ILOWT). It iterates
!  until it can bring the number of values in IHIGT to
!  exactly NORD, and then uses an insertion sort to rank
!  this set, since it is supposedly small.
!  Michel Olagnon - Feb. 2011
! __________________________________________________________
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Integer :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      Integer, Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      Integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      Integer :: IDEB, JDEB, IMIL, IFIN, NWRK, ICRS, IDCR, ILOW
      Integer :: JLM2, JLM1, JHM2, JHM1
!
      NDON = SIZE (XDONT)
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      If (NDON < 2) Then
         If (NORD >= 1) IRNGT (1) = 1
         Return
      End If
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of high values.
!
      If (XDONT(2) < XDONT(1)) Then
         ILOWT (1) = 2
         IHIGT (1) = 1
      Else
         ILOWT (1) = 1
         IHIGT (1) = 2
      End If
!
      If (NDON < 3) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         Return
      End If
! ---
      If (XDONT(3) > XDONT(ILOWT(1))) Then
         ILOWT (2) = ILOWT (1)
         If (XDONT(3) > XDONT(IHIGT(1))) Then
            ILOWT (1) = IHIGT (1)
            IHIGT (1) = 3
         Else
            ILOWT (1) = 3
         End If
      Else
         ILOWT (2) = 3
      End If
! ---
      If (NDON < 4) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         If (NORD >= 3) IRNGT (3) = ILOWT (2)
         Return
      End If
!
      If (XDONT(NDON) > XDONT(ILOWT(1))) Then
         ILOWT (3) = ILOWT (2)
         ILOWT (2) = ILOWT (1)
         If (XDONT(NDON) > XDONT(IHIGT(1))) Then
            ILOWT (1) = IHIGT (1)
            IHIGT (1) = NDON
         Else
            ILOWT (1) = NDON
         End If
      Else
         if (XDONT (NDON) > XDONT (ILOWT(2))) Then
            ILOWT (3) = ILOWT (2)
            ILOWT (2) = NDON
         else
            ILOWT (3) = NDON
         endif 
      End If
!
      If (NDON < 5) Then
         If (NORD >= 1) IRNGT (1) = IHIGT (1)
         If (NORD >= 2) IRNGT (2) = ILOWT (1)
         If (NORD >= 3) IRNGT (3) = ILOWT (2)
         If (NORD >= 4) IRNGT (4) = ILOWT (3)
         Return
      End If
! ---
      JDEB = 0
      IDEB = JDEB + 1
      JHIG = IDEB
      JLOW = 3
      XPIV = XDONT (IHIGT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                   (XDONT(ILOWT(3))-XDONT(IHIGT(IDEB)))
      If (XPIV >= XDONT(ILOWT(1))) Then
         XPIV = XDONT (IHIGT(IDEB)) + REAL(2*NORD)/REAL(NDON+NORD) * &
                                      (XDONT(ILOWT(2))-XDONT(IHIGT(IDEB)))
         If (XPIV >= XDONT(ILOWT(1))) &
             XPIV = XDONT (IHIGT(IDEB)) + REAL (2*NORD) / REAL (NDON+NORD) * &
                                          (XDONT(ILOWT(1))-XDONT(IHIGT(IDEB)))
      End If
      XPIV0 = XPIV
! ---
!  One puts values < pivot in the end and those >= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the ILOWT array as soon as we have more
!  than enough values in IHIGT.
!
!
      If (XDONT(NDON) < XPIV) Then
         ICRS = 3
         Do
            ICRS = ICRS + 1
            If (XDONT(ICRS) < XPIV) Then
               If (ICRS >= NDON) Exit
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
            Else
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
               If (JHIG >= NORD) Exit
            End If
         End Do
!
!  One restricts further processing because it is no use
!  to store more low values
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) >= XPIV) Then
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ICRS
               Else If (ICRS >= NDON) Then
                  Exit
               End If
            End Do
         End If
!
! ---
      Else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         Do ICRS = 4, NDON - 1
            If (XDONT(ICRS) < XPIV) Then
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
            Else
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
               If (JHIG >= NORD) Exit
            End If
         End Do
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) >= XPIV) Then
                  If (ICRS >= NDON) Exit
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ICRS
               End If
            End Do
         End If
      End If
! ---
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      Do
         if (JHIG == NORD) Exit
         If (JHM2 == JHIG .And. JLM2 == JLOW) Then
!
!   We are oscillating. Perturbate by bringing JHIG closer by one
!   to NORD
!
           If (NORD > JHIG) Then
                XMAX = XDONT (ILOWT(1))
                ILOW = 1
                Do ICRS = 2, JLOW
                   If (XDONT(ILOWT(ICRS)) > XMAX) Then
                      XMAX = XDONT (ILOWT(ICRS))
                      ILOW = ICRS
                   End If
                End Do
!
                JHIG = JHIG + 1
                IHIGT (JHIG) = ILOWT (ILOW)
                ILOWT (ILOW) = ILOWT (JLOW)
                JLOW = JLOW - 1
             Else
                IHIG = IHIGT (JHIG)
                XMIN = XDONT (IHIG)
                Do ICRS = 1, JHIG
                   If (XDONT(IHIGT(ICRS)) < XMIN) Then
                      IWRK = IHIGT (ICRS)
                      XMIN = XDONT (IWRK)
                      IHIGT (ICRS) = IHIG
                      IHIG = IWRK
                   End If
                End Do
                JHIG = JHIG - 1
             End If
         End If
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
! ---
!   We try to bring the number of values in the high values set
!   closer to NORD.
!
        Select Case (NORD-JHIG)
         Case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            Select Case (JLOW)
!!!!!           CASE DEFAULT
!!!!!              write (*,*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the low values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            Case (2)
               If (XDONT(ILOWT(1)) >= XDONT(ILOWT(2))) Then
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (1)
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (2)
               Else
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (2)
                  JHIG = JHIG + 1
                  IHIGT (JHIG) = ILOWT (1)
               End If
               Exit
! ---
            Case (3)
!
!
               IWRK1 = ILOWT (1)
               IWRK2 = ILOWT (2)
               IWRK3 = ILOWT (3)
               If (XDONT(IWRK2) > XDONT(IWRK1)) Then
                  ILOWT (1) = IWRK2
                  ILOWT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) < XDONT(IWRK3)) Then
                  ILOWT (3) = IWRK2
                  ILOWT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) > XDONT(ILOWT(1))) Then
                     ILOWT (2) = ILOWT (1)
                     ILOWT (1) = IWRK2
                  End If
               End If
               JLOW = 0
               Do ICRS = JHIG + 1, NORD
                  JLOW = JLOW + 1
                  IHIGT (ICRS) = ILOWT (JLOW)
               End Do
               JHIG = NORD
               Exit
! ---
            Case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JLOW
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = ILOWT (1)
               IWRK2 = ILOWT (2)
               IWRK3 = ILOWT (IFIN)
               If (XDONT(IWRK2) > XDONT(IWRK1)) Then
                  ILOWT (1) = IWRK2
                  ILOWT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) < XDONT(IWRK3)) Then
                  ILOWT (IFIN) = IWRK2
                  ILOWT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) > XDONT(IHIGT(1))) Then
                     ILOWT (2) = ILOWT (1)
                     ILOWT (1) = IWRK2
                  End If
               End If
!
               JDEB = JHIG
               NWRK = NORD - JHIG
               IWRK1 = ILOWT (1)
               JHIG = JHIG + 1
               IHIGT (JHIG) = IWRK1
               XPIV = XDONT (IWRK1) + REAL (NWRK) / REAL (NORD+NWRK) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(IWRK1))
!
!  One takes values >= pivot to IHIGT
!  Again, 2 parts, one where we take care of the remaining
!  low values because we might still need them, and the
!  other when we know that we will have more than enough
!  high values in the end.
! ---
               JLOW = 0
               Do ICRS = 2, IFIN
                  If (XDONT(ILOWT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     If (JHIG >= NORD) Exit
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(ILOWT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  End If
               End Do
           End Select
! ---
!
         Case (1)
!
!  Only 1 value is missing in high part
!
            XMAX = XDONT (ILOWT(1))
            ILOW = 1
            Do ICRS = 2, JLOW
               If (XDONT(ILOWT(ICRS)) > XMAX) Then
                  XMAX = XDONT (ILOWT(ICRS))
                  ILOW = ICRS
               End If
            End Do
!
            JHIG = JHIG + 1
            IHIGT (JHIG) = ILOWT (ILOW)
            Exit
!
!
         Case (0)
!
!  Low part is exactly what we want
!
            Exit
! ---
!
         Case (-5:-1)
!
!  Only few values too many in high part
!
            IRNGT (1) = IHIGT (1)
            Do ICRS = 2, NORD
               IWRK = IHIGT (ICRS)
               XWRK = XDONT (IWRK)
               Do IDCR = ICRS - 1, 1, - 1
                  If (XWRK > XDONT(IRNGT(IDCR))) Then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                     Exit
                  End If
               End Do
               IRNGT (IDCR+1) = IWRK
            End Do
!
            XWRK1 = XDONT (IRNGT(NORD))
            Do ICRS = NORD + 1, JHIG
               If (XDONT(IHIGT (ICRS)) > XWRK1) Then
                  XWRK = XDONT (IHIGT (ICRS))
                  Do IDCR = NORD - 1, 1, - 1
                     If (XWRK <= XDONT(IRNGT(IDCR))) Exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  End Do
                  IRNGT (IDCR+1) = IHIGT (ICRS)
                  XWRK1 = XDONT (IRNGT(NORD))
               End If
            End Do
!
            Return
!
!
         Case (:-6)
!
! last case: too many values in high part
! ---
            IDEB = JDEB + 1
            IMIL = (JHIG+IDEB) / 2
            IFIN = JHIG
! ---
!  One chooses a pivot from 1st, last, and middle values
!
            If (XDONT(IHIGT(IMIL)) > XDONT(IHIGT(IDEB))) Then
               IWRK = IHIGT (IDEB)
               IHIGT (IDEB) = IHIGT (IMIL)
               IHIGT (IMIL) = IWRK
            End If
            If (XDONT(IHIGT(IMIL)) < XDONT(IHIGT(IFIN))) Then
               IWRK = IHIGT (IFIN)
               IHIGT (IFIN) = IHIGT (IMIL)
               IHIGT (IMIL) = IWRK
               If (XDONT(IHIGT(IMIL)) > XDONT(IHIGT(IDEB))) Then
                  IWRK = IHIGT (IDEB)
                  IHIGT (IDEB) = IHIGT (IMIL)
                  IHIGT (IMIL) = IWRK
               End If
            End If
            If (IFIN <= 3) Exit
! ---
            XPIV = XDONT (IHIGT(1)) + REAL(NORD)/REAL(JHIG+NORD) * &
                                      (XDONT(IHIGT(IFIN))-XDONT(IHIGT(1)))
            If (JDEB > 0) Then
               If (XPIV <= XPIV0) &
                   XPIV = XPIV0 + REAL(2*NORD-JDEB)/REAL (JHIG+NORD) * &
                                  (XDONT(IHIGT(IFIN))-XPIV0)
            Else
               IDEB = 1
            End If
!
!  One takes values < XPIV to ILOWT
!  However, we do not process the first values if we have been
!  through the case when we did not have enough high values
! ---
            JLOW = 0
            JHIG = JDEB
! ---
            If (XDONT(IHIGT(IFIN)) < XPIV) Then
               ICRS = JDEB
               Do
                 ICRS = ICRS + 1
                  If (XDONT(IHIGT(ICRS)) < XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     If (ICRS >= IFIN) Exit
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                     If (JHIG >= NORD) Exit
                  End If
               End Do
! ---
               If (ICRS < IFIN) Then
                  Do
                     ICRS = ICRS + 1
                     If (XDONT(IHIGT(ICRS)) >= XPIV) Then
                        JHIG = JHIG + 1
                        IHIGT (JHIG) = IHIGT (ICRS)
                     Else
                        If (ICRS >= IFIN) Exit
                     End If
                  End Do
               End If
           Else
               Do ICRS = IDEB, IFIN
                  If (XDONT(IHIGT(ICRS)) < XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                     If (JHIG >= NORD) Exit
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(IHIGT(ICRS)) >= XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  End If
               End Do
            End If
!
         End Select
!
      End Do
! ---
!  Now, we only need to complete ranking of the 1:NORD set
!  Assuming NORD is small, we use a simple insertion sort
!
      IRNGT (1) = IHIGT (1)
      Do ICRS = 2, NORD
         IWRK = IHIGT (ICRS)
         XWRK = XDONT (IWRK)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK > XDONT(IRNGT(IDCR))) Then
               IRNGT (IDCR+1) = IRNGT (IDCR)
            Else
               Exit
            End If
         End Do
         IRNGT (IDCR+1) = IWRK
      End Do
     Return
!
!
End Subroutine I_rapknr
end module m_rapknr
