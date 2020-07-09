Module m_indnth
Integer, Parameter :: kdp = selected_real_kind(15)
public :: indnth
private :: kdp
private :: R_indnth, I_indnth, D_indnth
interface indnth
  module procedure d_indnth, r_indnth, i_indnth
end interface indnth
contains

Function D_indnth (XDONT, NORD) Result (INDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! __________________________________________________________
      Real (kind=kdp), Dimension (:), Intent (In) :: XDONT
      Integer :: INDNTH
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real (kind=kdp) :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      Integer, Dimension (NORD) :: IRNGT
      Integer, Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      Integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      Integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      Integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = NORD
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      If (NDON < 2) Then
         If (INTH == 1) INDNTH = 1
         Return
      End If
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
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
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         Return
      End If
!
      If (XDONT(3) < XDONT(IHIGT(1))) Then
         IHIGT (2) = IHIGT (1)
         If (XDONT(3) < XDONT(ILOWT(1))) Then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = 3
         Else
            IHIGT (1) = 3
         End If
      Else
         IHIGT (2) = 3
      End If
!
      If (NDON < 4) Then
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         If (INTH == 3) INDNTH = IHIGT (2)
         Return
      End If
!
      If (XDONT(NDON) < XDONT(IHIGT(1))) Then
         IHIGT (3) = IHIGT (2)
         IHIGT (2) = IHIGT (1)
         If (XDONT(NDON) < XDONT(ILOWT(1))) Then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = NDON
         Else
            IHIGT (1) = NDON
         End If
      Else
         IHIGT (3) = NDON
      End If
!
      If (NDON < 5) Then
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         If (INTH == 3) INDNTH = IHIGT (2)
         If (INTH == 4) INDNTH = IHIGT (3)
         Return
      End If
!

      JLOW = 1
      JHIG = 3
      XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                   (XDONT(IHIGT(3))-XDONT(ILOWT(1)))
      If (XPIV >= XDONT(IHIGT(1))) Then
         XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                      (XDONT(IHIGT(2))-XDONT(ILOWT(1)))
         If (XPIV >= XDONT(IHIGT(1))) &
             XPIV = XDONT (ILOWT(1)) + REAL (2*INTH) / REAL (NDON+INTH) * &
                                          (XDONT(IHIGT(1))-XDONT(ILOWT(1)))
      End If
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the IHIGT array as soon as we have more
!  than enough values in ILOWT.
!
!
      If (XDONT(NDON) > XPIV) Then
         ICRS = 3
         Do
            ICRS = ICRS + 1
            If (XDONT(ICRS) > XPIV) Then
               If (ICRS >= NDON) Exit
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            Else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               If (JLOW >= INTH) Exit
            End If
         End Do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) <= XPIV) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               Else If (ICRS >= NDON) Then
                  Exit
               End If
            End Do
         End If
!
!
      Else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         Do ICRS = 4, NDON - 1
            If (XDONT(ICRS) > XPIV) Then
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            Else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               If (JLOW >= INTH) Exit
            End If
         End Do
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) <= XPIV) Then
                  If (ICRS >= NDON) Exit
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               End If
            End Do
         End If
      End If
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      Do
         If (JLM2 == JLOW .And. JHM2 == JHIG) Then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             If (INTH > JLOW) Then
                XMIN = XDONT (IHIGT(1))
                IHIG = 1
                Do ICRS = 2, JHIG
                   If (XDONT(IHIGT(ICRS)) < XMIN) Then
                      XMIN = XDONT (IHIGT(ICRS))
                      IHIG = ICRS
                   End If
                End Do
!
                JLOW = JLOW + 1
                ILOWT (JLOW) = IHIGT (IHIG)
                IHIGT (IHIG) = IHIGT (JHIG)
                JHIG = JHIG - 1
             Else

                ILOW = ILOWT (1)
                XMAX = XDONT (ILOW)
                Do ICRS = 2, JLOW
                   If (XDONT(ILOWT(ICRS)) > XMAX) Then
                      IWRK = ILOWT (ICRS)
                      XMAX = XDONT (IWRK)
                      ILOWT (ICRS) = ILOW
                      ILOW = IWRK
                   End If
                End Do
                JLOW = JLOW - 1
             End If
         End If
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         Select Case (INTH-JLOW)
         Case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            Select Case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            Case (2)
               If (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
               Else
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
               End If
               Exit
!
            Case (3)
!
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (3)
               If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                  IHIGT (3) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  End If
               End If
               JHIG = 0
               Do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  ILOWT (ICRS) = IHIGT (JHIG)
               End Do
               JLOW = INTH
               Exit
!
            Case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (IFIN)
               If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                  IHIGT (IFIN) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  End If
               End If
!
               IWRK1 = IHIGT (1)
               JLOW = JLOW + 1
               ILOWT (JLOW) = IWRK1
               XPIV = XDONT (IWRK1) + 0.5 * (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
!
!  One takes values <= pivot to ILOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               Do ICRS = 2, IFIN
                  If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     If (JLOW >= INTH) Exit
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  End If
               End Do
            End Select
!
!
         Case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XDONT (IHIGT(1))
            IHIG = 1
            Do ICRS = 2, JHIG
               If (XDONT(IHIGT(ICRS)) < XMIN) Then
                  XMIN = XDONT (IHIGT(ICRS))
                  IHIG = ICRS
               End If
            End Do
!
            INDNTH = IHIGT (IHIG)
            Return
!
!
         Case (0)
!
!  Low part is exactly what we want
!
            Exit
!
!
         Case (-5:-1)
!
!  Only few values too many in low part
!
            IRNGT (1) = ILOWT (1)
            ILOW = 1 + INTH - JLOW
            Do ICRS = 2, INTH
               IWRK = ILOWT (ICRS)
               XWRK = XDONT (IWRK)
               Do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  If (XWRK < XDONT(IRNGT(IDCR))) Then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                     Exit
                  End If
               End Do
               IRNGT (IDCR+1) = IWRK
               ILOW = ILOW + 1
            End Do
!
            XWRK1 = XDONT (IRNGT(INTH))
            ILOW = 2*INTH - JLOW
            Do ICRS = INTH + 1, JLOW
               If (XDONT(ILOWT (ICRS)) < XWRK1) Then
                  XWRK = XDONT (ILOWT (ICRS))
                  Do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     If (XWRK >= XDONT(IRNGT(IDCR))) Exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  End Do
                  IRNGT (IDCR+1) = ILOWT (ICRS)
                  XWRK1 = XDONT (IRNGT(INTH))
               End If
               ILOW = ILOW + 1
            End Do
!
            INDNTH = IRNGT(INTH)
            Return
!
!
         Case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) Then
               IWRK = ILOWT (1)
               ILOWT (1) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
            End If
            If (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) Then
               IWRK = ILOWT (IFIN)
               ILOWT (IFIN) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
               If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) Then
                  IWRK = ILOWT (1)
                  ILOWT (1) = ILOWT (IMIL)
                  ILOWT (IMIL) = IWRK
               End If
            End If
            If (IFIN <= 3) Exit
!
            XPIV = XDONT (ILOWT(1)) + REAL(INTH)/REAL(JLOW+INTH) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))

!
!  One takes values > XPIV to IHIGT
!
            JHIG = 0
            JLOW = 0
!
            If (XDONT(ILOWT(IFIN)) > XPIV) Then
               ICRS = 0
               Do
                  ICRS = ICRS + 1
                  If (XDONT(ILOWT(ICRS)) > XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     If (ICRS >= IFIN) Exit
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     If (JLOW >= INTH) Exit
                  End If
               End Do
!
               If (ICRS < IFIN) Then
                  Do
                     ICRS = ICRS + 1
                     If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ILOWT (ICRS)
                     Else
                        If (ICRS >= IFIN) Exit
                     End If
                  End Do
               End If
            Else
               Do ICRS = 1, IFIN
                  If (XDONT(ILOWT(ICRS)) > XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     If (JLOW >= INTH) Exit
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  End If
               End Do
            End If
!
         End Select
!
      End Do
!
!  Now, we only need to find maximum of the 1:INTH set
!

      IWRK1 = ILOWT (1)
      XWRK1 =  XDONT (IWRK1)
      Do ICRS = 1+1, INTH
         IWRK = ILOWT (ICRS)
         XWRK = XDONT (IWRK)
         If (XWRK > XWRK1) Then
            XWRK1 = XWRK
            IWRK1 = IWRK
         End If
      End Do
      INDNTH = IWRK1
      Return
!
!
End Function D_indnth

Function R_indnth (XDONT, NORD) Result (INDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      Integer :: INDNTH
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Real :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      Integer, Dimension (NORD) :: IRNGT
      Integer, Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      Integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      Integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      Integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = NORD
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      If (NDON < 2) Then
         If (INTH == 1) INDNTH = 1
         Return
      End If
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
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
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         Return
      End If
!
      If (XDONT(3) < XDONT(IHIGT(1))) Then
         IHIGT (2) = IHIGT (1)
         If (XDONT(3) < XDONT(ILOWT(1))) Then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = 3
         Else
            IHIGT (1) = 3
         End If
      Else
         IHIGT (2) = 3
      End If
!
      If (NDON < 4) Then
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         If (INTH == 3) INDNTH = IHIGT (2)
         Return
      End If
!
      If (XDONT(NDON) < XDONT(IHIGT(1))) Then
         IHIGT (3) = IHIGT (2)
         IHIGT (2) = IHIGT (1)
         If (XDONT(NDON) < XDONT(ILOWT(1))) Then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = NDON
         Else
            IHIGT (1) = NDON
         End If
      Else
         IHIGT (3) = NDON
      End If
!
      If (NDON < 5) Then
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         If (INTH == 3) INDNTH = IHIGT (2)
         If (INTH == 4) INDNTH = IHIGT (3)
         Return
      End If
!

      JLOW = 1
      JHIG = 3
      XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                   (XDONT(IHIGT(3))-XDONT(ILOWT(1)))
      If (XPIV >= XDONT(IHIGT(1))) Then
         XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                      (XDONT(IHIGT(2))-XDONT(ILOWT(1)))
         If (XPIV >= XDONT(IHIGT(1))) &
             XPIV = XDONT (ILOWT(1)) + REAL (2*INTH) / REAL (NDON+INTH) * &
                                          (XDONT(IHIGT(1))-XDONT(ILOWT(1)))
      End If
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the IHIGT array as soon as we have more
!  than enough values in ILOWT.
!
!
      If (XDONT(NDON) > XPIV) Then
         ICRS = 3
         Do
            ICRS = ICRS + 1
            If (XDONT(ICRS) > XPIV) Then
               If (ICRS >= NDON) Exit
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            Else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               If (JLOW >= INTH) Exit
            End If
         End Do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) <= XPIV) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               Else If (ICRS >= NDON) Then
                  Exit
               End If
            End Do
         End If
!
!
      Else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         Do ICRS = 4, NDON - 1
            If (XDONT(ICRS) > XPIV) Then
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            Else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               If (JLOW >= INTH) Exit
            End If
         End Do
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) <= XPIV) Then
                  If (ICRS >= NDON) Exit
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               End If
            End Do
         End If
      End If
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      Do
         If (JLM2 == JLOW .And. JHM2 == JHIG) Then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             If (INTH > JLOW) Then
                XMIN = XDONT (IHIGT(1))
                IHIG = 1
                Do ICRS = 2, JHIG
                   If (XDONT(IHIGT(ICRS)) < XMIN) Then
                      XMIN = XDONT (IHIGT(ICRS))
                      IHIG = ICRS
                   End If
                End Do
!
                JLOW = JLOW + 1
                ILOWT (JLOW) = IHIGT (IHIG)
                IHIGT (IHIG) = IHIGT (JHIG)
                JHIG = JHIG - 1
             Else

                ILOW = ILOWT (1)
                XMAX = XDONT (ILOW)
                Do ICRS = 2, JLOW
                   If (XDONT(ILOWT(ICRS)) > XMAX) Then
                      IWRK = ILOWT (ICRS)
                      XMAX = XDONT (IWRK)
                      ILOWT (ICRS) = ILOW
                      ILOW = IWRK
                   End If
                End Do
                JLOW = JLOW - 1
             End If
         End If
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         Select Case (INTH-JLOW)
         Case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            Select Case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            Case (2)
               If (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
               Else
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
               End If
               Exit
!
            Case (3)
!
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (3)
               If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                  IHIGT (3) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  End If
               End If
               JHIG = 0
               Do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  ILOWT (ICRS) = IHIGT (JHIG)
               End Do
               JLOW = INTH
               Exit
!
            Case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (IFIN)
               If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                  IHIGT (IFIN) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  End If
               End If
!
               IWRK1 = IHIGT (1)
               JLOW = JLOW + 1
               ILOWT (JLOW) = IWRK1
               XPIV = XDONT (IWRK1) + 0.5 * (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
!
!  One takes values <= pivot to ILOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               Do ICRS = 2, IFIN
                  If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     If (JLOW >= INTH) Exit
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  End If
               End Do
            End Select
!
!
         Case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XDONT (IHIGT(1))
            IHIG = 1
            Do ICRS = 2, JHIG
               If (XDONT(IHIGT(ICRS)) < XMIN) Then
                  XMIN = XDONT (IHIGT(ICRS))
                  IHIG = ICRS
               End If
            End Do
!
            INDNTH = IHIGT (IHIG)
            Return
!
!
         Case (0)
!
!  Low part is exactly what we want
!
            Exit
!
!
         Case (-5:-1)
!
!  Only few values too many in low part
!
            IRNGT (1) = ILOWT (1)
            ILOW = 1 + INTH - JLOW
            Do ICRS = 2, INTH
               IWRK = ILOWT (ICRS)
               XWRK = XDONT (IWRK)
               Do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  If (XWRK < XDONT(IRNGT(IDCR))) Then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                     Exit
                  End If
               End Do
               IRNGT (IDCR+1) = IWRK
               ILOW = ILOW + 1
            End Do
!
            XWRK1 = XDONT (IRNGT(INTH))
            ILOW = 2*INTH - JLOW
            Do ICRS = INTH + 1, JLOW
               If (XDONT(ILOWT (ICRS)) < XWRK1) Then
                  XWRK = XDONT (ILOWT (ICRS))
                  Do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     If (XWRK >= XDONT(IRNGT(IDCR))) Exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  End Do
                  IRNGT (IDCR+1) = ILOWT (ICRS)
                  XWRK1 = XDONT (IRNGT(INTH))
               End If
               ILOW = ILOW + 1
            End Do
!
            INDNTH = IRNGT(INTH)
            Return
!
!
         Case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) Then
               IWRK = ILOWT (1)
               ILOWT (1) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
            End If
            If (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) Then
               IWRK = ILOWT (IFIN)
               ILOWT (IFIN) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
               If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) Then
                  IWRK = ILOWT (1)
                  ILOWT (1) = ILOWT (IMIL)
                  ILOWT (IMIL) = IWRK
               End If
            End If
            If (IFIN <= 3) Exit
!
            XPIV = XDONT (ILOWT(1)) + REAL(INTH)/REAL(JLOW+INTH) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))

!
!  One takes values > XPIV to IHIGT
!
            JHIG = 0
            JLOW = 0
!
            If (XDONT(ILOWT(IFIN)) > XPIV) Then
               ICRS = 0
               Do
                  ICRS = ICRS + 1
                  If (XDONT(ILOWT(ICRS)) > XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     If (ICRS >= IFIN) Exit
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     If (JLOW >= INTH) Exit
                  End If
               End Do
!
               If (ICRS < IFIN) Then
                  Do
                     ICRS = ICRS + 1
                     If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ILOWT (ICRS)
                     Else
                        If (ICRS >= IFIN) Exit
                     End If
                  End Do
               End If
            Else
               Do ICRS = 1, IFIN
                  If (XDONT(ILOWT(ICRS)) > XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     If (JLOW >= INTH) Exit
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  End If
               End Do
            End If
!
         End Select
!
      End Do
!
!  Now, we only need to find maximum of the 1:INTH set
!

      IWRK1 = ILOWT (1)
      XWRK1 =  XDONT (IWRK1)
      Do ICRS = 1+1, INTH
         IWRK = ILOWT (ICRS)
         XWRK = XDONT (IWRK)
         If (XWRK > XWRK1) Then
            XWRK1 = XWRK
            IWRK1 = IWRK
         End If
      End Do
      INDNTH = IWRK1
      Return
!
!
End Function R_indnth
Function I_indnth (XDONT, NORD) Result (INDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! __________________________________________________________
      Integer, Dimension (:), Intent (In)  :: XDONT
      Integer :: INDNTH
      Integer, Intent (In) :: NORD
! __________________________________________________________
      Integer :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      Integer, Dimension (NORD) :: IRNGT
      Integer, Dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      Integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      Integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      Integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = NORD
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      If (NDON < 2) Then
         If (INTH == 1) INDNTH = 1
         Return
      End If
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
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
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         Return
      End If
!
      If (XDONT(3) < XDONT(IHIGT(1))) Then
         IHIGT (2) = IHIGT (1)
         If (XDONT(3) < XDONT(ILOWT(1))) Then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = 3
         Else
            IHIGT (1) = 3
         End If
      Else
         IHIGT (2) = 3
      End If
!
      If (NDON < 4) Then
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         If (INTH == 3) INDNTH = IHIGT (2)
         Return
      End If
!
      If (XDONT(NDON) < XDONT(IHIGT(1))) Then
         IHIGT (3) = IHIGT (2)
         IHIGT (2) = IHIGT (1)
         If (XDONT(NDON) < XDONT(ILOWT(1))) Then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = NDON
         Else
            IHIGT (1) = NDON
         End If
      Else
         IHIGT (3) = NDON
      End If
!
      If (NDON < 5) Then
         If (INTH == 1) INDNTH = ILOWT (1)
         If (INTH == 2) INDNTH = IHIGT (1)
         If (INTH == 3) INDNTH = IHIGT (2)
         If (INTH == 4) INDNTH = IHIGT (3)
         Return
      End If
!

      JLOW = 1
      JHIG = 3
      XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                   (XDONT(IHIGT(3))-XDONT(ILOWT(1)))
      If (XPIV >= XDONT(IHIGT(1))) Then
         XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                      (XDONT(IHIGT(2))-XDONT(ILOWT(1)))
         If (XPIV >= XDONT(IHIGT(1))) &
             XPIV = XDONT (ILOWT(1)) + REAL (2*INTH) / REAL (NDON+INTH) * &
                                          (XDONT(IHIGT(1))-XDONT(ILOWT(1)))
      End If
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the IHIGT array as soon as we have more
!  than enough values in ILOWT.
!
!
      If (XDONT(NDON) > XPIV) Then
         ICRS = 3
         Do
            ICRS = ICRS + 1
            If (XDONT(ICRS) > XPIV) Then
               If (ICRS >= NDON) Exit
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            Else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               If (JLOW >= INTH) Exit
            End If
         End Do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) <= XPIV) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               Else If (ICRS >= NDON) Then
                  Exit
               End If
            End Do
         End If
!
!
      Else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         Do ICRS = 4, NDON - 1
            If (XDONT(ICRS) > XPIV) Then
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            Else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               If (JLOW >= INTH) Exit
            End If
         End Do
!
         If (ICRS < NDON-1) Then
            Do
               ICRS = ICRS + 1
               If (XDONT(ICRS) <= XPIV) Then
                  If (ICRS >= NDON) Exit
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               End If
            End Do
         End If
      End If
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      Do
         If (JLM2 == JLOW .And. JHM2 == JHIG) Then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             If (INTH > JLOW) Then
                XMIN = XDONT (IHIGT(1))
                IHIG = 1
                Do ICRS = 2, JHIG
                   If (XDONT(IHIGT(ICRS)) < XMIN) Then
                      XMIN = XDONT (IHIGT(ICRS))
                      IHIG = ICRS
                   End If
                End Do
!
                JLOW = JLOW + 1
                ILOWT (JLOW) = IHIGT (IHIG)
                IHIGT (IHIG) = IHIGT (JHIG)
                JHIG = JHIG - 1
             Else

                ILOW = ILOWT (1)
                XMAX = XDONT (ILOW)
                Do ICRS = 2, JLOW
                   If (XDONT(ILOWT(ICRS)) > XMAX) Then
                      IWRK = ILOWT (ICRS)
                      XMAX = XDONT (IWRK)
                      ILOWT (ICRS) = ILOW
                      ILOW = IWRK
                   End If
                End Do
                JLOW = JLOW - 1
             End If
         End If
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         Select Case (INTH-JLOW)
         Case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            Select Case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            Case (2)
               If (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) Then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
               Else
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
               End If
               Exit
!
            Case (3)
!
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (3)
               If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                  IHIGT (3) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  End If
               End If
               JHIG = 0
               Do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  ILOWT (ICRS) = IHIGT (JHIG)
               End Do
               JLOW = INTH
               Exit
!
            Case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (IFIN)
               If (XDONT(IWRK2) < XDONT(IWRK1)) Then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               End If
               If (XDONT(IWRK2) > XDONT(IWRK3)) Then
                  IHIGT (IFIN) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  If (XDONT(IWRK2) < XDONT(IHIGT(1))) Then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  End If
               End If
!
               IWRK1 = IHIGT (1)
               JLOW = JLOW + 1
               ILOWT (JLOW) = IWRK1
               XPIV = XDONT (IWRK1) + 0.5 * (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
!
!  One takes values <= pivot to ILOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               Do ICRS = 2, IFIN
                  If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     If (JLOW >= INTH) Exit
                  Else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(IHIGT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  End If
               End Do
            End Select
!
!
         Case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XDONT (IHIGT(1))
            IHIG = 1
            Do ICRS = 2, JHIG
               If (XDONT(IHIGT(ICRS)) < XMIN) Then
                  XMIN = XDONT (IHIGT(ICRS))
                  IHIG = ICRS
               End If
            End Do
!
            INDNTH = IHIGT (IHIG)
            Return
!
!
         Case (0)
!
!  Low part is exactly what we want
!
            Exit
!
!
         Case (-5:-1)
!
!  Only few values too many in low part
!
            IRNGT (1) = ILOWT (1)
            ILOW = 1 + INTH - JLOW
            Do ICRS = 2, INTH
               IWRK = ILOWT (ICRS)
               XWRK = XDONT (IWRK)
               Do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  If (XWRK < XDONT(IRNGT(IDCR))) Then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  Else
                     Exit
                  End If
               End Do
               IRNGT (IDCR+1) = IWRK
               ILOW = ILOW + 1
            End Do
!
            XWRK1 = XDONT (IRNGT(INTH))
            ILOW = 2*INTH - JLOW
            Do ICRS = INTH + 1, JLOW
               If (XDONT(ILOWT (ICRS)) < XWRK1) Then
                  XWRK = XDONT (ILOWT (ICRS))
                  Do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     If (XWRK >= XDONT(IRNGT(IDCR))) Exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  End Do
                  IRNGT (IDCR+1) = ILOWT (ICRS)
                  XWRK1 = XDONT (IRNGT(INTH))
               End If
               ILOW = ILOW + 1
            End Do
!
            INDNTH = IRNGT(INTH)
            Return
!
!
         Case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) Then
               IWRK = ILOWT (1)
               ILOWT (1) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
            End If
            If (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) Then
               IWRK = ILOWT (IFIN)
               ILOWT (IFIN) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
               If (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) Then
                  IWRK = ILOWT (1)
                  ILOWT (1) = ILOWT (IMIL)
                  ILOWT (IMIL) = IWRK
               End If
            End If
            If (IFIN <= 3) Exit
!
            XPIV = XDONT (ILOWT(1)) + REAL(INTH)/REAL(JLOW+INTH) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))

!
!  One takes values > XPIV to IHIGT
!
            JHIG = 0
            JLOW = 0
!
            If (XDONT(ILOWT(IFIN)) > XPIV) Then
               ICRS = 0
               Do
                  ICRS = ICRS + 1
                  If (XDONT(ILOWT(ICRS)) > XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     If (ICRS >= IFIN) Exit
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     If (JLOW >= INTH) Exit
                  End If
               End Do
!
               If (ICRS < IFIN) Then
                  Do
                     ICRS = ICRS + 1
                     If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ILOWT (ICRS)
                     Else
                        If (ICRS >= IFIN) Exit
                     End If
                  End Do
               End If
            Else
               Do ICRS = 1, IFIN
                  If (XDONT(ILOWT(ICRS)) > XPIV) Then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  Else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     If (JLOW >= INTH) Exit
                  End If
               End Do
!
               Do ICRS = ICRS + 1, IFIN
                  If (XDONT(ILOWT(ICRS)) <= XPIV) Then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  End If
               End Do
            End If
!
         End Select
!
      End Do
!
!  Now, we only need to find maximum of the 1:INTH set
!

      IWRK1 = ILOWT (1)
      XWRK1 =  XDONT (IWRK1)
      Do ICRS = 1+1, INTH
         IWRK = ILOWT (ICRS)
         XWRK = XDONT (IWRK)
         If (XWRK > XWRK1) Then
            XWRK1 = XWRK
            IWRK1 = IWRK
         End If
      End Do
      INDNTH = IWRK1
      Return
!
!
End Function I_indnth
end module m_indnth
