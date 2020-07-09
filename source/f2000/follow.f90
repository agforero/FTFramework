Program follow
!  Question From Colin Thefleau:
! 
! I have a small problem that my cloudy brain can't solve:
! Say I have a bunch of coordinates (realx(i), realy(i)) that form a circle 
! (or any closed form) if every point is plotted. I have to sort them so that 
! they "ride" the circle in one direction. For example beginning at one point 
! (the highest point for example), and go in the clock drive direction.
! Has someone an idea how this is to be done?  Or where can I find a sample 
! code for inspiration? I am really new to fortran and really can't find a 
! solution.
! --------------------------------------------------------------------------
! The following program is an attempt to answer that question for a
! "reasonable" profile. From the current point, it finds the "nearest"
! point in the set of remaining ones according to some weighted distance,
! weights penalizing the direction that one is coming from.
!
   integer, parameter :: nmax = 200
   real, dimension (nmax) :: xptst, yptst, xtmpt, ytmpt, xrndt
   integer, dimension (nmax) :: irndt
   real :: t, xtmp, ytmp, xunt, yunt, xori, yori, xvec, yvec, wdst, wdst0, &
           xlen, xang, xunt1, yunt1
   integer :: imin, imax, ipnt, inxt, itst
!
!  take a continuous curve and make the order random
!
   call random_number (xrndt)
   call mrgrnk (xrndt, irndt)
!
   do ipnt = 1, nmax
     t = 6.28318 * real (ipnt) / real (nmax)
     xtmpt (ipnt) = (5.+ 2 * cos (4.*t))*cos(t)
     ytmpt (ipnt) = -(5.+ 2 * cos (4.*t))*sin(t)
   enddo
   xptst = xtmpt (irndt)
   yptst = ytmpt (irndt)
!
! Bring starting point (Northmost) to first position
!
   imin = sum (maxloc(yptst))
   xtmp = xptst (1)
   ytmp = yptst (1)
   xptst (1) = xptst (imin)
   yptst (1) = yptst (imin)
   xptst (imin) = xtmp
   yptst (imin) = ytmp
!
! unit vector in the current direction (east)
!
   xunt = 1.
   yunt = 0.
!
! Find next point in line
!
   nextpoint: do inxt = 2, nmax-1
      xori = xptst (inxt-1)
      yori = yptst (inxt-1)
      wdst0 = huge(wdst)
      do itst = inxt, nmax
        xvec = xptst (itst) - xori
        yvec = yptst (itst) - yori
        xlen = sqrt (xvec*xvec+yvec*yvec)
        if (xlen < epsilon(1.0)) then
           imin = itst
           xunt1 = xunt
           yunt1 = xunt
           exit
        endif
!
!  Compute distance, weighted by a cosine function of the angle
!  with the last segment. Weight is 1 when straight ahead,
!  3 when going backwards, 2 if transverse. By using some
!  power of the cosine, one may increase or decrease the pressure
!  to go straight ahead with respect to transverse directions.
!
        xang = acos (0.9999*(xvec*xunt+yvec*yunt)/xlen)
        wdst = xlen * (3.0 - 2.0*cos(0.5*xang))
!
!  Retain minimum distance
!
        if (wdst <= wdst0) then
           wdst0 = wdst
           imin = itst
           xunt1 = xvec / xlen
           yunt1 = yvec / xlen
        endif
      enddo
!
!  Exchange retained point with current one
!
      xtmp = xptst (inxt)
      ytmp = yptst (inxt)
      xptst (inxt) = xptst (imin)
      yptst (inxt) = yptst (imin)
      xptst (imin) = xtmp
      yptst (imin) = ytmp
      xunt = xunt1
      yunt = yunt1
   enddo nextpoint 
!
! Output
!
   imax = sum (maxloc(ytmpt))
   do ipnt = 1, nmax
      write (*,*) ipnt,xptst (ipnt), yptst(ipnt), xtmpt (imax), ytmpt (imax)
      imax = mod (imax, nmax) + 1
   enddo 
contains
Subroutine MRGRNK (XVALT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
      Real, Dimension (:), Intent (In) :: XVALT
      Integer, Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      Real (Kind(XVALT)) :: XVALA, XVALB
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XVALT(IIND-1) <= XVALT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Mod(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XVALT (JWRKT(IINDA))
            XVALB = XVALT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XVALT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XVALT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
End Subroutine MRGRNK
end Program follow
