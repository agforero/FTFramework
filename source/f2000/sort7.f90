Program sort7
!
! This program is used to compare 2 algorithms for sorting
! every successive subset of 7 elements in an array.
! It shows how to call a subroutine with the correct interface
! block in the main program. Interfaces for ORDERPACK routines
! can be found in file interfaces.f90
!
      Interface
         Subroutine oldsub7 (xdont, xwrkt)
            Real, Dimension (:), Intent (In) :: xdont
            Real, Dimension (:), Intent (Out) :: xwrkt
         End Subroutine oldsub7
         Subroutine newsub7 (xdont, xwrkt)
            Real, Dimension (:), Intent (In) :: xdont
            Real, Dimension (:), Intent (Out) :: xwrkt
         End Subroutine newsub7
      End Interface
!
      Real, Dimension (35280) :: xvalt
      Real, Dimension (35280) :: xwrkt
      Real, Dimension (7) :: xval7t
      Real, Dimension (6) :: xval6t
      Real, Dimension (5) :: xval5t
      Real, Dimension (4) :: xval4t
      Real, Dimension (3) :: xval3t
      Real, Dimension (2) :: xval2t
      Real, Dimension (1) :: xval1t
      Real :: x1, x2, x3, x4, x5, x6, x7
      Integer :: i1, i2, i3, i4, i5, i6, ival
!
! generate all 5040 possible permutations of elements 1,...,7
!
      xval7t = (/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 /)
      ival = 0
      Do i1 = 1, 7
         x1 = xval7t (i1)
         xval6t = pack (xval7t, (xval7t /= x1))
         Do i2 = 1, 6
            x2 = xval6t (i2)
            xval5t = pack (xval6t, (xval6t /= x2))
            Do i3 = 1, 5
               x3 = xval5t (i3)
               xval4t = pack (xval5t, (xval5t /= x3))
               Do i4 = 1, 4
                  x4 = xval4t (i4)
                  xval3t = pack (xval4t, (xval4t /= x4))
                  Do i5 = 1, 3
                     x5 = xval3t (i5)
                     xval2t = pack (xval3t, (xval3t /= x5))
                     Do i6 = 1, 2
                        x6 = xval2t (i6)
                        xval1t = pack (xval2t, (xval2t /= x6))
                        x7 = xval1t (1)
                        xvalt (ival+1:ival+7) = (/ x1, x2, x3, x4, x5, &
                       & x6, x7 /)
                        ival = ival + 7
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
      Call newsub7 (xvalt, xwrkt)
      Do i1 = 1, 35280, 7
         If (any(xwrkt(i1:i1+6) /= xval7t)) Then
            Write (*,*) "newsub KO", i1, xwrkt (i1:i1+6)
            Exit
         End If
      End Do
      Call oldsub7 (xvalt, xwrkt)
      Do i1 = 1, 35280, 7
         If (any(xwrkt(i1:i1+6) /= xval7t)) Then
            Write (*,*) "oldsub KO", i1, xwrkt (i1:i1+6)
         End If
      End Do
End Program sort7
Subroutine oldsub7 (xdont, xwrkt)
!
! This subroutine is based on insertion sort, with a first
! pass of selection to bring the smallest value in first
! location, thus avoiding a test in the insertion phase
!
      Real, Dimension (:), Intent (In) :: xdont
      Real, Dimension (:), Intent (Out) :: xwrkt
      Real :: xwrk, xwrk1
      Integer :: ndon, ideb, iwrk, idcr, icrs
!
      ndon = size (xdont)
!
      Do ideb = 1, ndon - 6, 7
         idcr = ideb + 6
         If (xdont(ideb) < xdont(idcr)) Then
            xwrk = xdont (ideb)
            xwrkt (idcr) = xdont (idcr)
         Else
            xwrk = xdont (idcr)
            xwrkt (idcr) = xdont (ideb)
         End If
         Do iwrk = 1, 5
            idcr = idcr - 1
            xwrk1 = xdont (idcr)
            If (xwrk1 < xwrk) Then
               xwrkt (idcr) = xwrk
               xwrk = xwrk1
            Else
               xwrkt (idcr) = xwrk1
            End If
         End Do
         xwrkt (ideb) = xwrk
         Do icrs = ideb + 2, ideb + 6
            xwrk = xwrkt (icrs)
            If (xwrk < xwrkt(icrs-1)) Then
               xwrkt (icrs) = xwrkt (icrs-1)
               idcr = icrs - 1
               xwrk1 = xwrkt (idcr-1)
               Do
                  If (xwrk >= xwrk1) Exit
                  xwrkt (idcr) = xwrk1
                  idcr = idcr - 1
                  xwrk1 = xwrkt (idcr-1)
               End Do
               xwrkt (idcr) = xwrk
            End If
         End Do
      End Do
End Subroutine oldsub7
Subroutine newsub7 (xdont, xwrkt)
      Real, Dimension (:), Intent (In) :: xdont
      Real, Dimension (:), Intent (Out) :: xwrkt
      Real :: xwrk1, xwrk2, xwrk3, xwrk4, xwrk5, xwrk6
      Real :: xwrk, xwrki, xwrks
      Integer :: ndon, ideb, ideba
!
! This routine uses a sort of Shellsort in a first pass
! trying to bring the set to the following form:
! {minimum} {X1 . . } {X4 . .} with X1 <= X4
! then
! {minimum} {X1 X2 . } {X4 X5 .} with X2 <= X5
! then sort each of the subsets, and then
! {minimum} {second minimum X2 X3 } {X4 X5 X6=maximum} with X2 <= X3 and
!                                                           X4 <= X5
!  and lastly merge the two subsets.
!
!
      ndon = size (xdont)
      xwrkt = 0.0
!
      Do ideb = 1, ndon - 6, 7
         xwrk = xdont (ideb)
!
!  First pair
!
         ideba = ideb + 1
         If (xdont(ideba) < xdont(ideba+3)) Then
            xwrks = xdont (ideba+3)
            xwrki = xdont (ideba)
         Else
            xwrki = xdont (ideba+3)
            xwrks = xdont (ideba)
         End If
         If (xwrki < xwrk) Then
            If (xwrks < xwrk) Then
               xwrk1 = xwrks
               xwrk4 = xwrk
            Else
               xwrk1 = xwrk
               xwrk4 = xwrks
            End If
            xwrk = xwrki
         Else
            xwrk1 = xwrki
            xwrk4 = xwrks
         End If
!
!  Second pair
!
         ideba = ideba + 1
         If (xdont(ideba) < xdont(ideba+3)) Then
            xwrks = xdont (ideba+3)
            xwrki = xdont (ideba)
         Else
            xwrki = xdont (ideba+3)
            xwrks = xdont (ideba)
         End If
         If (xwrki < xwrk) Then
            xwrk2 = xwrk1
            If (xwrks < xwrk) Then
               xwrk1 = xwrks
               xwrks = xwrk
            Else
               xwrk1 = xwrk
            End If
            xwrk = xwrki
         Else
            If (xwrki >= xwrk1) Then
               xwrk2 = xwrki
            Else
               xwrk2 = xwrk1
               xwrk1 = xwrki
            End If
         End If
         If (xwrks >= xwrk4) Then
            xwrk5 = xwrks
         Else
            xwrk5 = xwrk4
            xwrk4 = xwrks
         End If
!
!  Third pair
!
         ideba = ideba + 1
         If (xdont(ideba) < xdont(ideba+3)) Then
            xwrks = xdont (ideba+3)
            xwrki = xdont (ideba)
         Else
            xwrki = xdont (ideba+3)
            xwrks = xdont (ideba)
         End If
         If (xwrki <= xwrk) Then
            xwrk3 = xwrk2
            xwrk2 = xwrk1
            If (xwrks < xwrk) Then
               xwrk1 = xwrks
               xwrks = xwrk
            Else
               xwrk1 = xwrk
            End If
            xwrk = xwrki
         Else
            If (xwrki >= xwrk2) Then
               xwrk3 = xwrki
            Else
               xwrk3 = xwrk2
               If (xwrki >= xwrk1) Then
                  xwrk2 = xwrki
               Else
                  xwrk2 = xwrk1
                  xwrk1 = xwrki
               End If
            End If
         End If
         If (xwrks >= xwrk5) Then
            xwrk6 = xwrks
         Else
            xwrk6 = xwrk5
            If (xwrks >= xwrk4) Then
               xwrk5 = xwrks
            Else
               xwrk5 = xwrk4
               xwrk4 = xwrks
            End If
         End If
!
! Merge the two subsets into their final location
!
         xwrkt (ideb) = xwrk
         xwrkt (ideb+1) = xwrk1
         If (xwrk4 >= xwrk3) Then
            xwrkt (ideb+2) = xwrk2
            xwrkt (ideb+3) = xwrk3
            xwrkt (ideb+4) = xwrk4
            xwrkt (ideb+5) = xwrk5
         Else
            If (xwrk4 >= xwrk2) Then
               xwrkt (ideb+2) = xwrk2
               xwrkt (ideb+3) = xwrk4
            Else
               xwrkt (ideb+2) = xwrk4
               xwrkt (ideb+3) = xwrk2
            End If
            If (xwrk3 <= xwrk5) Then
               xwrkt (ideb+4) = xwrk3
               xwrkt (ideb+5) = xwrk5
            Else
               xwrkt (ideb+4) = xwrk5
               xwrkt (ideb+5) = xwrk3
            End If
         End If
         xwrkt (ideb+6) = xwrk6
      End Do
End Subroutine newsub7
