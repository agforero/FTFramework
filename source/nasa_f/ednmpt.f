C$Procedure EDNMPT ( Ellipsoid normal vector to surface point )
 
      SUBROUTINE EDNMPT ( A, B, C, NORMAL, POINT )
 
C$ Abstract
C
C     Return the unique point on an ellipsoid's surface where the
C     outward normal direction is a given vector.
C
C$ Disclaimer
C
C     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
C     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
C     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
C     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
C     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
C     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
C     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
C     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
C     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
C     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.
C
C     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
C     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
C     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
C     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
C     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
C     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
C
C     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
C     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
C     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
C     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.
C
C$ Required_Reading
C
C     None.
C
C$ Keywords
C
C     ELLIPSOID
C     GEOMETRY
C     NORMAL
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      A
      DOUBLE PRECISION      B
      DOUBLE PRECISION      C
      DOUBLE PRECISION      NORMAL ( 3 )
      DOUBLE PRECISION      POINT  ( 3 )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     A          I   Length of the ellipsoid semi-axis along the x-axis.
C     B          I   Length of the ellipsoid semi-axis along the y-axis.
C     C          I   Length of the ellipsoid semi-axis along the z-axis.
C     NORMAL     I   Outward normal direction.
C     POINT      O   Point where outward normal is parallel to NORMAL.
C
C$ Detailed_Input
C
C     A          is the length of the semi-axis of the ellipsoid
C                that is parallel to the x-axis of the body-fixed
C                coordinate system.
C
C     B          is the length of the semi-axis of the ellipsoid
C                that is parallel to the y-axis of the body-fixed
C                coordinate system.
C
C     C          is the length of the semi-axis of the ellipsoid
C                that is parallel to the z-axis of the body-fixed
C                coordinate system.
C
C     NORMAL     is a non-zero vector. The unique point on the
C                ellipsoid at which NORMAL is an outward normal vector
C                is sought.
C
C$ Detailed_Output
C
C     POINT      is the unique point on the ellipsoid at which NORMAL
C                is an outward normal vector.
C
C                POINT is a 3-vector giving the bodyfixed coordinates
C                of a point on the ellipsoid. In bodyfixed coordinates,
C                the semi-axes of the ellipsoid are aligned with the x,
C                y, and z-axes of the coordinate system.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If any of the semi-axis lengths is non-positive, the error
C         SPICE(BADAXISLENGTH) will be signaled.
C
C     2)  If any of the semi-axis lengths underflows to zero when
C         divided by the largest semi-axis length, the error
C         SPICE(AXISUNDERFLOW) will be signaled.
C
C     3)  If NORMAL is the zero vector, the error SPICE(ZEROVECTOR)
C         will be signaled.
C
C     4)  If the input pass the above checks but lead to a 
C         divide-by-zero error or to an computing an invalid argument
C         of a fractional exponential expression, the error
C         SPICE(DEGENERATECASE) will be signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine can be used to determine the distance between an
C     ellipsoid and a non-intersecting plane. This distance computation
C     supports computation of terminator points on an ellipsoid.
C
C$ Examples
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine 
C     specific arithmetic implementation. 
C
C     1) Choose a triaxial ellipsoid with three unequal semi-axis
C        lengths. Pick several vectors; find the points on the
C        ellipsoid where the respective outward normals are parallel to
C        those vectors. 
C
C        Check the results: at each point, a computed outward normal
C        vector should have very small angular separation from the
C        input vector. Also, the point should be on the surface of the
C        ellipsoid. The ellipsoid can be thought of as a level surface
C        of the function
C
C                             2        2         2
C           f(x, y, z) = (x/a)  + (y/b)  +  (z/c)
C
C        where a, b, c are the semi-axis lengths of the ellipsoid.
C        Specifically, the ellipsoid is the set
C
C           { (x, y, z) : f(x, y, z)  =  1 }
C
C        We can evaluate f at a point to determine whether that point
C        is close to the ellipsoid's surface.
C
C
C        Example code begins here.
C 
C
C              PROGRAM EX1
C              IMPLICIT NONE
C        C
C        C     SPICELIB functions
C        C
C              DOUBLE PRECISION      VSEP
C        C
C        C     Local parameters
C        C
C              CHARACTER*(*)         FMT1
C              PARAMETER           ( FMT1 = '(A,3F15.9)' )
C        C
C        C     Local variables
C        C
C              DOUBLE PRECISION      A
C              DOUBLE PRECISION      B
C              DOUBLE PRECISION      C
C              DOUBLE PRECISION      NORMAL ( 3 )
C              DOUBLE PRECISION      POINT  ( 3 )
C              DOUBLE PRECISION      XNORML ( 3 )
C        C
C        C     Initialize the ellipsoid semi-axes.
C        C
C              A = 10.D0
C              B =  5.D0
C              C =  2.D0
C        C
C        C     Pick several vectors; find the points
C        C     on the ellipsoid where the respective
C        C     outward normals are parallel to those
C        C     vectors; check the results.
C        C
C              CALL VPACK  ( 0.D0, 0.D0, 3.D0, XNORML )
C              CALL EDNMPT ( A,    B,    C,    XNORML, POINT  )
C              CALL SURFNM ( A,    B,    C,    POINT,  NORMAL )
C
C              WRITE (*,*   ) ' '
C              WRITE (*,FMT1) 'Semi-axis lengths:   ', A, B, C
C              WRITE (*,FMT1) 'Input vector:        ', XNORML
C              WRITE (*,FMT1) 'Point:               ', POINT
C              WRITE (*,FMT1) 'Outward normal:      ', NORMAL
C              WRITE (*,FMT1) 'Angular error (rad): ', VSEP(NORMAL,
C             .                                          XNORML )
C              WRITE (*,FMT1) 'Off-surface error:   ',
C             .                 (POINT(1)/A)**2 + (POINT(2)/B)**2
C             .               + (POINT(3)/C)**2 - 1
C              WRITE (*,*) ' '
C
C
C              CALL VPACK  ( 15.D0, -7.D0, 3.D0, XNORML )
C              CALL EDNMPT ( A,      B,    C,    XNORML, POINT  )
C              CALL SURFNM ( A,      B,    C,    POINT,  NORMAL )
C
C              WRITE (*,FMT1) 'Semi-axis lengths:   ', A, B, C
C              WRITE (*,FMT1) 'Input vector:        ', XNORML
C              WRITE (*,FMT1) 'Point:               ', POINT
C              WRITE (*,FMT1) 'Outward normal:      ', NORMAL
C              WRITE (*,FMT1) 'Angular error (rad): ', VSEP(NORMAL,
C             .                                             XNORML )
C              WRITE (*,FMT1) 'Off-surface error:   ',
C             .                 (POINT(1)/A)**2 + (POINT(2)/B)**2
C             .               + (POINT(3)/C)**2 - 1
C              WRITE (*,*) ' '
C
C              CALL VPACK  ( 15.D0, -7.D0, 3.D0, XNORML )
C              CALL EDNMPT ( A,      B,    C,    XNORML, POINT  )
C              CALL SURFNM ( A,      B,    C,    POINT,  NORMAL )
C
C              WRITE (*,FMT1) 'Semi-axis lengths:   ', A, B, C
C              WRITE (*,FMT1) 'Input vector:        ', XNORML
C              WRITE (*,FMT1) 'Point:               ', POINT
C              WRITE (*,FMT1) 'Outward normal:      ', NORMAL
C              WRITE (*,FMT1) 'Angular error (rad): ', VSEP(NORMAL,
C             .                                             XNORML )
C              WRITE (*,FMT1) 'Off-surface error:   ',
C             .                 (POINT(1)/A)**2 + (POINT(2)/B)**2
C             .               + (POINT(3)/C)**2 - 1
C              WRITE (*,*) ' '
C
C              CALL VPACK  ( A/2, B/2,  C/2, XNORML )
C              CALL EDNMPT ( A,   B,    C,   XNORML, POINT  )
C              CALL SURFNM ( A,   B,    C,   POINT,  NORMAL )
C
C              WRITE (*,FMT1) 'Semi-axis lengths:   ', A, B, C
C              WRITE (*,FMT1) 'Input vector:        ', XNORML
C              WRITE (*,FMT1) 'Point:               ', POINT
C              WRITE (*,FMT1) 'Outward normal:      ', NORMAL
C              WRITE (*,FMT1) 'Angular error (rad): ', VSEP(NORMAL,
C             .                                             XNORML )
C              WRITE (*,FMT1) 'Off-surface error:   ',
C             .                 (POINT(1)/A)**2 + (POINT(2)/B)**2
C             .               + (POINT(3)/C)**2 - 1
C              WRITE (*,*) ' '
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was:
C 
C
C     Semi-axis lengths:      10.000000000    5.000000000    2.000000000
C     Input vector:            0.000000000    0.000000000    3.000000000
C     Point:                   0.000000000    0.000000000    2.000000000
C     Outward normal:          0.000000000    0.000000000    1.000000000
C     Angular error (rad):     0.000000000
C     Off-surface error:       0.000000000
C
C     Semi-axis lengths:      10.000000000    5.000000000    2.000000000
C     Input vector:           15.000000000   -7.000000000    3.000000000
C     Point:                   9.731032027   -1.135287070    0.077848256
C     Outward normal:          0.891657447   -0.416106809    0.178331489
C     Angular error (rad):     0.000000000
C     Off-surface error:       0.000000000
C
C     Semi-axis lengths:      10.000000000    5.000000000    2.000000000
C     Input vector:           15.000000000   -7.000000000    3.000000000
C     Point:                   9.731032027   -1.135287070    0.077848256
C     Outward normal:          0.891657447   -0.416106809    0.178331489
C     Angular error (rad):     0.000000000
C     Off-surface error:       0.000000000
C
C     Semi-axis lengths:      10.000000000    5.000000000    2.000000000
C     Input vector:            5.000000000    2.500000000    1.000000000
C     Point:                   9.694128639    1.211766080    0.077553029
C     Outward normal:          0.880450906    0.440225453    0.176090181
C     Angular error (rad):     0.000000000
C     Off-surface error:       0.000000000
C
C
C$ Restrictions
C
C     None.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     W.L. Taber      (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 17-MAY-2016 (NJB) (WLT)
C
C-&
 
C$ Index_Entries
C
C     point on an ellipsoid having given surface normal 
C
C-&
 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      TOUCHD

      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local variables
C
      DOUBLE PRECISION      ARG
      DOUBLE PRECISION      LAMBDA
      DOUBLE PRECISION      NA2
      DOUBLE PRECISION      NB2
      DOUBLE PRECISION      NC2
      DOUBLE PRECISION      SA
      DOUBLE PRECISION      SB
      DOUBLE PRECISION      SC
      DOUBLE PRECISION      SCALE

C
C     Use discovery check-in. We will test RETURN though, since
C     we could have invalid inputs to our arithmetic computations
C     if a SPICE error condition exists upon entry.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Perform some preliminary checks.
C     
C     We need a non-degenerate ellipsoid to start with.
C
      IF (      ( A .LE. 0.D0 ) 
     .     .OR. ( B .LE. 0.D0 ) 
     .     .OR. ( C .LE. 0.D0 )  ) THEN

         CALL CHKIN  ( 'EDNMPT'                               )
         CALL SETMSG ( 'All ellipsoid semi-axis lengths must '
     .   //            'be strictly positive. Lengths were: '
     .   //            'A = #; B = #; C = #'                  )
         CALL ERRDP  ( '#',  A                                )
         CALL ERRDP  ( '#',  B                                )
         CALL ERRDP  ( '#',  C                                )
         CALL SIGERR ( 'SPICE(BADAXISLENGTH)'                 )
         CALL CHKOUT ( 'EDNMPT'                               )
         RETURN

      END IF

C
C     We'll work with scaled copies of the input semi-axis
C     lengths.
C
      SCALE = MAX ( A, B, C )

      SA    = TOUCHD( A / SCALE )
      SB    = TOUCHD( B / SCALE )
      SC    = TOUCHD( C / SCALE )
C
C     If any of the scaled-semi axes underflowed to zero,
C     we can't continue.
C
      IF (      ( SA .LE. 0.D0 ) 
     .     .OR. ( SB .LE. 0.D0 ) 
     .     .OR. ( SC .LE. 0.D0 )  ) THEN

         CALL CHKIN  ( 'EDNMPT'                                   )
         CALL SETMSG ( 'Scaled semi-axis lengths must be '
     .   //            'strictly positive. Scaled lengths were: '
     .   //            'SA = #; SB = #; SC = #'                   )
         CALL ERRDP  ( '#',  SA                                   )
         CALL ERRDP  ( '#',  SB                                   )
         CALL ERRDP  ( '#',  SC                                   )
         CALL SIGERR ( 'SPICE(AXISUNDERFLOW)'                     )
         CALL CHKOUT ( 'EDNMPT'                                   )
         RETURN

      END IF

C
C     The normal vector can't be the zero vector.
C
      IF ( VZERO(NORMAL) ) THEN

         CALL CHKIN  ( 'EDNMPT'                                       )
         CALL SETMSG ( 'The input normal vector was the zero vector. '
     .   //            'There is no solution.'                        )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                            )
         CALL CHKOUT ( 'EDNMPT'                                       )
         RETURN
         
      END IF

C
C     In general, given an ellipsoid E with axes a, b, c,
C     and given a normal vector
C
C        n = ( n1, n2, n3 )
C
C     this problem is equivalent to finding a point (x,y,z)
C     on E such that
C
C             2     2     2
C        ( x/a,  y/b,  z/c ) = lambda * n
C
C     which implies
C                                           2     2     2
C        ( x, y, z )         = lambda ( n1*a, n2*b, n3*c )
C
C
C     Then since the vector on the left side is on the surface
C     of E, we must have
C
C              2       2   2        2   2       2   2
C        lambda  ( ( n1 * a  ) + (n2 * b ) + (n3 + c ) ) = 1
C
C
C     Requiring lambda to be positive, we have
C
C                               2   2        2   2       2   2
C        lambda = 1 / sqrt( ( n1 * a  ) + (n2 * b ) + (n3 + c ) )
C
C
      NA2 = NORMAL(1) * SA * SA
      NB2 = NORMAL(2) * SB * SB
      NC2 = NORMAL(3) * SC * SC

      ARG = TOUCHD(   ( NA2 * NORMAL(1) ) + ( NB2 * NORMAL(2) ) 
     .              + ( NC2 * NORMAL(3) )                        )

      IF ( ARG .LE. 0.D0 ) THEN

         CALL CHKIN  ( 'EDNMPT'                                   )
         CALL SETMSG ( 'Scale factor LAMBDA must be positive, '
     .   //            'but reciprocal of square of LAMBDA is #.' )
         CALL ERRDP  ( '#', ARG                                   )
         CALL SIGERR ( 'SPICE(DEGENERATECASE)'                    )
         CALL CHKOUT ( 'EDNMPT'                                   )
         RETURN
 
      END IF

C
C     Compute LAMBDA as above, and scale it too. This will place
C     POINT on the original ellipsoid.
C
      LAMBDA   = ( ARG ** (-0.5D0) ) * SCALE

      POINT(1) = LAMBDA * NA2
      POINT(2) = LAMBDA * NB2
      POINT(3) = LAMBDA * NC2

      END

      
