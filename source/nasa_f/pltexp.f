C$Procedure PLTEXP ( Plate expander )
 
      SUBROUTINE PLTEXP ( IVERTS, DELTA, OVERTS )     
 
C$ Abstract
C
C     Expand a triangular plate by a specified amount. The expanded
C     plate is co-planar with, and has the same orientation as, the
C     original. The centroids of the two plates coincide.
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
C     DSK
C
C$ Keywords
C
C     MATH 
C     GEOMETRY
C     TOPOGRAPHY
C
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION      IVERTS ( 3, 3 )
      DOUBLE PRECISION      DELTA
      DOUBLE PRECISION      OVERTS ( 3, 3 )
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     IVERTS     I   Vertices of the plate to be expanded.
C     DELTA      I   Fraction by which the plate is to be expanded.
C     OVERTS     O   Vertices of the expanded plate.
C
C$ Detailed_Input
C
C     IVERTS     is an array containing three vertices of a triangular
C                plate. Each vertex is a three-dimensional vector. The
C                elements
C
C                   IVERTS(J,I), J = 1, 3
C
C                are, respectively, the X, Y, and Z components of the
C                Ith vertex.
C
C
C     DELTA      is a fraction by which the plate is to be scaled.
C                Scaling is done so that the scaled plate has the
C                following properties:
C
C                   -  it is co-planar with the input plate
C
C                   -  its centroid coincides with that of the input
C                      plate
C
C                   -  its sides remain parallel to the corresponding
C                      sides of the input plate
C
C                   -  the distance of each vertex from the centroid is
C                      (1+DELTA) times the corresponding distance for
C                      the input plate
C
C$ Detailed_Output
C
C     OVERTS     is an array containing three vertices of the triangular
C                plate resulting from scaling the input plate.
C
C                If CTROID is the centroid (the average of the vertices)
C                of the input plate, then the Ith vertex of OVERTS
C
C                   OVERTS(J,I), J = 1, 3
C
C                is equal to 
C
C                   CTROID(J) + (1+DELTA)*( IVERTS(J,I) - CTROID(J) ),
C
C                   J = 1, 3
C
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     Error free.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports "greedy" ray-plate intercept algorithms.
C     Such algorithms attempt to ensure that false negatives---in which
C     an intersection is not found due to round-off error---do not
C     occur. In such an algorithm, the plate of interest is expanded
C     slightly before the intersection test is performed.
C     
C$ Examples
C
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as input
C     (if any), the compiler and supporting libraries, and the machine
C     specific arithmetic implementation.
C     
C     1) Expand an equilateral triangle that lies in the plane 
C        
C           { (x,y,z) : z = 7 }
C
C        Use an expansion fraction of 1.D0; this doubles the size of
C        the plate.
C
C        Example code begins here. 
C
C
C              PROGRAM EX1
C              IMPLICIT NONE
C
C              DOUBLE PRECISION      DELTA
C              DOUBLE PRECISION      IVERTS ( 3, 3 )
C              DOUBLE PRECISION      OVERTS ( 3, 3 )
C              DOUBLE PRECISION      S
C
C              INTEGER               I
C
C              S = SQRT(3.D0)/2
C
C              CALL VPACK (    S,  -0.5D0,  7.D0, IVERTS(1,1) )
C              CALL VPACK ( 0.D0,    1.D0,  7.D0, IVERTS(1,2) )
C              CALL VPACK (   -S,  -0.5D0,  7.D0, IVERTS(1,3) )
C
C              DELTA = 1.D0
C
C              CALL PLTEXP ( IVERTS, DELTA, OVERTS )
C
C              WRITE (*,*) ' '
C              WRITE (*,*) 'Vertices of input plate: '
C
C              WRITE (*, '(1X,A,3(3E20.12))' ) ' I1 = ',
C             .          (IVERTS(I,1), I = 1, 3)
C              WRITE (*, '(1X,A,3(3E20.12))' ) ' I2 = ',
C             .          (IVERTS(I,2), I = 1, 3)
C              WRITE (*, '(1X,A,3(3E20.12))' ) ' I3 = ',
C             .          (IVERTS(I,3), I = 1, 3)
C
C              WRITE (*,*) ' '
C              WRITE (*,*) 'Vertices of output plate: '
C
C              WRITE (*, '(1X,A,3(3E20.12))' ) ' O1 = ',
C             .          (OVERTS(I,1), I = 1, 3)
C              WRITE (*, '(1X,A,3(3E20.12))' ) ' O2 = ',
C             .          (OVERTS(I,2), I = 1, 3)
C              WRITE (*, '(1X,A,3(3E20.12))' ) ' O3 = ',
C             .          (OVERTS(I,3), I = 1, 3)
C              WRITE (*,*) ' '
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran/64-bit
C     platform, the output was:
C
C
C     Vertices of input plate:
C      I1 =   0.866025403784E+00 -0.500000000000E+00  0.700000000000E+01
C      I2 =   0.000000000000E+00  0.100000000000E+01  0.700000000000E+01
C      I3 =  -0.866025403784E+00 -0.500000000000E+00  0.700000000000E+01
C
C     Vertices of output plate:
C      O1 =   0.173205080757E+01 -0.100000000000E+01  0.700000000000E+01
C      O2 =   0.000000000000E+00  0.200000000000E+01  0.700000000000E+01
C      O3 =  -0.173205080757E+01 -0.100000000000E+01  0.700000000000E+01
C
C
C     Note that the height of the plate is unchanged, but the vectors
C     from the centroid to the vertices have doubled in length.
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
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 29-FEB-2016 (NJB)
C           
C        Based on original version 28-MAY-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     expand triangular plate
C
C-&
 
C
C     Local variables
C
      DOUBLE PRECISION      SCLCTR ( 3 )
      DOUBLE PRECISION      D
      DOUBLE PRECISION      S

C
C     Compute the centroid of the input vertices. Scale the centroid
C     by DELTA, since we'll only use the scaled form.
C
C     Unroll all loops to avoid loop overhead.
C
      D = DELTA / 3.D0

      SCLCTR(1) = D * ( IVERTS(1,1) + IVERTS(1,2) + IVERTS(1,3) ) 
      SCLCTR(2) = D * ( IVERTS(2,1) + IVERTS(2,2) + IVERTS(2,3) ) 
      SCLCTR(3) = D * ( IVERTS(3,1) + IVERTS(3,2) + IVERTS(3,3) ) 
      
C
C     Compute the offsets of the vertices from the centroid CTROID;
C     scale each offset by (1+DELTA). The Ith expanded vertex is
C
C        CTROID + (1+DELTA) * ( IVERTS(*,I) - CTROID )
C
C     which can be re-written as
C
C        ( (1+DELTA) * IVERTS(*,I) )  -  ( DELTA * CTROID )
C
C     or
C
C        ( (1+DELTA) * IVERTS(*,I) )  -  SCLCTR
C
C
C
      S = 1.D0 + DELTA


      OVERTS(1,1) = ( S * IVERTS(1,1) ) - SCLCTR(1) 
      OVERTS(2,1) = ( S * IVERTS(2,1) ) - SCLCTR(2) 
      OVERTS(3,1) = ( S * IVERTS(3,1) ) - SCLCTR(3) 

      OVERTS(1,2) = ( S * IVERTS(1,2) ) - SCLCTR(1) 
      OVERTS(2,2) = ( S * IVERTS(2,2) ) - SCLCTR(2) 
      OVERTS(3,2) = ( S * IVERTS(3,2) ) - SCLCTR(3) 

      OVERTS(1,3) = ( S * IVERTS(1,3) ) - SCLCTR(1) 
      OVERTS(2,3) = ( S * IVERTS(2,3) ) - SCLCTR(2) 
      OVERTS(3,3) = ( S * IVERTS(3,3) ) - SCLCTR(3) 

      END
