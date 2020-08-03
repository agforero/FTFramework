C$Procedure PLTNRM ( DSK, compute outward normal of plate )
 
      SUBROUTINE PLTNRM ( V1, V2, V3, NORMAL )

C$ Abstract
C
C     Compute an outward normal vector of a triangular plate.
C     The vector does not necessarily have unit length.
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
C     DSK
C     FILES
C     TOPOGRAPHY
C   
C$ Declarations

      IMPLICIT NONE

      DOUBLE PRECISION      V1     ( 3 )
      DOUBLE PRECISION      V2     ( 3 )
      DOUBLE PRECISION      V3     ( 3 )
      DOUBLE PRECISION      NORMAL ( 3 )
      
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     V1,
C     V2,
C     V3         I   Vertices of a plate.
C     NORMAL     O   Plate's outward normal vector.
C
C$ Detailed_Input
C
C     V1,
C     V2,
C     V3             are vertices of a triangular plate.
C                    
C$ Detailed_Output
C
C     NORMAL         is an outward normal vector of the plate defined by
C                    the input vertices. The order of the vertices is
C                    used to determine the choice of normal direction:
C                    the normal vector is
C
C                       ( V2 - V1 ) x ( V3 - V2 )
C
C                    
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) The input plate may be degenerate: it may be a line segment
C        or a point. These are not considered to be erroneous inputs.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine saves computation time by not scaling the output
C     vector to unit length. The caller can scale the vector using
C     the routine VHAT.
C
C$ Examples
C
C
C     1) Compute an upward normal of an equilateral triangle lying 
C        in the X-Y plane and centered at the origin.
C
C        Example code begins here.
C
C
C           PROGRAM EX1
C           IMPLICIT NONE
C
C           DOUBLE PRECISION      NORMAL ( 3 )
C           DOUBLE PRECISION      S
C           DOUBLE PRECISION      V1     ( 3 )
C           DOUBLE PRECISION      V2     ( 3 )
C           DOUBLE PRECISION      V3     ( 3 )
C
C
C           S = SQRT(3.D0)/2
C
C           CALL VPACK (    S,  -0.5D0,  0.D0, V1 )
C           CALL VPACK ( 0.D0,    1.D0,  0.D0, V2 )
C           CALL VPACK (   -S,  -0.5D0,  0.D0, V3 )
C
C
C           CALL PLTNRM ( V1, V2, V3, NORMAL )
C
C           WRITE (*, '(1X,A,3(3E20.12))' ) 'NORMAL = ', NORMAL
C
C           END
C
C        When run on a PC/Linux/gfortran/64-bit platform, the output
C        from this program was:
C
C
C  NORMAL =   0.000000000000E+00  0.000000000000E+00  0.259807621135E+01
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
C-    SPICELIB Version 1.0.0, 26-JAN-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     compute normal vector of triangular plate from vertices
C
C-&

C
C     SPICELIB functions
C

C
C     Local variables
C     
      DOUBLE PRECISION      EDGE1  ( 3 )
      DOUBLE PRECISION      EDGE2  ( 3 )

C
C     This routine is error-free.
C

C
C     Type 2 plate vertices are ordered in the positive
C     (right-handed) sense about the outward normal.
C     
      CALL VSUB ( V2, V1, EDGE1 )
      CALL VSUB ( V3, V2, EDGE2 )

      CALL VCRSS ( EDGE1, EDGE2, NORMAL )

      END
