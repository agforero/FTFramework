C$Procedure      NPSGPT ( Nearest point on line segment )
 
      SUBROUTINE NPSGPT ( EP1, EP2, POINT, PNEAR, DIST )
 
C$ Abstract
C
C     Find the nearest point on a line segment to a given point.
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
C     GEOMETRY
C     MATH
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      EP1   ( 3 )
      DOUBLE PRECISION      EP2   ( 3 )
      DOUBLE PRECISION      POINT ( 3 )
      DOUBLE PRECISION      PNEAR ( 3 )
      DOUBLE PRECISION      DIST

 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     EP1,
C     EP2        I   Endpoints of a line segment.
C     POINT      I   A point in 3-dimensional space.
C     PNEAR      O   Nearest point on the line segment to POINT.
C     DIST       O   Distance between PNEAR and POINT.
C
C$ Detailed_Input
C
C     EP1,
C     EP2        are the endpoints of a line segment in 3-dimensional 
C                space. EP1 and EP2 need not be distinct.
C
C     POINT      is an arbitrary point in 3-dimensional space.
C          
C$ Detailed_Output
C
C     PNEAR      is the closest point on the line segment to POINT.
C
C     DIST       is the distance between POINT and PNEAR.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) The input segment is allowed to be degenerate: it may be
C        a single point.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     None.
C
C$ Examples
C
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as input
C     (if any), the compiler and supporting libraries, and the machine
C     specific arithmetic implementation.
C
C     1) Compute the nearest point on a line segment to a given
C        point in a simple case for which the results can easily be
C        checked.
C
C
C        Example code begins here. 
C
C 
C              PROGRAM EX1
C              IMPLICIT NONE
C        C
C        C     Local parameters
C        C
C              CHARACTER*(*)         FMT1
C              PARAMETER           ( FMT1 = '(A,3F13.8)' )
C        C
C        C     Local variables
C        C
C              DOUBLE PRECISION      DIST
C              DOUBLE PRECISION      ENDPT1 ( 3 )
C              DOUBLE PRECISION      ENDPT2 ( 3 )
C              DOUBLE PRECISION      PNEAR  ( 3 )
C              DOUBLE PRECISION      POINT  ( 3 )
C
C        C
C        C     Initialize the line segment's endpoints.
C        C
C              CALL VPACK ( 1.D0, -2.D0, 3.D0, ENDPT1 )
C              CALL VPACK ( 1.D0,  2.D0, 3.D0, ENDPT2 )
C        C
C        C     Set the input point.
C        C
C              CALL VPACK ( 1.D0,  0.D0, 0.D0, POINT )
C        C
C        C     Find the near point on the segment.
C        C
C              CALL NPSGPT ( ENDPT1, ENDPT2, POINT, PNEAR, DIST )
C
C              WRITE (*,*) ' '
C              WRITE (*,FMT1) 'Endpoint 1:  ', ENDPT1
C              WRITE (*,FMT1) 'Endpoint 2:  ', ENDPT2
C              WRITE (*,FMT1) 'Point:       ', POINT
C              WRITE (*,*) ' '
C              WRITE (*,FMT1) 'Near point:  ', PNEAR
C              WRITE (*,FMT1) 'Distance:    ', DIST
C              WRITE (*,*) ' '
C
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran/64-bit
C     platform, the output was:
C
C
C        Endpoint 1:     1.00000000  -2.00000000   3.00000000
C        Endpoint 2:     1.00000000   2.00000000   3.00000000
C        Point:          1.00000000   0.00000000   0.00000000
C
C        Near point:     1.00000000   0.00000000   3.00000000
C        Distance:       3.00000000
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
C     N.J. Bachman   (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 02-FEB-MAR-2016 (NJB)
C
C        Updated from DSKLIB Version 1.0.0, 20-MAR-2015 (NJB)
C
C-&
 
C$ Index_Entries
C
C     nearest point on line segment
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      VDIST
      DOUBLE PRECISION      VDOT

      LOGICAL               FAILED
      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     Local variables
C
      DOUBLE PRECISION      SEG    ( 3 )
      DOUBLE PRECISION      SEGDOT
      DOUBLE PRECISION      LNEAR  ( 3 )
      DOUBLE PRECISION      OFFSET ( 3 )
      DOUBLE PRECISION      OFFDOT

C
C     Use discovery check-in.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Find a direction vector defined by the endpoints.
C     
      CALL VSUB ( EP2, EP1, SEG )

      IF ( VZERO(SEG) ) THEN
C
C        The endpoints coincide, and both coincide with the 
C        near point.
C        
         CALL VEQU ( EP1, PNEAR )
         
         DIST = VDIST( EP1, POINT )

         RETURN

      END IF

C
C     Find the nearest point to POINT on the line defined by
C     EP1 and SEG.
C      
      CALL NPLNPT ( EP1, SEG, POINT, LNEAR, DIST )

      IF ( FAILED() ) THEN
         RETURN
      END IF

C
C     Determine whether LNEAR is on the segment, "before" EP1, or
C     "after" EP2, where SEG points in the "increasing" direction.
C
      CALL VSUB ( LNEAR, EP1, OFFSET )

      OFFDOT = VDOT( OFFSET, SEG )

      IF ( OFFDOT .LT. 0.D0 ) THEN
C
C        The nearest point on the line precedes the first endpoint.
C        The closest point on the segment is the first endpoint.
C
         CALL VEQU ( EP1, PNEAR )

         DIST = VDIST ( EP1, POINT )         
         
      ELSE
C
C        See whether OFFSET is past the second endpoint. Compare
C        the dot product of OFFSET with SEG to that of SEG with
C        itself, since SEG is the offset of EP2 from EP1.
C
         SEGDOT = VDOT( SEG, SEG )

         IF ( OFFDOT .GT. SEGDOT ) THEN
C
C           The nearest point on the line follows the last endpoint.
C           The closest point on the segment is the last endpoint.
C
            CALL VEQU ( EP2, PNEAR )

            DIST = VDIST ( EP2, POINT )

         ELSE
C
C           The near point is on the segment. LNEAR is actually the
C           solution.
C
            CALL VEQU ( LNEAR, PNEAR )
            
C
C           DIST was correctly set by the call to NPLNPT.
C            
         END IF

      END IF

      RETURN
      END





