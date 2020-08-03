C$Procedure ZZINREC ( DSK, in rectangular element? )

      SUBROUTINE ZZINREC ( P, BOUNDS, MARGIN, EXCLUD, INSIDE )

C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Test a point represented by a set of rectangular coordinates for
C     inclusion in a specified rectangular volume element. The test is
C     performed using margins for the element.
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
C     GEOMETRY
C     INTERSECTION
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      P      ( 3 )
      DOUBLE PRECISION      BOUNDS ( 2, 3 )
      DOUBLE PRECISION      MARGIN
      INTEGER               EXCLUD
      LOGICAL               INSIDE

      INTEGER               NONE
      PARAMETER           ( NONE   = 0 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     P          I   Input point.
C     BOUNDS     I   Coordinate bounds of element.
C     MARGIN     I   Element expansion margin.
C     EXCLUD     I   Index of coordinate to exclude from test.
C     INSIDE     O   Flag indicating whether point is in element.
C     NONE       P   Meaning: exclude nothing.
C
C$ Detailed_Input
C
C     P          is a point expressed in Cartesian coordinates. The
C                point is to be checked to determine whether it is
C                inside the volume element specified by BOUNDS.
C
C
C     BOUNDS     is an 2x3 array containing the bounds of a volume
C                element expressed in rectangular coordinates. BOUNDS
C                defines the volume element used in the comparison. In
C                the element
C
C                   BOUNDS(I,J) 
C
C                J is the coordinate index. J is one of
C
C                   { 1, 2, 3 }
C
C                I is the bound index.
C
C                   I = 1   ->   lower bound
C                   I = 2   ->   upper bound                
C
C
C     MARGIN     is a scale factor used to expand the volume element
C                (box) described by the input bounds. Each side of the
C                element is scaled by MARGIN, and the coordinates of
C                the corners of the element are moved outward by the
C                resulting distance to create an expanded element. The
C                input point is tested for containment within this
C                expanded element.
C
C
C     EXCLUD     is either a coordinate index or the parameter NONE.
C
C                If EXCLUD is set to one of
C
C                   { 1, 2, 3 }
C
C                then the indicated coordinate is excluded from
C                comparison with the corresponding volume element
C                boundaries.
C
C                If EXCLUD is set to NONE, all coordinates are
C                compared.
C
C                Exclusion of coordinates is used in cases where a
C                point is known to be on a level surface of a given
C                coordinate. For example, if a point is on the plane
C                equal to the upper Z bound, the point's Z coordinate
C                need not be used in the comparison and in fact can't
C                be meaningfully compared, due to round-off errors.
C
C$ Detailed_Output
C
C     INSIDE     is a logical flag that is set to .TRUE. if and
C                only if the input coordinates represent a 
C                point inside or on the surface of the volume
C                element, according to the comparisons that are
C                performed.
C
C                The value of INSIDE is not affected by the value
C                of any excluded coordinate.
C
C$ Parameters
C
C     NONE       when used as a value of the input argument EXCLUD,
C                indicates that no coordinates should be excluded from
C                comparison.
C
C$ Exceptions
C
C     1)  If MARGIN is negative, the error SPICE(VALUEOUTOFRANGE)
C         is signaled.
C
C     2) If EXCLUD is less than NONE or greater than 3, the error
C        SPICE(INDEXOUTOFRANGE) is signaled.
C
C     3) If any coordinate upper bound is less than the corresponding
C        lower bound, the error SPICE(BOUNDSOUTOFORDER) is signaled.
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
C     See usage in ZZRYTREC.
C
C$ Restrictions
C
C     This is a private routine. It is meant to be used only by the DSK
C     subsystem.
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
C-    SPICELIB Version 1.0.0, 16-MAY-2016 (NJB) 
C
C        Original version 03-OCT-2014 (NJB)
C
C-&
 
C$ Index_Entries
C
C     test point against rectangular element using margin
C
C-&


C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Element boundary indices:
C     
      INTEGER               LOWER
      PARAMETER           ( LOWER  = 1 )

      INTEGER               UPPER
      PARAMETER           ( UPPER  = 2 )

C
C     Local variables
C
      DOUBLE PRECISION      AMNCOR ( 3 )
      DOUBLE PRECISION      AMXCOR ( 3 )
      DOUBLE PRECISION      DELTA  ( 3 )
      DOUBLE PRECISION      MAXCOR ( 3 )
      DOUBLE PRECISION      MINCOR ( 3 )
      DOUBLE PRECISION      L      ( 3 )
      
      INTEGER               I


C
C     Check-in is discovery style.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Assume the point is outside to start. This allows us
C     to skip setting INSIDE when we find a boundary test
C     failure.
C
      INSIDE = .FALSE.


C
C     Reject negative margins.
C
      IF ( MARGIN .LT. 0.D0 ) THEN

         CALL CHKIN  ( 'ZZINREC'                                )
         CALL SETMSG ( 'Margin must be non-negative but was #.' )
         CALL ERRDP  ( '#', MARGIN                              )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'                 )
         CALL CHKOUT ( 'ZZINREC'                                )
         RETURN

      END IF

C
C     Check the exclusion index.
C
      IF ( ( EXCLUD .LT. NONE ) .OR. ( EXCLUD .GT. 3 ) ) THEN
         
         CALL CHKIN  ( 'ZZINREC'                             )
         CALL SETMSG ( 'EXCLUD was #; allowed range is 0:3.' )
         CALL ERRINT ( '#', EXCLUD                           )
         CALL SIGERR ( 'SPICE(INDEXOUTOFRANGE)'              )
         CALL CHKOUT ( 'ZZINREC'                             )
         RETURN

      END IF

C
C     Get local copies of the coordinate bounds.
C 
      DO I = 1, 3   

         MINCOR(I) = BOUNDS( LOWER, I)
         MAXCOR(I) = BOUNDS( UPPER, I)
         L     (I) = MAXCOR(I) - MINCOR(I)

         IF ( L(I) .LT. 0.D0 ) THEN

            CALL CHKIN  ( 'ZZINREC'                              )
            CALL SETMSG ( 'Bounds are out of order for index #; '
     .      //            'bounds are #:#.'                      )
            CALL ERRDP  ( '#', BOUNDS(LOWER,I)                   )
            CALL ERRDP  ( '#', BOUNDS(UPPER,I)                   )
            CALL SIGERR ( 'SPICE(BOUNDSOUTOFORDER)'              )
            CALL CHKOUT ( 'ZZINREC'                              )
            RETURN

         END IF

      END DO

C
C     Compare coordinates to adjusted coordinate
C     boundaries.
C
      DO I = 1, 3

         IF ( EXCLUD .NE. I ) THEN
C
C           Create adjusted bounds for the Ith coordinate.
C 
            DELTA(I)  = MARGIN * ABS( L(I) )

            AMNCOR(I) = MINCOR(I) - DELTA(I)
            AMXCOR(I) = MAXCOR(I) + DELTA(I)

            IF (      ( P(I) .LT. AMNCOR(I) ) 
     .           .OR. ( P(I) .GT. AMXCOR(I) ) ) THEN

                RETURN

            END IF

         END IF

      END DO

C
C     All tests that were commanded have been passed. The input
C     point is considered to be contained in the expanded volume
C     element.
C
      INSIDE = .TRUE.

      RETURN
      END 
