C$Procedure ZZRYTREC ( DSK, ray touches rectangular element )
 
      SUBROUTINE ZZRYTREC ( VERTEX, RAYDIR, BOUNDS, 
     .                      MARGIN, NXPTS,  XPT    )
 
C$ Abstract
C
C     SPICE Private routine intended solely for the support of SPICE
C     routines. Users should not call this routine directly due to the
C     volatile nature of this routine.
C
C     Find nearest intersection to ray's vertex of ray and
C     rectangular volume element. If the vertex is inside
C     the element, the vertex is considered to be the solution.
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
C     GEOMETRY
C     INTERCEPT
C     INTERSECTION
C     RAY
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      BOUNDS ( 2, 3 )
      DOUBLE PRECISION      MARGIN
      INTEGER               NXPTS
      DOUBLE PRECISION      XPT    ( 3 )
 

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     VERTEX     I   Ray's vertex.
C     RAYDIR     I   Ray's direction vector.
C     BOUNDS     I   Bounds of rectangular volume element.
C     MARGIN     I   Margin used for element expansion.
C     NXPTS      O   Number of intercept points.
C     XPT        O   Intercept.
C
C$ Detailed_Input
C
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                the ray to be used in the intercept computation. 
C
C                Both the vertex and ray direction must be represented
C                in the reference frame of the segment to which the
C                volume element boundaries correspond. The vertex is
C                considered to be an offset from the center of the
C                reference frame associated with the segment.
C 
C     BOUNDS     is a 2x3 array containing the bounds of a rectangular
C                volume element. Normally this is the coverage boundary
C                of a DSK segment. In the element
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
C$ Detailed_Output
C
C     XPT        is the intercept of the ray on the boundary of the
C                input volume element, if such an intercept exists. If
C                the ray's vertex is inside the element, XPT is set
C                equal to the vertex. XPT is valid if and only if FOUND
C                is .TRUE.

C                XPT is expressed in the reference frame associated
C                with the inputs VERTEX and RAYDIR. XPT represents
C                an offset from the origin of the coordinate system.
C
C                XPT is valid only if NXPTS is set to 1.
C
C
C     NXPTS      is the number of intercept points of the ray and
C                the volume element. 
C
C                Currently there are only two possible values for
C                NXPTS:
C
C                   1 for an intersection
C                   0 for no intersection
C    
C                If the vertex is inside the element, NXPTS is
C                set to 1.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1) If the minimum value of any coordinate exceeds the maximum,
C        the error will be signaled by a routine in the call tree of
C        this routine.
C          
C$ Files
C
C     None. However, the input segment boundaries normally have
C     been obtained from a loaded DSK file.
C     
C$ Particulars
C
C     This routine sits on top of data DSK type-specific ray-segment
C     intercept routines such as DSKX02.
C
C$ Examples
C
C     See usage in ZZDSKBUX.
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
C-    SPICELIB Version 1.0.0, 01-MAR-2016 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     find intercept of ray on rectangular volume element
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               NONE
      PARAMETER           ( NONE   = 0 )

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
      DOUBLE PRECISION      BOXORI ( 3 )
      DOUBLE PRECISION      EXTENT ( 3 )
      DOUBLE PRECISION      DELTA  ( 3 )
      DOUBLE PRECISION      L      ( 3 )

      INTEGER               I

      LOGICAL               FOUND
      LOGICAL               INSIDE


      IF ( RETURN() ) THEN
         RETURN
      END IF

C
C     Compute the original volume edge lengths from the coordinate
C     bounds.
C 
      DO I = 1, 3   

         L(I) = BOUNDS(UPPER, I) - BOUNDS(LOWER, I)

         IF ( L(I) .LE. 0.D0 ) THEN

            CALL CHKIN  ( 'ZZRYTREC'                               )
            CALL SETMSG ( 'Coordinate # bounds were #:#; bounds '
     .      //            'must be strictly increasing.'           )
            CALL ERRINT ( '#',  I                                  )
            CALL ERRDP  ( '#',  BOUNDS(LOWER,I)                    )
            CALL ERRDP  ( '#',  BOUNDS(UPPER,I)                    )
            CALL SIGERR ( 'SPICE(BADCOORDBOUNDS)'                  )
            CALL CHKOUT ( 'ZZRYTREC'                               )
            RETURN

         END IF

      END DO

C
C     Determine whether the vertex is inside the element.
C     Use double the margin for this test, since we don't
C     want to have false negative tests for rays having
C     vertices lying on the expanded element boundary.
C
      NXPTS = 0

      CALL ZZINREC ( VERTEX, BOUNDS, 2*MARGIN, NONE, INSIDE )

      IF ( INSIDE ) THEN
C
C        We know the answer.
C
         NXPTS = 1

         CALL VEQU ( VERTEX, XPT )

         RETURN

      END IF
 
C
C     Expand the box using the specified margin.
C
      DO I = 1, 3

         DELTA(I)  = MARGIN * ABS( L(I) ) 

         BOXORI(I) = BOUNDS(LOWER,I) - DELTA(I)

         EXTENT(I) = L(I) + ( 2*DELTA(I) )

      END DO

C
C     Find the ray-surface intercept on the expanded element,
C     if the intercept exists.
C
      CALL ZZRAYBOX ( VERTEX, RAYDIR, BOXORI, EXTENT, XPT, FOUND )

      IF ( FOUND ) THEN
         NXPTS = 1
      END IF

      RETURN
      END 
