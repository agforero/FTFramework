C$Procedure RC2GRD ( DSKBRIEF, rectangles to pixel grid )

      SUBROUTINE RC2GRD ( NREC,   BNDS1,  BNDS2, MAXGRD, MAXORD, 
     .                    VALUE,  ORD1,   ORD2,  CIVOR1, CIVOR2,  
     .                    PXMAP1, PXMAP2, NROWS, NCOLS,  GRID   )

C$ Abstract
C
C     Map rectangle list to pixel grid. Each pixel is entirely covered
C     by at least one input rectangle or not covered by any rectangle,
C     except on its boundary. The pixels covered by at least one
C     rectangle are marked with a specified logical value.
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
C     DSKBRIEF
C
C$ Declarations
 
      IMPLICIT NONE

      INTEGER               NREC
      DOUBLE PRECISION      BNDS1  ( 2, NREC )
      DOUBLE PRECISION      BNDS2  ( 2, NREC )
      INTEGER               MAXGRD
      INTEGER               MAXORD
      LOGICAL               VALUE
      INTEGER               ORD1   ( * )
      INTEGER               ORD2   ( * )
      INTEGER               CIVOR1 ( * )
      INTEGER               CIVOR2 ( * )
      INTEGER               PXMAP1 ( * )
      INTEGER               PXMAP2 ( * )
      INTEGER               NROWS
      INTEGER               NCOLS
      LOGICAL               GRID   ( * )
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NREC       I   Number of rectangles.
C     BNDS1      I   Bounds of the first coordinates of rectangles.
C     BNDS2      I   Bounds of the second coordinates of rectangles.
C     MAXGRD     I   Maximum size of output grid.
C     MAXORD     I   Maximum size of output order vectors.
C     VALUE      I   Value identifying inclusion in rectangles.
C     ORD1       O   Order vector of first coordinates.
C     ORD2       O   Order vector of second coordinates.
C     CIVOR1     O   Compressed inverse order vector of first 
C                    coordinates.
C     CIVOR2     O   Compressed inverse order vector of second 
C                    coordinates.
C     PXMAP1     O   Mapping from pixel first coordinates to
C                    rectangle first coordinates.
C     PXMAP2     O   Mapping from pixel second coordinates to
C                    rectangle second coordinates.
C     NROWS      O   Number of rows in the output grid.
C     NCOLS      O   Number of columns in the output grid.
C     GRID       O   Output grid.
C     
C$ Detailed_Input
C
C     NREC       is the number of rectangles in the input set.
C                These rectangles represent a spatial region to
C                be represented as a union of rectangular components 
C                that overlap only at their boundaries.
C
C                In DSKBRIEF, the rectangles represent a region having
C                no spatial coverage, also called a "gap region."
C
C     BNDS1,
C     BNDS2      are, respectively, the bounds of the first and second
C                coordinates of the rectangles. The elements
C     
C                   BNDS1(J,I), J = 1, 2
C                   BNDS2(J,I), J = 1, 2
C
C                are the lower and upper bounds of coordinates 1 and for
C                the Ithe rectangle.
C
C     MAXGRD     is the maximum size of the output grid.
C
C     MAXORD     is the maximum size of the output order vectors and
C                compressed inverse order vectors.
C
C     VALUE      is the logical value used in the pixel grid to mark
C                pixels that are in the region covered by the input
C                rectangles.
C
C$ Detailed_Output
C
C     ORD1,
C     ORD2       are, respectively, order vectors for 1-dimensional
C                arrays of all values of the first and second
C                coordinates. Upper and lower bounds are merged in
C                these arrays.
C         
C     CIVOR1,
C     CIVOR2     are, respectively, compressed inverse order vectors
C                for 1-dimensional arrays of all values of the first
C                and second coordinates. 
C
C     PXMAP1,
C     PXMAP2     are, respectively, arrays mapping pixel boundaries
C                to rectangle bounds. PXMAP1 maps pixel X-boundaries
C                to values in BNDS1; PXMAP2 maps pixel Y-boundaries
C                to values in BNDS2. 
C
C                Pixel X-coordinates range from 1 to NCOLS+1; pixel
C                Y-coordinate range from 1 to NROWS+1. 
C     
C     NROWS      is the number of rows in the output pixel grid.
C
C     NCOLS      is number of columns in the output pixel grid.
C
C     GRID       is a pixel grid. GRID is an array of logical values;
C                the dimensions of GRID are NROWS x NCOLS. Each pixel
C                corresponds to a rectangle in the input coordinate
C                space that is either entirely in one of the input
C                rectangles, or is, with the exception of its edges,
C                outside of all of them.
C
C                The elements of GRID having the value VALUE are those
C                pixels corresponding to rectangles in the input
C                coordinate space that are entirely contained in an
C                input rectangle.
C
C                DSKBRIEF uses the marked pixels of GRID to denote
C                pixels belonging to regions covered by a DSK
C                segment.
C   
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If NREC is less than 1, the error SPICE(VALUEOUTOFRANGE) is
C         signaled.
C
C     2)  If MAXGRD is less than the required size of the output grid,
C         the error SPICE(VALUEOUTOFRANGE) is signaled.
C
C     3)  If MAXORD is less than 1, the error SPICE(VALUEOUTOFRANGE) is
C         signaled.
C
C     4)  No input rectangle may have an edge of zero length. If any
C         pair of lower and upper coordinate bounds is not strictly
C         increasing, the error SPICE(INVALIDBOUNDS) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports determination of DSK spatial coverage gaps.
C
C     The coverage "gap region" of a set of DSK segments is that
C     region, within the rectangle defined by the global extrema of the
C     extents of the segments' first and second coordinates, where
C     there is no spatial coverage.
C
C     The pixel grid created by this routine has the property that each
C     pixel corresponds to a rectangle in the input coordinate space
C     that is either entirely in one of the input rectangles, or is,
C     with the exception of its edges, outside of all of them.
C 
C     The pixel grid simplifies the task of efficiently locating
C     "components" of the gap region, where components are rectangular
C     regions in the input coordinate space that overlap at most at
C     their edges.
C
C     Since the edges of pixels correspond to boundaries of the input
C     rectangles, it's easy to map a rectangular union of pixels
C     to a rectangle in the coordinate space of the input rectangles.
C     Each image rectangle is a component of the gap region.
C
C$ Examples
C
C     See usage in DSKBRIEF.
C
C$ Restrictions
C
C     1) For use only within program DSKBRIEF.
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
C     DSKBRIEF Version 2.0.0, 06-OCT-2016 (NJB)
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local variables
C      
      INTEGER               BSIZE
      INTEGER               COL
      INTEGER               I
      INTEGER               J
      INTEGER               K
      INTEGER               MAXPXX 
      INTEGER               MAXPXY 
      INTEGER               MINPXX 
      INTEGER               MINPXY 
      INTEGER               NGRID
      INTEGER               RNGMAX
      INTEGER               ROW

C
C     Saved values
C
      SAVE


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'RC2GRD' )

C
C     Check input size arguments for obvious initialization errors.
C
      IF ( NREC .LT. 1 ) THEN

         CALL SETMSG ( 'NREC is #; must be positive.' )
         CALL ERRINT ( '#', NREC                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'       )
         CALL CHKOUT ( 'RC2GRD'                       )
         RETURN

      END IF

      IF ( MAXGRD .LT. 1 ) THEN

         CALL SETMSG ( 'MAXGRD is #; must be positive.' )
         CALL ERRINT ( '#', MAXGRD                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'         )
         CALL CHKOUT ( 'RC2GRD'                         )
         RETURN

      END IF

      IF ( MAXORD .LT. 1 ) THEN

         CALL SETMSG ( 'MAXORD is #; must be positive.' )
         CALL ERRINT ( '#', MAXORD                      )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'         )
         CALL CHKOUT ( 'RC2GRD'                         )
         RETURN

      END IF

C
C     All input rectangle heights and widths must be strictly
C     positive.
C
      DO I = 1, NREC

         IF ( BNDS1(2,I) .LE. BNDS1(1,I) ) THEN

            CALL SETMSG ( 'BNDS1(2,#) = #; BNDS1(1,#) = #. Rectangle '
     .      //            'widths (and heights) must be positive.'    )
            CALL ERRINT ( '#', I                                      )
            CALL ERRDP  ( '#', BNDS1(2,I)                             )
            CALL ERRINT ( '#', I                                      )
            CALL ERRDP  ( '#', BNDS1(1,I)                             )
            CALL SIGERR ( 'SPICE(INVALIDBOUNDS)'                      )
            CALL CHKOUT ( 'RC2GRD'                                    )
            RETURN

         END IF

         IF ( BNDS2(2,I) .LE. BNDS2(1,I) ) THEN

            CALL SETMSG ( 'BNDS2(2,#) = #; BNDS2(1,#) = #. Rectangle '
     .      //            'heights (and widths) must be positive.'    )
            CALL ERRINT ( '#', I                                      )
            CALL ERRDP  ( '#', BNDS2(2,I)                             )
            CALL ERRINT ( '#', I                                      )
            CALL ERRDP  ( '#', BNDS2(1,I)                             )
            CALL SIGERR ( 'SPICE(INVALIDBOUNDS)'                      )
            CALL CHKOUT ( 'RC2GRD'                                    )
            RETURN

         END IF

 
      END DO

C
C     Find the order of the array of X bounds. We treat the array as a
C     one-dimensional array of length 2*NREC.
C
C     Produce the corresponding "compressed" inverse order vector. By
C     "compressed" we mean: suppose the set of input values were sorted
C     and compressed so that it contained no duplicates. For each
C     member of the original value array, map the member's index to the
C     index of the member in the compressed, sorted array. The
C     compressed inverse order vector contains this mapping.
C
      BSIZE = 2 * NREC

      CALL IOVCMP ( BNDS1, BSIZE, ORD1, CIVOR1, RNGMAX )
      
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'RC2GRD' )
         RETURN
      END IF

      DO I = 1, BSIZE

         J = ORD1(I)

         K = ((J-1)/2) + 1

      END DO
             
C
C     The width of the pixel grid is one less than the number of
C     distinct X bound values.
C     
      NCOLS = RNGMAX - 1

C
C     Get the order vector and compressed inverse order vector of
C     the Y bounds. (Note we have the same number of X and Y 
C     bounds.)
C
      CALL IOVCMP ( BNDS2, BSIZE, ORD2, CIVOR2, RNGMAX )
      
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'RC2GRD' )
         RETURN
      END IF

C
C     The height of the pixel grid is one less than the number of
C     distinct Y bound values.
C     
      NROWS = RNGMAX - 1

C
C     Check the grid size again, now that we know how large it
C     needs to be.
C
      NGRID = NROWS * NCOLS

      IF ( MAXGRD .LT. NGRID ) THEN

         CALL SETMSG ( 'MAXGRD is #; must be have size '
     .   //            'at least # in order to hold '
     .   //            'pixels for current set of '
     .   //            'rectangles.'                    )
         CALL ERRINT ( '#', MAXGRD                      )
         CALL ERRINT ( '#', NGRID                       )
         CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'         )
         CALL CHKOUT ( 'RC2GRD'                         )
         RETURN

      END IF
      
C
C     A program using this routine normally will need to map pixel
C     coordinates back to their corresponding d.p. values. Create
C     arrays to represent these mappings.
C     
C     Note that, since the bounds arrays are generally larger than the
C     corresponding pixel grid dimensions, the mappings we're about to
C     perform (which map the integers 1:BIZE into the ranges of the
C     compressed inverse order vectors) are not 1-1. They're still
C     valid; the process is just a bit ungainly because it can involve
C     overwriting elements of the output array. Each time this happens,
C     the affected output array element gets overwritten with the same
C     value it already had.
C
C     We'll store the mappings in the arrays PXMAP1 and PXMAP2.
C     The pixel coordinates
C
C        ( I, J ) 
C
C     correspond to the double precision coordinates
C
C        BNDS1( PXMAP1(I) )
C        BNDS2( PXMAP2(J) )
C
C     where we're treating BNDS1 and BNDS2 as one-dimensional
C     arrays of length BSIZE.
C

      DO I = 1, BSIZE

         PXMAP1( CIVOR1(I) ) = I
         PXMAP2( CIVOR2(I) ) = I

      END DO

C
C     Now map all rectangles to the integer indices of their
C     bounds in pixel space. Note that the pixel grid has
C     dimensions
C
C        ( NROWS, NCOLS )
C
C     and the ranges of the integer coordinates of the 
C     rectangle boundaries are
C
C        1 : NROWS + 1
C        1 : NCOLS + 1
C
C     We'll fill in the pixel grid to indicate which pixels are
C     covered by rectangles, and which ones lie in gaps.
C
C     Initialize the grid to indicate that it consists of one 
C     large gap.
C
      DO I = 1, NGRID

         GRID(I) = .NOT. VALUE

      END DO

C
C     For each input rectangle, mark the corresponding pixels
C     covered by the rectangle. Note that maximum pixel indices 
C     are less by one than those of the corresponding rectangle
C     upper bound indices.
C
      DO I = 1, NREC
C
C        Compute the bounds of the current rectangle in pixel
C        space. Recall that the all bounds for a given coordinate
C        (X or Y) are combined in a sequence of size 2*NREC.
C
         J      = 2*(I-1) + 1

         MINPXX = CIVOR1( J )
         MINPXY = CIVOR2( J )

         K      = 2*I

         MAXPXX = CIVOR1( K )
         MAXPXY = CIVOR2( K )
                  

         DO COL = MINPXX, MAXPXX-1

            DO ROW = MINPXY, MAXPXY-1
C
C              Mark the pixel at indices (ROW, COL) as
C              covered. 
C
               J       = NROWS*(COL-1) + ROW

               GRID(J) = VALUE

            END DO

         END DO

      END DO

      CALL CHKOUT ( 'RC2GRD' )
      RETURN
      END




