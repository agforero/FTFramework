C$Procedure FNDCMP ( DSKBRIEF, find rectangular components )

      SUBROUTINE FNDCMP ( NROWS, NCOLS,  VALUE,  MAXN, 
     .                    GRID,  VSET,   MRKSET, TMPSET,
     .                    NCOMP, MINPXX, MAXPXX, MINPXY, MAXPXY )

C$ Abstract
C
C     Find rectangular components in a 2-d logical grid.
C
C     The result is a list of rectangles. For each output rectangle,
C     the grid value is equal to the input VALUE at each pixel in that
C     rectangle. The rectangles are described by the bounds of their X
C     and Y pixel coordinates.
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

      INTEGER               LBCELL
      PARAMETER           ( LBCELL = -5 )

      INTEGER               NROWS
      INTEGER               NCOLS
      LOGICAL               GRID   ( NROWS, NCOLS )
      LOGICAL               VALUE
      INTEGER               MAXN
      INTEGER               VSET   ( LBCELL : * )
      INTEGER               MRKSET ( LBCELL : * )
      INTEGER               TMPSET ( LBCELL : * )
      INTEGER               NCOMP
      INTEGER               MINPXX   ( * )
      INTEGER               MAXPXX   ( * )
      INTEGER               MINPXY   ( * )
      INTEGER               MAXPXY   ( * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     NROWS      I   Number of rows in the input grid.
C     NCOLS      I   Number of columns in the input grid.
C     VALUE      I   Value marking pixels to aggregate.
C     MAXN       I   Maximum number of components.
C     GRID      I-O  Input grid.
C     VSET      I-O  Set of pixels set to VALUE.
C     MRKSET    I-O  Set of marked pixels.
C     TMPSET    I-O  Workspace set.
C     NCOMP      O   Number of rectangular components.
C     MINPXX     O   Minimum X coordinates of components.
C     MAXPXX     O   Maximum X coordinates of components.
C     MINPXY     O   Minimum Y coordinates of components.
C     MAXPXY     O   Maximum Y coordinates of components.
C     
C$ Detailed_Input
C
C     NROWS      is the number of rows in the input pixel grid.
C
C     NCOLS      is number of columns in the input pixel grid.
C
C     VALUE      is the logical value used in the pixel grid to mark
C                pixels that are to be grouped into rectangular
C                components.
C
C     MAXN       is the maximum number of rectangular components that
C                can be returned.
C
C     GRID       is, on input, a pixel grid. GRID is an array of
C                logical values; the dimensions of GRID are NROWS x
C                NCOLS.
C
C                Caution: GRID is modified by this routine.
C
C     VSET       is, on input, an initialized integer set. VSET will
C                be filled with the 1-dimensional indices of pixels
C                in GRID that are set to the value VALUE.
C
C     MRKSET     is, on input, an initialized integer set. MRKSET is a
C                workspace variable: it will be filled with the
C                1-dimensional indices of pixels in GRID that belong to
C                a rectangular component that under construction.
C
C     TMPSET     is, on input, an initialized integer set. TMPSET is a
C                temporary variable.
C
C$ Detailed_Output
C
C     GRID       is an overwritten version of the input pixel grid. 
C                Pixels belonging to components are marked by negating
C                those pixels' values.
C
C     VSET       is a modified integer set. The contents of VSET are
C                undefined.
C
C     MRKSET     is a modified integer set. The contents of MRKSET are
C                undefined.
C
C     TMPSET     is a modified integer set. The contents of TMPSET are
C                undefined.
C     
C     NCOMP      is the number of rectangular components making up
C                the set of marked pixels.
C
C     MINPXX     is an array of lower bounds of the X coordinates of
C                the rectangular components. The Ith element of MINPXX
C                is the lower X bound of the Ith component.
C
C     MAXPXX     is an array of upper bounds of the X coordinates of
C                the rectangular components. The Ith element of MAXPXX
C                is the upper X bound of the Ith component.
C
C     MINPXY     is an array of lower bounds of the Y coordinates of
C                the rectangular components. The Ith element of MINPXY
C                is the lower Y bound of the Ith component.
C
C     MAXPXY     is an array of upper bounds of the Y coordinates of
C                the rectangular components. The Ith element of MAXPXY
C                is the upper Y bound of the Ith component.
C        
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If MAXN is smaller than the number of rectangular
C         components found, the error SPICE(ARRAYTOOSMALL) is
C         signaled.
C 
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports determination of DSK spatial coverage gaps.
C     It can just as easily be used to determine regions of spatial
C     coverage.
C
C     The marked pixels in the input pixel grid are grouped by this
C     routine into disjoint, rectangular sets called "components." The
C     output coordinate bounds returned by this routine define the
C     boundaries of these components.
C
C     This routine is designed to operate efficiently. While
C     constructing a component, it accumulates pixels in increasing
C     order of the pixels' 1-dimensional indices. The allows the
C     routine to avoid sorting the indices in order to make the
C     component's pixel indices into a set. However, this algorithm
C     constrains the way components are identified: their pixels are
C     accumulated in column-major order.
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
C-    DSKBRIEF Version 1.0.0, 30-JAN-2017 (NJB)
C
C        Original version 03-SEP-2016 (NJB)
C-&


C
C     SPICELIB functions
C
      INTEGER               CARDI

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local variables
C
      INTEGER               COL
      INTEGER               COLSIZ
      INTEGER               I
      INTEGER               ID
      INTEGER               J
      INTEGER               MAXCOL
      INTEGER               MAXROW
      INTEGER               MINCOL
      INTEGER               MINROW
      INTEGER               REMAIN
      INTEGER               ROW

      LOGICAL               FOUND


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'FNDCMP' )

      CALL SCARDI ( 0, VSET   )
      CALL SCARDI ( 0, MRKSET )
      CALL SCARDI ( 0, TMPSET )

C
C     First step: make a pass through the grid, and store the ID
C     of each pixel matching VALUE.
C
C     Proceed in column-major order.
C
      DO COL = 1, NCOLS

         DO ROW = 1, NROWS

            IF ( GRID(ROW,COL) .EQV. VALUE ) THEN
C
C              It's a match.
C
C              ID is the one-dimensional index of the current pixel.
C
               ID = (COL-1)*NROWS + ROW

C
C              Since we're traversing the grid in increasing ID 
C              order, the elements of VSET will automatically be
C              in increasing order. We don't need to sort them.
C
               CALL APPNDI ( ID, VSET )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( 'FNDCMP' )
                  RETURN
               END IF

            END IF

         END DO

      END DO

C
C     Now find rectangular sets of pixels equal to VALUE.
C
      CALL SCARDI ( 0, MRKSET )

      REMAIN = CARDI( VSET )

      NCOMP  = 0

      DO WHILE ( REMAIN .GT. 0 )
C
C        Get the row and column coordinates of the first pixel in VSET.
C
         ID     = VSET(1) 

         COL    = ( ( ID - 1 ) / NROWS ) + 1

         ROW    = ID - ( (COL-1)*NROWS )

         MINROW = ROW
         MINCOL = COL
C
C        We'll extend the component in the direction of higher row
C        indices as far as possible, then in the direction of higher
C        column indices  as far as possible. The reason for this is
C        that we want to accumulate pixels in increasing order of ID.
C         
         MAXROW = NROWS
         MAXCOL = COL
         FOUND  = .TRUE.
         
         DO WHILE (  ( COL .LE. NCOLS ) .AND. FOUND )
C
C           COL is a valid column number at the top of the loop.
C           We increment COL at the bottom of the loop.
C
C           Initialize ROW for a pass through the current column.
C           
            ROW = MINROW - 1

C
C           Caution: the value of MAXROW in the loop termination
C           condition below changes during loop execution! The
C           value is NROWS on the first pass; then it changes
C           to the maximum row number of the first column of the
C           component.
C
            DO WHILE ( ( ROW .LT. MAXROW ) .AND. FOUND )
C
C              Note the .LT. operator in the loop termination
C              condition. We increment ROW at the top of the
C              loop, so the value of ROW is correct after
C              loop termination.
C              
               ROW   = ROW + 1

               FOUND = GRID(ROW,COL) .EQV. VALUE 

            END DO

            IF ( COL .EQ. MINCOL ) THEN
C
C              The index of the last row that matched becomes the
C              maximum row index of this component.
C
               IF ( FOUND ) THEN
C
C                 The row index reached NROWS.
C                 
                  MAXROW = NROWS

               ELSE
C
C                 The last matching row was the one preceding ROW.
C
                  MAXROW = ROW - 1
C
C                 Set FOUND to .TRUE. so we'll go on to look at the
C                 next column.
C                 
                  FOUND = .TRUE.

               END IF
C
C              Now we know the size of the columns of the component.
C              
               COLSIZ = MAXROW - MINROW + 1
C
C              Always go on to look at the next column. FOUND is 
C              .TRUE. at this point.
C
               MAXCOL = COL
               
            ELSE 
C
C              After we process the first column of the component,
C              we don't adjust MAXROW again. It's set to the highest
C              row number of the first column of the component.
C
               IF ( .NOT. FOUND ) THEN
C
C                 The current column fails to match in some row of the
C                 current column. This column can't be included in the
C                 component.
C
                  MAXCOL = COL - 1
                  
               ELSE
C
C                 The current column matches from row indices MINROW
C                 to MAXROW. This column is part of the component.
C
                  MAXCOL = COL
C
C                 Set FOUND to .TRUE. so we'll go on to look at the
C                 next column.
C                 
                  FOUND = .TRUE.                 

               END IF

            END IF

            
            IF ( FOUND ) THEN
C
C              We've found a column of matching pixels in the 
C              current component. 
C
C              Add the pixels of the column to the marked set.
C
C              Let ID be the ID of the first pixel of the column.
C
               ID = (MAXCOL-1)*NROWS + MINROW

               DO I = 1, COLSIZ

                  J = ID - 1 + I

                  CALL APPNDI ( J, MRKSET )

                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( 'FNDCMP' )
                     RETURN
                  END IF

C
C                 Fill in the matching pixels so they won't
C                 match again.
C
                  GRID( MINROW-1+I, MAXCOL ) = .NOT. VALUE

               END DO
C
C              Note that we've added the IDs to MRKSET in increasing
C              order, so MRKSET remains a set. We don't need to sort
C              its contents.
C
C              Prepare to examine the next column. 
C
               REMAIN = REMAIN - COLSIZ
               COL    = COL    + 1

            END IF

         END DO
C
C        We've finished building a component.
C
         NCOMP = NCOMP + 1

C
C        Update VSET: subtract the pixels of the new component.
C        
C        Note that subtracting one set from another should be an
C        efficient process, if done correctly. We trust DIFFI to manage
C        this.
C 
         CALL DIFFI ( VSET,   MRKSET, TMPSET )
         CALL COPYI ( TMPSET,         VSET   )

         CALL SCARDI ( 0, MRKSET )
         CALL SCARDI ( 0, TMPSET )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'FNDCMP' )
            RETURN
         END IF


C        The bounds of the component we just found are given by
C
C           MINROW, MAXROW, MINCOL, MAXCOL
C
         IF ( NCOMP .LE. MAXN ) THEN

            MINPXX(NCOMP) = MINCOL
            MAXPXX(NCOMP) = MAXCOL
            MINPXY(NCOMP) = MINROW
            MAXPXY(NCOMP) = MAXROW

         ELSE
C
C           We're out of room.
C           
            CALL SETMSG ( 'There are more output rectangles than '
     .      //            'can be accommodated in the output '
     .      //            'rectangle boundary arrays. So far, '
     .      //            '# components have been found; the '
     .      //            'maximum supported number is #.'        )
            CALL ERRINT ( '#', NCOMP                              )
            CALL ERRINT ( '#', MAXN                               )
            CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                  )
            CALL CHKOUT ( 'FNDCMP'                                )
            RETURN

         END IF

      END DO
      

      CALL CHKOUT ( 'FNDCMP' )
      RETURN
      END




