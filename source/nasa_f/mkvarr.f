C$Procedure MKVARR ( MKDSK: make vertex array from grid data )

      SUBROUTINE MKVARR ( INFILE, AUNITS, DUNITS, ROWMAJ, TOPDWN, 
     .                    LEFTRT, CORSYS, CORPAR, REFVAL, HSCALE,
     .                    NCOLS,  NROWS,  LFTCOR, TOPCOR, COLSTP,
     .                    ROWSTP, MAXNV,  VERTS                  )
      
C$ Abstract
C
C     Make normalized vertex array from a regular height grid. Elements
C     of the output array are in row-major, top-down, left-right order.
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
C     MKDSK
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'

      CHARACTER*(*)         INFILE
      CHARACTER*(*)         AUNITS
      CHARACTER*(*)         DUNITS
      LOGICAL               ROWMAJ
      LOGICAL               TOPDWN
      LOGICAL               LEFTRT
      INTEGER               CORSYS
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      REFVAL
      DOUBLE PRECISION      HSCALE
      INTEGER               NCOLS
      INTEGER               NROWS
      DOUBLE PRECISION      LFTCOR
      DOUBLE PRECISION      TOPCOR
      DOUBLE PRECISION      COLSTP
      DOUBLE PRECISION      ROWSTP
      INTEGER               MAXNV
      DOUBLE PRECISION      VERTS  ( 3, * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     INFILE     I   Name of input file.
C     AUNITS     I   Angular units.
C     DUNITS     I   Distance units.
C     ROWMAJ     I   Flag indicting whether grid is row-major.
C     TOPDWN     I   Flag indicting whether grid is top-down.
C     LEFTRT     I   Flag indicting whether grid is left-right.
C     CORSYS     I   Coordinate system.
C     CORPAR     I   Coordinate parameters.
C     REFVAL     I   Height reference value.
C     HSCALE     I   Height scale.
C     NROWS      I   Number of rows in grid.
C     NCOLS      I   Number of columns in grid.
C     LFTCOR     I   Coordinate of leftmost column of grid.
C     TOPCOR     I   Coordinate of top row of grid.
C     COLSTP     I   Column step size.
C     ROWSTP     I   Row step size.
C     MAXNV      I   Maximum number of vertices to return.
C     VERTS      O   Vertex array.
C
C$ Detailed_Input
C 
C     INFILE     is the name of an input data file containing height
C                grid data.
C
C
C     AUNITS     is the name of the angular unit associated with the
C                grid coordinates, if the grid coordinate system is
C                latitudinal or planetodetic. AUNITS must be supported
C                by the SPICELIB routine CONVRT.
C
C     DUNITS     is the name of the distance unit associated with the
C                grid coordinates. DUNITS must be supported by the
C                SPICELIB routine CONVRT.
C
C     ROWMAJ     is a logical flag that is set by the caller to .TRUE.
C                when the input grid data are organized in row-major
C                order, and that is set to .FALSE. when the data are in
C                column-major order.
C
C                "Row-major" means that the Nth consecutive sequence of
C                NCOLS tokens in the file represents the Nth row of
C                data in the output grid.     
C
C                "Column-major" means that the Nth consecutive sequence
C                of NROWS tokens in the file represents the Nth column
C                of data in the output grid.
C
C                Note that the mapping from a token's position in the
C                input file to its position in the output grid is not
C                defined by ROWMAJ alone: the values of TOPDWN and
C                LEFTRT are needed as well to fully define the mapping.
C
C
C     TOPDWN     is a logical flag that is set by the caller to .TRUE.
C                when the input file contains data in top-down order,
C                and is set to .FALSE. when the input data are in
C                bottom-up order. Here "top" means "having the highest
C                value of the second coordinate."
C
C                TOPDWN is true if and only if the datum of the top
C                row and Nth column in the output grid precedes the
C                datum for any other row and Nth column in the input
C                file. In other words, the data from the input file
C                fill in the rows of the output grid in "top-down"
C                order.
C
C                When the input file is row-major and TOPDWN is .TRUE.,
C                the top row of the output grid contains the data from
C                the first NCOLS tokens in the input file. When the
C                input file is column-major and TOPDWN is .TRUE., the
C                columns themselves are in top-down order: the first
C                element of each sequence of NROWS tokens belongs to
C                the top row of the output grid.
C
C     LEFTRT     is a logical flag that is set by the caller to .TRUE.
C                when the input file contains data in left-right order,
C                and is set to .FALSE. when the input data are in
C                right-left order. Here "left" means "having the lowest
C                value of the first coordinate."
C
C                LEFTRT is true if and only if the datum of the left
C                column and Nth row in the output grid precedes the
C                datum for any other column and Nth row in the input
C                file. In other words, the data from the input file
C                fill in the columns of the output grid in "left-right"
C                order.
C
C                When the input file is column-major and LEFTRT is
C                .TRUE., the leftmost column of the output grid
C                contains the data from the first NROWS tokens in the
C                input file. When the input file is row-major and LEFT
C                is .TRUE., the rows themselves are in left-right
C                order: the first element of each sequence of NCOLS
C                tokens belongs to the left column of the output grid.
C
C     CORSYS     is a DSK subsystem code designating the coordinate
C                system of the input coordinates.
C
C     CORPAR     is an array containing parameters associated with
C                the input coordinate system. The contents of the
C                array are as described in the DSK include file
C                dskdsc.inc.
C
C     REFVAL     is a reference value to be added to the input height.
C                REFVAL is used only for latitudinal and rectangular
C                coordinates.
C
C                REFVAL must be non-negative.
C
C                Units are km.
C
C     HSCALE     is a conversion factor to be applied to height data. 
C                Multiplying a datum by HSCALE converts the datum to
C                kilometers. HSCALE need not correspond to a standard
C                unit.
C
C     NROWS      is the number of rows in the output grid.
C
C     NCOLS      is the number of columns in the output grid.
C
C     LFTCOR     is the coordinate of the leftmost column of a
C                rectangular data grid. If the coordinate system is
C                latitudinal, this is the minimum longitude on the
C                grid. If the system is rectangular, this is the
C                minimum X value on the grid.
C     
C                Units are given by AUNITS or DUNITS, as applicable.
C
C     TOPCOR     is the coordinate of the top row of a rectangular data
C                grid. If the coordinate system is latitudinal, this is
C                the maximum latitude on the grid. If the system is
C                rectangular, this is the maximum Y value on the grid.
C   
C                Units are given by AUNITS or DUNITS, as applicable.
C
C     COLSTP     is the uniform separation between columns of the grid.
C
C                Units are given by AUNITS or DUNITS, as applicable.
C
C     ROWSTP     is the uniform separation between rows of the grid.
C
C                Units are given by AUNITS or DUNITS, as applicable.
C    
C     MAXNV      is the maximum number of vertices that can be placed
C                in the output array.
C
C
C$ Detailed_Output
C
C     VERTS      is an array of 3-vectors corresponding to the height
C                data in the input file. The data in VERTS always 
C                represent a grid organized in
C
C                   row-major
C                   top-down
C                   left-right
C
C                order, regardless of the organization of the input
C                file.
C                  
C                VERTS should be declared by the caller as:
C     
C                   DOUBLE PRECISION VERTS ( 3, MAXNV )
C
C                Units are always km.
C
C$ Parameters
C
C     The input file must have a maximum line length of LNSIZE.
C
C     See mkdsk.inc.
C
C$ Exceptions
C
C     1)  If the number of data values read does not match the value
C
C             NROWS * NCOLS
C          
C         the error SPICE(INVALIDDATACOUNT) is signaled.
C
C     2)  If an error occurs while converting angular units to radians
C         or distance units to km, the error will be diagnosed by a
C         routine in the call tree of this routine.
C
C     3)  If MAXNV indicates the output buffer is too small to hold the
C         output array, the error SPICE(BUFFERTOOSMALL) is signaled.
C
C     4)  If REFVAL is negative, the error is diagnosed by routines
C         in the call tree of this routine.
C     
C     5)  If the input coordinate parameters in CORPAR are invalid, the
C         error is diagnosed by routines in the call tree of this
C         routine.
C
C     6)  If an error occurs while reading the input file, the
C         error is diagnosed by routines in the call tree of this
C         routine.
C
C$ Files
C
C     The file specified by INFILE can have any of the attributes (one
C     choice from each row below):
C
C        row-major  or column-major
C        top-down   or bottom-up
C        left-right or right-left
C
C     The number of tokens per line may vary. The number need have no
C     particular relationship to the row or column dimensions of the
C     output grid.
C
C     The file must contain only tokens that can be read as double
C     precision values. No non-printing characters can be present in
C     the file.
C
C     Tokens can be delimited by blanks or commas. Tokens must not be
C     split across lines.
C
C     Blank lines are allowed; however, their use is discouraged
C     because they'll cause line numbers in diagnostic messages to
C     be out of sync with actual line numbers in the file.
C
C     The file must end with a line terminator.
C
C$ Particulars
C
C     The "output grid" is the two-dimensional array of 3-dimensional
C     vertices; the row count is NROWS and the column count is NCOLS.
C
C     The output grid is "normalized" in the sense that it always is 
C
C        row-major
C        top-down
C        left-right
C
C$ Examples
C
C     See usage in the MKDSK routine MKGRID.
C
C$ Restrictions
C
C     1) For use only within program MKDSK.
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
C-    MKDSK Version 1.0.0, 25-FEB-2017 (NJB)
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               BUFSIZ
      PARAMETER           ( BUFSIZ = 1000 )

C
C     Local variables
C
      DOUBLE PRECISION      ASCALE
      DOUBLE PRECISION      COORDS ( 2 )
      DOUBLE PRECISION      CORSCL
      DOUBLE PRECISION      DSCALE
      DOUBLE PRECISION      HEIGHT
      DOUBLE PRECISION      HREF
      DOUBLE PRECISION      LEFTCO
      DOUBLE PRECISION      HSTEP
      DOUBLE PRECISION      TOPCO
      DOUBLE PRECISION      VALUES ( BUFSIZ )
      DOUBLE PRECISION      VSTEP

      INTEGER               COL
      INTEGER               I
      INTEGER               J
      INTEGER               K
      INTEGER               N
      INTEGER               NC
      INTEGER               NR
      INTEGER               NV
      INTEGER               ROW
      INTEGER               TOTAL

      LOGICAL               DONE

C
C     Saved variables
C
      SAVE                  VALUES

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'MKVARR' )

C
C     Compute unit scales.
C
      CALL CONVRT ( 1.D0, AUNITS, 'RADIANS', ASCALE )
      CALL CONVRT ( 1.D0, DUNITS, 'KM',      DSCALE )
      
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'MKVARR' )
         RETURN
      END IF

C
C     Convert grid boundary and step values to radians and km,
C     as needed. Convert the reference value as well.
C
      IF ( CORSYS .EQ. RECSYS ) THEN
         
         CORSCL = DSCALE
      ELSE
         CORSCL = ASCALE
      END IF

      LEFTCO = CORSCL * LFTCOR
      TOPCO  = CORSCL * TOPCOR
      HSTEP  = CORSCL * COLSTP
      VSTEP  = CORSCL * ROWSTP

      HREF   = DSCALE * REFVAL
      
C
C     Presume the row and column counts are correct.
C
      NV = NROWS * NCOLS
C
C     Make sure we have room for the output array.
C
      IF ( NV .GT. MAXNV ) THEN

         CALL SETMSG ( 'Room for # vertices is needed; '
     .   //            'amount available is #.'         )
         CALL ERRINT ( '#', NV                          )
         CALL ERRINT ( '#', MAXNV                       )
         CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'           )
         CALL CHKOUT ( 'MKVARR'                         )
         RETURN

      END IF

C
C     Fetch data from the input file; distribute it to the
C     vertex array. Grab the first buffer of data:
C
      CALL RDGRD5 ( INFILE, BUFSIZ, N, VALUES, DONE )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'MKVARR' )
         RETURN
      END IF


      TOTAL = 0

      DO WHILE ( ( TOTAL .LT. NV ) .AND. ( N .GT. 0 )  )
C
C        Distribute the data we've buffered.
C
         DO I = 1, N

            TOTAL = TOTAL + 1
C
C           Compute the row and column of the current
C           vertex, then covert those indices to a one-
C           dimensional index.
C
            IF ( ROWMAJ ) THEN
C
C              We fill in the vertex grid one row at a time.
C
               IF ( TOPDWN ) THEN
C
C                 We fill in rows from the top down.
C
                  ROW = ( ( TOTAL - 1 )/NCOLS  ) + 1

                  IF ( LEFTRT ) THEN
C
C                    We fill in rows from left to right.
C
                     COL = TOTAL - ( (ROW-1)*NCOLS )

                  ELSE
C
C                    We fill in rows from right to left.
C
                     J   = TOTAL - ( (ROW-1)*NCOLS )

                     COL = NCOLS + 1 - J
                     
                  END IF


               ELSE
C
C                 We fill in rows from the bottom up.
C
                  NR  = ( ( TOTAL - 1 )/NCOLS  ) + 1

                  ROW = NROWS + 1 - NR


                  IF ( LEFTRT ) THEN
C
C                    We fill in rows from left to right.
C
                     COL = TOTAL - ( (NR-1)*NCOLS )

                  ELSE
C
C                    We fill in rows from right to left.
C
                     J   = TOTAL - ( (NR-1)*NCOLS )

                     COL = NCOLS + 1 - J
                     
                  END IF


               END IF

        
            ELSE
C
C              We fill in the vertex grid one column at a time.
C
               IF ( LEFTRT ) THEN
C
C                 We fill in columns from left to right.
C
                  COL = ( ( TOTAL - 1 )/ NROWS ) + 1


                  IF ( TOPDWN ) THEN
C
C                    We fill in columns from the top down.
C
                     ROW = TOTAL - ( (COL-1)*NROWS )

                  ELSE
C
C                    We fill in columns from the bottom up.
C
                     J   = TOTAL - ( (COL-1)*NROWS )

                     ROW = NROWS + 1 - J
                     
                  END IF


               ELSE
C
C                 We fill in columns from right to left.
C
                  NC  = ( ( TOTAL - 1 )/ NROWS ) + 1

                  COL = NCOLS + 1 - NC

                  IF ( TOPDWN ) THEN
C
C                    We fill in columns from the top down.
C
                     ROW = TOTAL - ( (NC-1)*NROWS )

                  ELSE
C
C                    We fill in columns from the bottom up.
C
                     J   = TOTAL - ( (NC-1)*NROWS )

                     ROW = NROWS + 1 - J
                     
                  END IF

               END IF

            END IF

C
C           Compute the domain coordinates for this vertex.
C
            CALL RC2COR ( LEFTCO, TOPCO, HSTEP, 
     .                    VSTEP,  COL,   ROW,   COORDS )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'MKVARR' )
               RETURN
            END IF

C
C           Compute the 1-D vertex index.
C
            K = ( ( ROW - 1 ) * NCOLS )  +  COL

C
C           Compute the vertex and insert it into the vertex array.
C
            HEIGHT = VALUES(I) * HSCALE

            CALL MAKVTX ( CORSYS, CORPAR, COORDS,  
     .                    HREF,   HEIGHT, VERTS(1,K) )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'MKVARR' )
               RETURN
            END IF

         END DO

C
C        Try to read more data, but only if the EOF hasn't been
C        reached. Note that reading after EOF has been reached
C        will initiate a new read of the file, starting with the
C        first line.
C
         IF ( .NOT. DONE ) THEN
            CALL RDGRD5 ( INFILE, BUFSIZ, N, VALUES, DONE )
         ELSE
            N = 0
         END IF

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'MKVARR' )
            RETURN
         END IF

      END DO

C
C     We've read the data. Make sure we got the number of vertices
C     we were expecting.
C
      IF ( TOTAL .NE. NV ) THEN

         CALL SETMSG ( 'Column count = #; row count = #. Expected '
     .   //            'height values for # vertices but found '
     .   //            'data for #. Setup file and data file '
     .   //            ' are inconsistent.'                        )
         CALL ERRINT ( '#', NCOLS                                  )
         CALL ERRINT ( '#', NROWS                                  )
         CALL ERRINT ( '#', NV                                     )
         CALL ERRINT ( '#', TOTAL                                  )
         CALL SIGERR ( 'SPICE(INVALIDDATACOUNT)'                   )
         CALL CHKOUT ( 'MKVARR'                                    )
         RETURN

      END IF

      CALL CHKOUT ( 'MKVARR' )
      RETURN
      END
