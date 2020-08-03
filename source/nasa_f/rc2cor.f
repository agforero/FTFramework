C$Procedure RC2COR ( MKDSK: map grid row and column to coordinates )

      SUBROUTINE RC2COR ( LFTCOR, TOPCOR, COLSTP, 
     .                    ROWSTP, COL,    ROW,    COORDS )

C$ Abstract
C
C     Map grid row and column indices to standard coordinates.
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
C     MKDSK
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'

      DOUBLE PRECISION      LFTCOR
      DOUBLE PRECISION      TOPCOR
      DOUBLE PRECISION      COLSTP
      DOUBLE PRECISION      ROWSTP
      INTEGER               COL
      INTEGER               ROW
      DOUBLE PRECISION      COORDS ( 2 )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     LFTCOR     I   Name of input file.
C     TOPCOR     I   Maximum number of values to return.
C     COLSTP     I   Column step size.
C     ROWSTP     I   Row step size.
C     COL        I   Column index.
C     ROW        I   Row index.
C     COORDS     O   Coordinates of input grid point.
C
C$ Detailed_Input
C 
C     LFTCOR     is the coordinate of the leftmost column of a
C                rectangular data grid. If the coordinate system is
C                latitudinal, this is the minimum longitude on the
C                grid. If the system is rectangular, this is the
C                minimum X value on the grid.
C     
C                Units are caller-defined but must be consistent
C                with those of COLSTP.
C
C     TOPCOR     is the coordinate of the top row of a rectangular data
C                grid. If the coordinate system is latitudinal, this is
C                the maximum latitude on the grid. If the system is
C                rectangular, this is the maximum Y value on the grid.
C   
C                Units are caller-defined but must be consistent
C                with those of ROWSTP.
C
C     COLSTP     is the uniform separation between columns of the grid.
C
C                Units are caller-defined but must be consistent
C                with those of LFTCOR.

C     ROWSTP     is the uniform separation between rows of the grid.
C
C                Units are caller-defined but must be consistent
C                with those of TOPCOR.
C
C
C     COL        is the index of a column in the grid. Indices are 
C                1-based: the first coordinate of the first column
C                is LFTCOR.
C
C     ROW        is the index of a row in the grid. Indices are 
C                1-based: the second coordinate of the top row
C                is TOPCOR.
C
C     
C$ Detailed_Output
C
C     COORDS     is an array containing the coordinates corresponding
C                to ROW and COL. The first element of COORDS is the
C                coordinate corresponding to COL; the second element
C                corresponds to ROW.
C
C                For example, if the coordinate system is latitudinal,
C                COORDS(1) is the longitude corresponding to COL and
C                COORDS(2) is the latitude corresponding to ROW.
C     
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If either the row or column step is not strictly positive,
C         the error SPICE(INVALIDSTEP) is signaled.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     The computation performed by this routine refers to a rectangular
C     grid in a latitudinal, planetodetic, or rectangular coordinate
C     system. The coordinate system is not specified but must be
C     the same for all inputs.
C
C     All rows have equal separation in latitude or Y-coordinate,
C     depending on the coordinate system used. All columns have equal
C     separation in longitude or X-coordinate.
C
C     Column 1 of the grid corresponds to the coordinate value specified
C     by LFTCOR.
C
C     Row 1 of the grid corresponds to the coordinate value specified
C     by TOPCOR. Thus the second coordinate of the rows decreases with
C     increasing row index.
C
C$ Examples
C
C     See usage in the MKDSK routine MKVARR.
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
      LOGICAL               RETURN
      
C
C     Local variables
C
       
      IF ( RETURN() ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'RC2COR' )

C
C     Make sure the step sizes are valid.
C
      IF ( COLSTP .LE. 0.D0 ) THEN

         CALL SETMSG ( 'Column step was #; must be strictly positive.' )
         CALL ERRDP  ( '#', COLSTP                                     )
         CALL SIGERR ( 'SPICE(INVALIDSTEP)'                            )
         CALL CHKOUT ( 'RC2COR'                                        )
         RETURN

      END IF

      IF ( ROWSTP .LE. 0.D0 ) THEN

         CALL SETMSG ( 'Row step was #; must be strictly positive.' )
         CALL ERRDP  ( '#', ROWSTP                                  )
         CALL SIGERR ( 'SPICE(INVALIDSTEP)'                         )
         CALL CHKOUT ( 'RC2COR'                                     )
         RETURN

      END IF

C
C     In a lon/lat system, the leftmost column is at minimum longitude;
C     the top row is at maximum latitude.
C
C     In the rectangular system, the leftmost column is at minimum
C     X; the top row is at maximum Y.
C         
C     For all systems, the computation is the same.
C
      COORDS(1) = LFTCOR + ( ( COL - 1 )*COLSTP )
      COORDS(2) = TOPCOR - ( ( ROW - 1 )*ROWSTP )

      CALL CHKOUT ( 'RC2COR' )
      RETURN
      END 
