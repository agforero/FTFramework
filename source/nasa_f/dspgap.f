C$Procedure   DSPGAP ( DSKBRIEF, display spatial coverage gaps )
 
      SUBROUTINE DSPGAP ( GAP, CORSYS, NSIG, NCOMP, BDS1, BDS2 )
 
C$ Abstract
C
C     Display spatial coverage gaps of a collection of DSK segments.
C     The gaps are expressed as a union of rectangles in a specified
C     coordinate system.
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
C     DSKBRIEF
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dskdsc.inc'

      LOGICAL               GAP
      INTEGER               CORSYS
      INTEGER               NSIG
      INTEGER               NCOMP
      DOUBLE PRECISION      BDS1   ( 2, * )
      DOUBLE PRECISION      BDS2   ( 2, * )
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     GAP        I   Flag commanding display of summary or warning.
C     CORSYS     I   Coordinate system code.
C     NSIG       I   Number of significant digits.
C     NCOMP      I   Number of rectangles in gap region.
C     BDS1       I   Bounds of first coordinates of rectangles.
C     BDS2       I   Bounds of second coordinates of rectangles.
C
C$ Detailed_Input
C
C     GAP        is a logical flag indicating whether a gap summary or
C                warning is to be displayed. The summary is displayed
C                if and only if GAP is .TRUE. and coverage gaps exist.
C
C     CORSYS     is an integer code indicating the coordinate 
C                system in which the coverage gaps are represented.
C
C     NSIG       is the number of significant digits to display 
C                for floating point values.
C
C     NCOMP      is the number of coordinate rectangles---also 
C                called "components"---comprising the coverage
C                gap region.
C
C     BDS1       is an array of bounds of the first coordinates of the
C                gap rectangles. If the coordinate system is
C                latitudinal or planetodetic, the first coordinate is
C                longitude. If the coordinate system is rectangular,
C                the first coordinate is X.
C
C                The first coordinate of the Ith rectangle ranges from
C
C                   BDS1(1,I) to BDS1(2,I)
C
C                Units of longitude are radians; rectangular units are
C                km.
C
C     BDS2       is an array of bounds of the second coordinates of the
C                gap rectangles. If the coordinate system is
C                latitudinal or planetodetic, the second coordinate is
C                latitude. If the coordinate system is rectangular, the
C                second coordinate is Y.
C
C                The second coordinate of the Ith rectangle ranges from
C
C                   BDS2(1,I) to BDS2(2,I)
C
C                Units of latitude are radians; rectangular units are
C                km.
C
C$ Detailed_Output
C
C     None.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If an invalid coordinate system is encountered, this routine
C         signals the error SPICE(NOTSUPPORTED).
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine displays summary information for spatial coverage
C     gaps of a set of DSK segments. The gaps are represented as union
C     of rectangles in a specified coordinate system.
C
C     The display is written to standard output.
C
C     The coverage "gap region" of a set of DSK segments is that
C     region, within the rectangle defined by the global extrema of the
C     extents of the segments' first and second coordinates, where
C     there is no spatial coverage.
C
C     For example, if a set of DSK segments uses latitudinal
C     coordinates, the global extrema of the sets coordinates are the
C     minimum and maximum longitudes and latitudes, where the extrema
C     are taken over all segments in the set. Any point having
C     longitude between the set's longitude extrema and latitude
C     between the set's latitude extrema, such that point's longitude
C     and latitude are not within the coverage of any segment of the
C     set, lies in a gap. The coverage gap of the set is the union of
C     all such points. This region can be represented as a union of
C     rectangles.
C
C     In order for the concept of spatial coverage of a set of DSK
C     segments to make sense, the segments must have the common
C     attributes:
C
C        Body
C        Surface
C        Frame
C
C        Coordinate system
C
C           If the coordinate system is planetodetic, the
C           system parameters must match exactly.
C
C        Data type
C        Data class
C
C$ Examples
C
C     See usage in DSKBRIEF. 
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
C-    DSKBRIEF Version 1.0.0, 10-FEB-2017 (NJB)
C
C        Original version 07-OCT-2016 (NJB)
C
C-&
 
C$ Index_Entries
C
C     display spatial coverage gaps of dsk segment set
C
C-&

C
C     SPICELIB functions
C     
      DOUBLE PRECISION      DPR

      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               LNSIZE
      PARAMETER           ( LNSIZE = 132 )

      INTEGER               BUFSIZ
      PARAMETER           ( BUFSIZ = 1000 )

      INTEGER               NCOLS
      PARAMETER           ( NCOLS  = 4 )

C
C     Local variables
C
      CHARACTER*(LNSIZE)    CO1STR
      CHARACTER*(LNSIZE)    CO2STR
      CHARACTER*(1)         LABELS ( BUFSIZ )
      CHARACTER*(LNSIZE)    OUTLIN
      CHARACTER*(LNSIZE)    TABLE  ( BUFSIZ )

      DOUBLE PRECISION      SCALE
      DOUBLE PRECISION      VALUES ( NCOLS, BUFSIZ )

      INTEGER               FROM
      INTEGER               I
      INTEGER               NLINES
      INTEGER               REMAIN
      INTEGER               START1
      INTEGER               STARTS ( NCOLS )


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSPGAP' )


      DO I = 1, BUFSIZ
         LABELS(I) = ' '
         TABLE (I) = ' '
      END DO

C
C     Indent the first label by 2 spaces relative to the label above
C     it.
C
      START1 = 7

C
C     Set up coordinate system-specific column titles
C     and scale factors.
C
      IF (  GAP  .AND.  ( NCOMP .GT. 0 ) ) THEN

         CALL TOSTDO ( '    Coverage gaps:' )

         IF (     ( CORSYS .EQ. LATSYS ) 
     .       .OR. ( CORSYS .EQ. PDTSYS ) ) THEN

            SCALE  = DPR()

            CO1STR = 'Longitude range (deg)'
            CO2STR = 'Latitude range (deg)'

         ELSE IF ( CORSYS .EQ. RECSYS ) THEN

            SCALE  = 1.D0

            CO1STR = 'X range (km)'
            CO2STR = 'Y range (km)'

         ELSE
 
            CALL SETMSG ( 'Coordinate system # is not '
     .      //            'currently supported.'        )
            CALL ERRINT ( '#', CORSYS                   )
            CALL SIGERR ( 'SPICE(NOTSUPPORTED)'         )
            CALL CHKOUT ( 'DSPGAP'                      )
            RETURN

         END IF                       
C
C        Create the gap table for the first batch
C        of data. (We need to get the column positions
C        before we output the header, hence the order
C        of operations conducted here.)
C
         REMAIN = NCOMP

         NLINES = MIN( BUFSIZ, REMAIN )

         DO I = 1, NLINES

            VALUES(1,I) = BDS1(1,I) * SCALE
            VALUES(2,I) = BDS1(2,I) * SCALE
            VALUES(3,I) = BDS2(1,I) * SCALE
            VALUES(4,I) = BDS2(2,I) * SCALE

         END DO


         CALL CORTAB ( NLINES, LABELS, START1, NSIG,
     .                 NCOLS,  VALUES, STARTS, TABLE )
         IF ( FAILED() ) THEN
            CALL CHKOUT( 'DSPGAP' )
            RETURN
         END IF

         REMAIN = REMAIN - NLINES

C
C        Display the header.
C
         OUTLIN             = ' '
         OUTLIN(START1   :) = CO1STR
         OUTLIN(STARTS(3):) = CO2STR

         CALL TOSTDO ( OUTLIN )       

C
C        Display the  first batch of gap information.
C
         DO I = 1, NLINES
                  
            CALL TOSTDO ( TABLE(I) )

         END DO

C
C        Display the remaining gap information.
C
         FROM = NLINES + 1

         DO WHILE ( REMAIN .GT. 0 )

            NLINES = MIN( BUFSIZ, REMAIN )

            DO I = 1, NLINES

               VALUES(1,I) = BDS1(1,FROM) * SCALE
               VALUES(2,I) = BDS1(2,FROM) * SCALE
               VALUES(3,I) = BDS2(1,FROM) * SCALE
               VALUES(4,I) = BDS2(2,FROM) * SCALE

            END DO

            FROM = FROM + NLINES

            CALL CORTAB ( NLINES, LABELS, START1, NSIG,
     .                    NCOLS,  VALUES, STARTS, TABLE )
  
            IF ( FAILED() ) THEN
               CALL CHKOUT( 'DSPGAP' )
               RETURN
            END IF

            DO I = 1, NLINES

               CALL TOSTDO ( TABLE(I) )

            END DO

            REMAIN = REMAIN - NLINES

         END DO       
      

      ELSE IF ( NCOMP .GT. 0 ) THEN
C
C        Gap display is not enabled, but there are gaps.
C
         CALL TOSTDO ( '    ***Coverage has gaps. Use the -gaps ' 
     .   //            'option to display them.***'             )

      END IF
       
      CALL CHKOUT ( 'DSPGAP' )
      RETURN
      END
