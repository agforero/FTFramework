C$Procedure CORTAB ( DSKBRIEF, display coordinate table )

      SUBROUTINE CORTAB ( N,     LABELS, START1, NSIG,  
     .                    NCOLS, VALUES, STARTS, TABLE )

C$ Abstract
C
C     Display a table of coordinates or other double precision
C     values.
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

      INTEGER               N
      CHARACTER*(*)         LABELS ( * )
      INTEGER               START1
      INTEGER               NSIG
      INTEGER               NCOLS
      DOUBLE PRECISION      VALUES ( NCOLS, N )
      INTEGER               STARTS ( NCOLS )
      CHARACTER*(*)         TABLE  ( * )

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     N          I   is the number of table rows.
C     LABELS     I   is an array of row labels.
C     START1     I   is the position of the first data column.
C     NSIG       I   is the number of significant digits to display.
C     NCOLS      I   is the number of data columns.
C     VALUES     I   is an array of data values.
C     STARTS     O   is an array of data column positions.
C     TABLE      O   is the output table.
C
C$ Detailed_Input
C
C     N          is the number of rows in the table to be created.
C
C     LABELS     is an array of labels for the output rows. The labels
C                are located at the left side of the table.
C
C     START1     is the character position of the start of the first
C                data column in the table.
C
C     NSIG       is the number of significant digits to display for
C                the output values. NSIG ranges from 7 to 16.
C
C     NCOLS      is the number of data columns in the table.
C
C     VALUES     is an array of double precision values to be placed in 
C                the table. VALUES has dimensions
C
C                   NCOLS x N
C
C$ Detailed_Output
C
C     STARTS     is an array containing the starting character positions
C                of the output data columns. The first element of STARTS
C                is START1.
C
C     TABLE      is the output table. The Ith row of the table starts 
C                with 
C
C                   LABELS(I)
C
C                The values 
C
C                   VALUES(1,I), ..., VALUES(NCOLS,I)
C
C                follow, with the Jth value starting at index STARTS(J).
C
C                Values having magnitudes in the range
C
C                   [1, 1e6)
C
C                as well as zero, are expressed in fixed-point format.
C                Non-zero values having magnitudes outside this range
C                are expressed in scientific notation.
C
C                In each data column, the decimal points of the values
C                are aligned. 
C
C$ Parameters
C
C     ANGMRG     See the description in dsktol.inc.
C
C$ Exceptions
C
C     None.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine supports coverage gap display.
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
C-    DSKBRIEF Version 1.0.0, 15-MAR-2017 (NJB)
C     
C        Original version 18-NOV-2016 (NJB)
C
C-&




C
C     SPICELIB functions
C     
      LOGICAL               RETURN

C
C     Local parameters
C
      INTEGER               FMTLEN
      PARAMETER           ( FMTLEN = 30 )

C
C     Local variables
C
      CHARACTER*(FMTLEN)    FMT1

      DOUBLE PRECISION      DPVAL
      DOUBLE PRECISION      LOGV

      INTEGER               ADJ
      INTEGER               BCOL
      INTEGER               COLIX
      INTEGER               DPIX
      INTEGER               I
      INTEGER               MAXA
      INTEGER               MAXP
      INTEGER               MAXT
      INTEGER               POS1
      INTEGER               POS2
      INTEGER               WE
      INTEGER               WI
      INTEGER               WF
      INTEGER               WT


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'CORTAB' )

      DO I = 1, N
         TABLE(I) = LABELS(I)
      END DO


C
C     BCOL is the character index of the first column.
C
      STARTS(1) = START1
      BCOL      = STARTS(1)


      DO COLIX = 1, NCOLS
C
C        Find maximum character counts of integer and fractional parts
C        of the values in the current column.
C
C           MAXA is the maximum adjustment value for all rows.
C
C           MAXP is the maximum decimal point index for all rows.
C           The index is relative to the first character of the row.
C
C           MAXT is the maximum total field width for any value in
C           the current column, which has index COLIX.
C
         MAXA = 0
         MAXP = 0
         MAXT = 0

         DO I = 1, N

            DPVAL = VALUES(COLIX,I)

            IF (      ( ABS(DPVAL) .GE. 1.D0 ) 
     .          .AND. ( ABS(DPVAL) .LT. 1.D6 ) ) THEN
C
C              We can represent the value in fixed-point format.
C
C              Let WI be the width of the integer part; WF is the width
C              of the fractional part. WT is the total width of the
C              string representing the value.
C           
               LOGV = LOG10 ( ABS(DPVAL) )

C
C              Include room for the sign in WI.
C
               WI   = INT( LOGV ) + 1
               WF   = MAX( 0,  NSIG - WI )
               WT   = WI + WF + 2

               FMT1 = '(F@W.@F)'
               CALL REPMI ( FMT1, '@W', WT, FMT1 )
               CALL REPMI ( FMT1, '@F', WF, FMT1 )


            ELSE IF ( DPVAL .EQ. 0.D0 ) THEN

               WI   = 1
               WF   = MAX( 0,  NSIG - WI )
               WT   = WI + WF + 1

               FMT1 = '(F@W.@F)'
               CALL REPMI ( FMT1, '@W', WT, FMT1 )
               CALL REPMI ( FMT1, '@F', WF, FMT1 )

            ELSE
C
C              Use exponential notation. WF includes room for the 
C              exponent, which includes the 'E' symbol, a sign, and
C              three digits.
C
               WE   = 5
               WI   = 1
               WF   = MAX( 0,  NSIG - WI )
               WT   = WI   + WF + 2 + WE 

               FMT1 = '(1PE@W.@F)'
               CALL REPMI ( FMT1, '@W', WT, FMT1 )
               CALL REPMI ( FMT1, '@F', WF, FMT1 )

            END IF
C
C           Write the value to the output table, starting at
C           index BCOL.
C
            WRITE ( TABLE(I)(BCOL:), FMT1 ) DPVAL           
C
C           Find the offset of the decimal point within the string
C           value just written.
C
            DPIX = INDEX( TABLE(I)(BCOL:), '.' )

            MAXP = MAX( MAXP, DPIX )
            MAXT = MAX( MAXT, WT )

            IF ( DPVAL .EQ. 0.D0 ) THEN
C
C              Adjust the output string to have a leading zero if it's
C              needed. Fortran 77 says leading zeros are optional for
C              fixed-point representations of numbers having magnitude
C              less than 1. Zero qualifies.
C               
               IF ( DPIX .EQ. 1 ) THEN
                  
                  CALL PREFIX ( '0', 0, TABLE(I)(BCOL:BCOL-1+WT) )

               ELSE IF ( DPIX .GT. 1 ) THEN

                  POS1 = DPIX - 1
                  POS2 = BCOL - 1 + POS1

                  IF ( TABLE(I)(BCOL:POS2) .EQ. ' ' ) THEN

                     TABLE(I)(POS2:POS2) = '0'

                  END IF

               ELSE
C
C                 Backstop case.
C
                  CALL SETMSG ( 'No decimal point found in string '
     .            //            'representing zero. I = #; '
     .            //            'TABLE(I) = #.'                   )
                  CALL ERRINT ( '#', I                            )
                  CALL ERRCH  ( '#', TABLE(I)                     )
                  CALL SIGERR ( 'SPICE(BUG)'                      )
                  CALL CHKOUT ( 'CORTAB'                          )
                  RETURN

               END IF

            END IF


         END DO

C
C        Align the values just written. 
C
         DO I = 1, N

            DPIX = INDEX( TABLE(I)(BCOL:), '.' )

            ADJ  = MAXP - DPIX

            CALL SHIFTC ( TABLE(I)(BCOL:), 'R',     ADJ,
     .                    ' ',             TABLE(I)(BCOL:) )
            
            MAXA = MAX( MAXA, ADJ )

         END DO

         BCOL            = BCOL + MAXT + MAXA + 2

         IF ( COLIX .EQ. 2 ) THEN
            BCOL = BCOL + 4
         END IF

C
C        Set the start value for the next column.
C
         IF ( COLIX .LT. NCOLS ) THEN
            STARTS(COLIX+1) = BCOL
         END IF

      END DO

      CALL CHKOUT ( 'CORTAB' )
      RETURN
      END
