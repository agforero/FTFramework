C$Procedure DSPDSC ( Display DSK segment descriptor )

      SUBROUTINE DSPDSC ( DSKDSC, N, ITEMS, NSIG ) 

C$ Abstract
C
C     Display a specified set of attributes from a DSK segment
C     descriptor.
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

      INCLUDE 'dskdsc.inc'
      INCLUDE 'srftrn.inc'

      DOUBLE PRECISION      DSKDSC ( DSKDSZ )
      INTEGER               N
      INTEGER               ITEMS  ( * )
      INTEGER               NSIG

C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     DSKDSC     I   DSK segment descriptor.
C     N          I   Number of descriptor items to display.
C     ITEMS      I   Array of descriptor item specifiers.
C     NSIG       I   Number of significant digits in floating point
C                    output.
C
C$ Detailed_Input
C
C     DSKDSC     is a DSK segment descriptor for which summary
C                information is to be displayed.
C
C     N          is the number of DSK segment attribute specifiers
C                in the ITEMS array.
C
C     ITEMS      is an array of DSK segment attribute specifiers.
C                These specifiers identify the descriptor information
C                to display. Information is displayed in the order
C                of the corresponding elements of ITEMS.
C
C                Each specifier is the index in the DSK descriptor
C                of an attribute. For example, if ITEMS(1) is
C                set to the value
C
C                   FRMIDX
C
C                then the first item displayed is the DSK segment's
C                reference frame.
C
C                Pairs of coordinate bounds and time bounds are
C                indicated by the index for the lower bound alone. For
C                example, both the lower and upper latitude bounds are
C                displayed as the Nth item if ITEMS(N) is set to MN2IDX
C                and if the coordinate system is latitudinal.
C
C                When the coordinate system name is displayed, any
C                associated coordinate system parameters are displayed
C                as well.
C
C     NSIG       is the number of significant digits in floating point
C                numeric output. The range of NSIG is 6:17.
C
C$ Detailed_Output
C
C     None. This routine operates by side effects.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the coordinate system descriptor in DSKDSC is not
C         recognized, the error SPICE(NOTSUPPORTED) is signaled.
C
C     2)  If NSIG is outside of the range 6:17, it is replaced with the
C         closest value in this range.
C
C$ Files
C
C     None.
C
C$ Particulars
C
C     This routine writes to standard output.
C
C     ID codes of bodies, surfaces, and reference frames are
C     displayed along with the corresponding names, whenever
C     the required ID-name mappings are available.
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
C        Previous version 04-OCT-2016 (NJB)
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
      CHARACTER*(*)         FMT1
      PARAMETER           ( FMT1 = '(1PE24.16)' )
      INTEGER               LNSIZE
      PARAMETER           ( LNSIZE = 132 )

      INTEGER               BDNMLN
      PARAMETER           ( BDNMLN = 36 )

      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )

      INTEGER               DTYPLN
      PARAMETER           ( DTYPLN = 80 )

      INTEGER               TIMLEN
      PARAMETER           ( TIMLEN = 35 )

      INTEGER               NTYPE
      PARAMETER           ( NTYPE  = 4 )

      INTEGER               VC1IX
      PARAMETER           ( VC1IX  = 38 )
      
C
C     Local variables
C
      CHARACTER*(TIMLEN)    BEGTIM
      CHARACTER*(BDNMLN)    BODY
      CHARACTER*(TIMLEN)    ENDTIM
      CHARACTER*(FRNMLN)    FRAME
      CHARACTER*(LNSIZE)    LABELS ( 3 )
      CHARACTER*(LNSIZE)    OUTLIN
      CHARACTER*(SFNMLN)    SRFNAM
      CHARACTER*(LNSIZE)    TABLE  ( 3 )
      CHARACTER*(DTYPLN)    TYPLST ( NTYPE )
      CHARACTER*(DTYPLN)    TYPNAM

      DOUBLE PRECISION      F
      DOUBLE PRECISION      MAXALT
      DOUBLE PRECISION      MAXLAT
      DOUBLE PRECISION      MAXLON
      DOUBLE PRECISION      MAXRAD
      DOUBLE PRECISION      MAXX
      DOUBLE PRECISION      MAXY
      DOUBLE PRECISION      MAXZ
      DOUBLE PRECISION      MINALT
      DOUBLE PRECISION      MINLAT
      DOUBLE PRECISION      MINLON
      DOUBLE PRECISION      MINRAD
      DOUBLE PRECISION      MINX
      DOUBLE PRECISION      MINY
      DOUBLE PRECISION      MINZ
      DOUBLE PRECISION      VALCOL ( 3 )
      DOUBLE PRECISION      VALUES ( 2, 3 )

      INTEGER               BODYID
      INTEGER               CORSYS
      INTEGER               DCLASS
      INTEGER               FRMCDE
      INTEGER               I
      INTEGER               J
      INTEGER               SRFACE
      INTEGER               STARTS ( 3 )
      INTEGER               TYPCDE

      LOGICAL               FND
      LOGICAL               ISNAME
      

C
C     Saved variables
C
      SAVE                  TYPLST

C
C     Initial values
C     
      DATA TYPLST /  '<Not implemented>',
     .               'Shape model using triangular plates',
     .               '<Not implemented>',
     .               'Packed integer DEM'                 /



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSPDSC' )

C
C     Display items in the order they're listed.
C
      DO I = 1, N

         IF ( ITEMS(I) .EQ. CTRIDX ) THEN
C
C           For historical reasons, the index of the body ID
C           is named CTRIDX.
C
            BODYID = NINT( DSKDSC(CTRIDX) )
C
C           Show body ID.
C
            CALL BODC2N ( BODYID, BODY, FND )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'DSPDSC' )
               RETURN
            END IF

            IF ( .NOT. FND ) THEN
               BODY = 'Name not available'
            END IF

            OUTLIN = 'Body:                               # (#)'

            CALL REPMI  ( OUTLIN, '#', BODYID, OUTLIN )
            CALL REPMC  ( OUTLIN, '#', BODY,   OUTLIN )
            CALL TOSTDO ( OUTLIN )

         END IF


         IF ( ITEMS(I) .EQ. SRFIDX ) THEN
C
C           Display the surface. 
C
            SRFACE = NINT( DSKDSC(SRFIDX) )
            BODYID = NINT( DSKDSC(CTRIDX) )

            CALL SRFC2S ( SRFACE, BODYID, SRFNAM, ISNAME )

            IF ( .NOT. ISNAME ) THEN
               SRFNAM = '  Name not available'
            END IF

            OUTLIN = '  Surface:                          # (#)'

            CALL REPMI  ( OUTLIN, '#', SRFACE, OUTLIN )
            CALL REPMC  ( OUTLIN, '#', SRFNAM, OUTLIN )
            CALL TOSTDO ( OUTLIN )

         END IF


         IF ( ITEMS(I) .EQ. FRMIDX ) THEN
C
C           Display the reference frame.
C
            FRMCDE = NINT( DSKDSC(FRMIDX ) )

            CALL FRMNAM ( FRMCDE, FRAME )

            IF ( FRAME .NE. ' ' ) THEN

               OUTLIN = '  Reference frame:                  #'
               CALL REPMC  ( OUTLIN, '#', FRAME,  OUTLIN )
            ELSE
               OUTLIN = '  Reference frame name N/A; ID code: #'
               CALL REPMI  ( OUTLIN, '#', FRMCDE, OUTLIN )
            END IF

            CALL TOSTDO ( OUTLIN )

         END IF


         IF ( ITEMS(I) .EQ. TYPIDX ) THEN
C
C           Display the data type.
C
            TYPCDE = NINT( DSKDSC(TYPIDX) ) 

            IF ( ( TYPCDE .GT. 0 ) .AND. ( TYPCDE .LE. NTYPE ) ) THEN

               TYPNAM = TYPLST(TYPCDE) 
            ELSE
               TYPNAM = '  Data type description not available'
            END IF

            OUTLIN = '  Data type:                        # (#)'
            CALL REPMI  ( OUTLIN, '#', TYPCDE,    OUTLIN )
            CALL REPMC  ( OUTLIN, '#', TYPNAM, OUTLIN )
            CALL TOSTDO ( OUTLIN )

         END IF


         IF ( ITEMS(I) .EQ. CLSIDX ) THEN
            
            DCLASS = NINT( DSKDSC(CLSIDX) )
C
C           Display the data class.
C
            IF ( DCLASS .EQ. SVFCLS ) THEN

               OUTLIN = '  Data class:                       '
     .         //       '1 (Single-valued surface)'

            ELSE IF ( DCLASS .EQ. GENCLS ) THEN
               OUTLIN = '  Data class:                       '
     .         //       '2 (General surface)'

            ELSE
               OUTLIN = '  Data class:                       '
     .         //       'unknown'
            END IF

            CALL TOSTDO ( OUTLIN )
 
         END IF



         CORSYS = NINT(DSKDSC(SYSIDX))
        

         IF ( ITEMS(I) .EQ. SYSIDX ) THEN
C
C           Display the coordinate system and coordinate system
C           parameters, if applicable.
C
C           Display the coordinate bounds as well.
C
            OUTLIN = '  Coordinate system:                #'

            IF ( CORSYS .EQ. LATSYS ) THEN
C
C              The system is LATITUDINAL.
C           
               CALL REPMC ( OUTLIN, '#', 
     .                      'Planetocentric Latitudinal', OUTLIN )
               CALL TOSTDO( OUTLIN )


            ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
C
C              The system is PLANETODETIC.
C           
               CALL REPMC ( OUTLIN, '#', 'Planetodetic', OUTLIN )
               CALL TOSTDO( OUTLIN )

               LABELS(1)   =  '   Equatorial radius (km):'
               LABELS(2)   =  '   Polar radius      (km):'
               LABELS(3)   =  '   Flattening coefficient:'

               VALCOL(1) = DSKDSC(PARIDX)
               F         = DSKDSC(PARIDX+1)
               VALCOL(2) = (1.D0 - F) * VALCOL(1)
               VALCOL(3) = F

               CALL CORTAB ( 3, LABELS, VC1IX-1,  NSIG, 
     .                       1, VALCOL, STARTS,   TABLE   )

               DO J = 1, 3
                  CALL TOSTDO( TABLE(J) )
               END DO


            ELSE IF ( CORSYS .EQ. RECSYS ) THEN
C
C              The system is RECTANGULAR.
C           
               CALL REPMC ( OUTLIN, '#', 'Rectangular',    OUTLIN )
               CALL TOSTDO ( OUTLIN )

            ELSE

               CALL SETMSG ( 'The coordinate system code # is not '
     .         //            'recognized.'                          )
               CALL ERRINT ( '#',  CORSYS                           )
               CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                  )
               CALL CHKOUT ( 'DSPDSC'                               )
               RETURN

            END IF

         END IF


         IF ( CORSYS .EQ. LATSYS ) THEN

 
            MINLON = DSKDSC(MN1IDX)
            MAXLON = DSKDSC(MX1IDX)
            MINLAT = DSKDSC(MN2IDX)
            MAXLAT = DSKDSC(MX2IDX)
            MINRAD = DSKDSC(MN3IDX)
            MAXRAD = DSKDSC(MX3IDX)
  

            LABELS(1) =  '    Min, max longitude  (deg):'
            LABELS(2) =  '    Min, max latitude   (deg):'
            LABELS(3) =  '    Min, max radius      (km):'

            VALUES(1,1) = MINLON * DPR()
            VALUES(1,2) = MINLAT * DPR()
            VALUES(1,3) = MINRAD
            VALUES(2,1) = MAXLON * DPR()
            VALUES(2,2) = MAXLAT * DPR()
            VALUES(2,3) = MAXRAD

            CALL CORTAB ( 3, LABELS, VC1IX,  NSIG, 
     .                    2, VALUES, STARTS, TABLE   )


            IF ( ITEMS(I) .EQ. MN1IDX ) THEN
C
C              Show longitude bounds.
C
               CALL TOSTDO ( TABLE(1) )

            END IF


            IF ( ITEMS(I) .EQ. MN2IDX ) THEN
C
C              Show latitude bounds.
C
               CALL TOSTDO ( TABLE(2) )

            END IF

            IF ( ITEMS(I) .EQ. MN3IDX ) THEN
C
C              Show radius bounds.
C
               CALL TOSTDO ( TABLE(3) )

           END IF


         ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
C
C           The system is PLANETODETIC.
C           
            MINLON = DSKDSC(MN1IDX)
            MAXLON = DSKDSC(MX1IDX)
            MINLAT = DSKDSC(MN2IDX)
            MAXLAT = DSKDSC(MX2IDX)
            MINALT = DSKDSC(MN3IDX)
            MAXALT = DSKDSC(MX3IDX)

            LABELS(1) =  '    Min, max longitude  (deg):'
            LABELS(2) =  '    Min, max latitude   (deg):'
            LABELS(3) =  '    Min, max altitude    (km):'

            VALUES(1,1) = MINLON * DPR()
            VALUES(1,2) = MINLAT * DPR()
            VALUES(1,3) = MINALT
            VALUES(2,1) = MAXLON * DPR()
            VALUES(2,2) = MAXLAT * DPR()
            VALUES(2,3) = MAXALT

            CALL CORTAB ( 3, LABELS, VC1IX,  NSIG, 
     .                    2, VALUES, STARTS, TABLE   )


            IF ( ITEMS(I) .EQ. MN1IDX ) THEN
C
C              Show longitude bounds.
C
               CALL TOSTDO ( TABLE(1) )

            END IF


            IF ( ITEMS(I) .EQ. MN2IDX ) THEN
C
C              Show latitude bounds.
C
               CALL TOSTDO ( TABLE(2) )

            END IF


            IF ( ITEMS(I) .EQ. MN3IDX ) THEN
C
C              Show altitude bounds.
C
               CALL TOSTDO ( TABLE(3) )

            END IF


         ELSE IF ( CORSYS .EQ. RECSYS ) THEN
C
C           The system is RECTANGULAR.
C           
            MINX = DSKDSC(MN1IDX)
            MAXX = DSKDSC(MX1IDX)
            MINY = DSKDSC(MN2IDX)
            MAXY = DSKDSC(MX2IDX)
            MINZ = DSKDSC(MN3IDX)
            MAXZ = DSKDSC(MX3IDX)

            LABELS(1) =  '   Min, max X coordinate (km):'
            LABELS(2) =  '   Min, max Y coordinate (km):'
            LABELS(3) =  '   Min, max Z coordinate (km):'

            VALUES(1,1) = MINX
            VALUES(1,2) = MINY
            VALUES(1,3) = MINZ
            VALUES(2,1) = MAXX
            VALUES(2,2) = MAXY
            VALUES(2,3) = MAXZ

            CALL CORTAB ( 3, LABELS, VC1IX,  NSIG, 
     .                    2, VALUES, STARTS, TABLE   )
            IF ( ITEMS(I) .EQ. MN1IDX ) THEN
C
C              Show X bounds.
C
               CALL TOSTDO ( TABLE(1) )

            END IF

            IF ( ITEMS(I) .EQ. MN2IDX ) THEN

               CALL TOSTDO ( TABLE(2) )

            END IF


            IF ( ITEMS(I) .EQ. MN3IDX ) THEN

               CALL TOSTDO ( TABLE(3) )

            END IF


         ELSE

            CALL SETMSG ( 'The coordinate system code # is not '
     .      //            'recognized.'                          )
            CALL ERRINT ( '#',  CORSYS                           )
            CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                  )
            CALL CHKOUT ( 'DSPDSC'                               )
            RETURN

         END IF


         IF ( ITEMS(I) .EQ. BTMIDX ) THEN
C
C           Display both the start and stop times.
C
            CALL ETCAL ( DSKDSC(BTMIDX), BEGTIM )
            CALL ETCAL ( DSKDSC(ETMIDX), ENDTIM )
            CALL SUFFIX ( 'TDB', 1, BEGTIM )
            CALL SUFFIX ( 'TDB', 1, ENDTIM )            

            OUTLIN = '  Start time:                       '//BEGTIM
            CALL TOSTDO ( OUTLIN )

            WRITE (OUTLIN, '(A,'//FMT1//')')
     .               '    Seconds past J2000 TDB:             ', 
     .               DSKDSC(BTMIDX)
            CALL TOSTDO ( OUTLIN )


            OUTLIN = '  Stop time:                        '//ENDTIM
            CALL TOSTDO ( OUTLIN )

            WRITE (OUTLIN, '(A,'//FMT1//')')
     .               '    Seconds past J2000 TDB:             ', 
     .               DSKDSC(ETMIDX)
            CALL TOSTDO ( OUTLIN )

         END IF

      END DO

      CALL CHKOUT ( 'DSPDSC' )
      RETURN
      END


