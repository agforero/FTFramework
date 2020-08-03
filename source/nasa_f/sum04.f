C$Procedure   SUM04 ( DSKBRIEF, summarize type 4 segment )
 
      SUBROUTINE SUM04 ( HANDLE, DLADSC, NSIG )
 
C$ Abstract
C
C     Display type 4-specific summary of contents of a SPICE 
C     Digital Shape Kernel (DSK).
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
 
      INCLUDE 'dla.inc'
      INCLUDE 'dskdsc.inc'
      INCLUDE 'dsk04.inc'

      INTEGER               HANDLE
      INTEGER               DLADSC ( DLADSZ )
      INTEGER               NSIG

 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     HANDLE     I   Handle of DSK file.
C     DLADSC     I   DLA descriptor of segment.
C     NSIG       I   Number of significant digits in floating point
C                    output.
C
C$ Detailed_Input
C
C     HANDLE     is the handle of a DSK file containing a segment
C                to be summarized.
C
C     DLADSC     is the DLA descriptor of a segment to be summarized.
C
C     NSIG       is the number of significant digits in floating point
C                numeric output. 
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
C     1)  If an unrecognized coordinate system is encountered, this
C         routine signals the error SPICE(NOTSUPPORTED).
C
C     2)  If a map projection not compatible with the input coordinate
C         system is encountered, this routine signals the error
C         SPICE(NOTSUPPORTED).
C
C$ Files
C
C     See the input HANDLE.
C
C$ Particulars
C
C     This routine displays detailed summary information for a
C     specified type 4 DSK segment. The display is written to 
C     standard output.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     1) The expected range of NSIG is 6:17.
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
C-    DSKBRIEF Version 1.0.0, 06-OCT-2016 (NJB)
C
C        Version 4.0.0 04-OCT-2016 (NJB)
C        Version 4.0.0 04-OCT-2013 (NJB)
C        Version 3.2.0 14-NOV-2012 (NJB)
C        Version 3.1.0 11-NOV-2012 (NJB)
C        Version 3.0.0 04-OCT-2012 (NJB)
C        Version 2.0.0 22-SEP-2012 (NJB)
C        Version 1.0.0 17-AUG-2012 (NJB)
C
C-&
 
C$ Index_Entries
C
C     summarize type 4 dsk segment
C
C-&

C
C     SPICELIB functions
C     
      DOUBLE PRECISION      DPR

C
C     Local parameters
C 
      INTEGER               LNSIZE
      PARAMETER           ( LNSIZE = 132 )

C
C     Local variables
C
      CHARACTER*(LNSIZE)    ACCSTR
      CHARACTER*(LNSIZE)    FMT1
      CHARACTER*(LNSIZE)    ITPSTR
      CHARACTER*(LNSIZE)    OUTLIN
      CHARACTER*(LNSIZE)    XSTR

      DOUBLE PRECISION      ACCDSC ( MAXDSZ )
      DOUBLE PRECISION      CENT1
      DOUBLE PRECISION      CENT2
      DOUBLE PRECISION      CO1PSZ
      DOUBLE PRECISION      CO2PSZ
      DOUBLE PRECISION      DPBUFF ( 1 )
      DOUBLE PRECISION      DPNFMT
      DOUBLE PRECISION      DSKDSC ( DSKDSZ )
      DOUBLE PRECISION      HSCALE
      DOUBLE PRECISION      ITPDSC ( MAXDSZ )
      DOUBLE PRECISION      NULVAL
      DOUBLE PRECISION      NVDSC  ( MAXDSZ )
      DOUBLE PRECISION      PRJDSC ( MAXDSZ )
      DOUBLE PRECISION      REFDSC ( MAXDSZ )
      DOUBLE PRECISION      XDSC   ( MAXDSZ )

      INTEGER               B
      INTEGER               CORSYS
      INTEGER               DSKFMT
      INTEGER               E
      INTEGER               GRDIMS ( 2, DMAX )
      INTEGER               I
      INTEGER               IBASE
      INTEGER               N
      INTEGER               NC
      INTEGER               NF
      INTEGER               NGDIM
      INTEGER               NR
      INTEGER               NUMFMT
      INTEGER               NW

      LOGICAL               NULLOK     

C
C     Saved variables
C
 
C
C     Initial values
C     

      CALL CHKIN ( 'SUM04' )

C
C     Create d.p. format string.
C
      NW   = NSIG + 7
      NF   = NSIG - 1

      FMT1 = '(1PE@W.@F)'

      CALL REPMI ( FMT1, '@W', NW, FMT1 )
      CALL REPMI ( FMT1, '@F', NF, FMT1 )

C
C     Display type 4 parameters.
C
      CALL TOSTDO ( ' ' )
 
      CALL TOSTDO ( 'Type 4 parameters' )
      CALL TOSTDO ( '-----------------' )
 
C
C     Check the DSK format version.
C
      CALL DSKD04 ( HANDLE, DLADSC, KWFMTV, 1, 1, N, DPBUFF )

      DSKFMT = NINT ( DPBUFF(1) )

      IF ( DSKFMT .NE. 3 ) THEN

         CALL SETMSG ( 'DSK format version was expected to '
     .   //            'be 3 but was #. Only version 3 is '
     .   //            'supported by this software version.' )
         CALL ERRINT ( '#', DSKFMT                           )
         CALL SIGERR( 'SPICE(VERSIONMISMATCH)'               )
         CALL CHKOUT( 'SUM04'                                )
         RETURN

      END IF

      CALL DSKB04 ( HANDLE,  DLADSC,  MAXDSZ,  MAXDSZ,
     .              MAXDSZ,  NR,      NC,      CO1PSZ,
     .              CO2PSZ,  CENT1,   CENT2,   NULLOK,
     .              NULVAL,  ITPDSC,  XDSC,    NVDSC  )

C
C     Get and display the numeric data format.
C
      OUTLIN = '   Height data numeric format:         #'

      CALL DSKD04 ( HANDLE, DLADSC, KWNUMF, 1, 1, N, DPNFMT )

      NUMFMT = NINT( DPNFMT )

      IF ( NUMFMT .EQ. NUMI16 ) THEN

         CALL REPMC  ( OUTLIN, '#', '16-bit integer', OUTLIN )
         CALL TOSTDO ( OUTLIN )

      ELSE IF ( NUMFMT .EQ. NUMI32 ) THEN

         CALL REPMC  ( OUTLIN, '#', '32-bit integer', OUTLIN )
         CALL TOSTDO ( OUTLIN )

      ELSE

         CALL REPMC  ( OUTLIN, '#', 'Not supported', OUTLIN )
         CALL TOSTDO ( OUTLIN )

      END IF
      
C
C     Get the map projection descriptor; display the projection.
C
      CALL DSKD04 ( HANDLE, DLADSC, KWPRJD, 1, MAXDSZ, N, PRJDSC )

      OUTLIN = '   Map projection:                     #'

      IF ( PRJDSC( PRJXCD ) .EQ. PRJEQR ) THEN

         CALL REPMC  ( OUTLIN, '#', 'Equirectangular', OUTLIN )
         CALL TOSTDO ( OUTLIN )

      ELSE IF ( PRJDSC( PRJXCD ) .EQ. PRJSTE ) THEN

         CALL REPMC  ( OUTLIN, '#', 'Stereographic',    OUTLIN )
         CALL TOSTDO ( OUTLIN )

      ELSE
          CALL REPMC ( OUTLIN, '#', 'Not supported',   OUTLIN )   
      END IF

C
C     Get the reference surface; display the surface parameters.
C
      CALL DSKD04 ( HANDLE, DLADSC, KWREFD, 1, MAXDSZ, N, REFDSC )

      OUTLIN = '   Reference surface:                  #'

      IF ( REFDSC( REFXCD ) .EQ. REFINH ) THEN

         CALL REPMC  ( OUTLIN, '#', 
     .                 'From coordinate system (above)', 
     .                 OUTLIN                                     )
         CALL TOSTDO ( OUTLIN )

      ELSE IF ( REFDSC( PRJXCD ) .EQ. REFSPH ) THEN

         CALL REPMC  ( OUTLIN, '#', 'Sphere',               OUTLIN )
         CALL TOSTDO ( OUTLIN )

         OUTLIN = '   Radius (km):                        #'
         CALL REPMD  ( OUTLIN, '#', REFDSC( REFXPA ), 9,    OUTLIN )
         CALL TOSTDO ( OUTLIN )


      ELSE
          CALL REPMC ( OUTLIN, '#', 'Not supported',   OUTLIN )   
      END IF


C
C     Get the DSK descriptor.
C 
      CALL DSKGD( HANDLE, DLADSC, DSKDSC )

C
C     Get the coordinate system.
C
      CORSYS = NINT( DSKDSC(SYSIDX) )

C
C     Show row and column counts.
C
      OUTLIN = '   Number of pixel grid rows:          #'
      CALL REPMI  ( OUTLIN, '#', NR, OUTLIN )
      CALL TOSTDO ( OUTLIN )
 
      OUTLIN = '   Number of pixel grid columns:       #'
      CALL REPMI  ( OUTLIN, '#', NC, OUTLIN )
      CALL TOSTDO ( OUTLIN )

C
C     Show nested grid dimensions if there are at least two levels.
C
C
C     Fetch the number of nested grid dimensions and the
C     dimensions themselves.
C
      IBASE = DLADSC( IBSIDX )

      B = IBASE + IXNDIM

      CALL DASRDI ( HANDLE, B, B, NGDIM )

      OUTLIN = '   Integer grid nesting levels:        #'

      CALL REPMI  ( OUTLIN, '#', NGDIM, OUTLIN )
      CALL TOSTDO ( OUTLIN )


      IF ( NGDIM .GE. 2 ) THEN

         B = IBASE + IXGDIM
         E = B     + 2*NGDIM - 1

         CALL DASRDI ( HANDLE, B, E, GRDIMS )

         DO I = 1, NGDIM

            OUTLIN = '      Grid dimensions at level #:      '
     .      //       '# rows x # columns'

            CALL REPMI  ( OUTLIN, '#', I,           OUTLIN )
            CALL REPMI  ( OUTLIN, '#', GRDIMS(1,I), OUTLIN )
            CALL REPMI  ( OUTLIN, '#', GRDIMS(2,I), OUTLIN )
            CALL TOSTDO ( OUTLIN )

         END DO
 
      END IF



 
C
C     Show pixel dimensions and grid center coordinates.
C
C     Note: the code below is applicable only to the
C     equirectangular projection.
C
      IF ( PRJDSC( PRJXCD ) .EQ. PRJEQR ) THEN

         IF ( CORSYS .EQ. LATSYS ) THEN
C
C           Show coordinate 1 pixel dimension.
C            
            OUTLIN = '   Longitude pixel dimension (deg):    #'
 
            CALL REPMF ( OUTLIN, '#', CO1PSZ*DPR(), NSIG, 'E', OUTLIN )
            CALL TOSTDO( OUTLIN )

C
C           Show coordinate 2 pixel dimension. 
C
            OUTLIN = '   Latitude pixel dimension (deg):     #'
            CALL REPMF ( OUTLIN, '#', CO2PSZ*DPR(), NSIG, 'E', OUTLIN )
            CALL TOSTDO( OUTLIN )

C
C           Show grid center coordinates.
C
            OUTLIN = '   Pixel grid center longitude (deg):  #'
            CALL REPMF ( OUTLIN, '#', CENT1*DPR(), NSIG, 'F', OUTLIN )
            CALL TOSTDO( OUTLIN )

            OUTLIN = '   Pixel grid center latitude (deg):   #'
            CALL REPMF ( OUTLIN, '#', CENT2*DPR(), NSIG, 'F', OUTLIN )
            CALL TOSTDO( OUTLIN )


         ELSE IF ( CORSYS .EQ. PDTSYS ) THEN
C
C           Show coordinate 1 pixel dimension.
C
            OUTLIN = '   Longitude pixel dimension   (deg):  #'
            CALL REPMF ( OUTLIN, '#', CO1PSZ*DPR(), NSIG, 'E', OUTLIN )
            CALL TOSTDO( OUTLIN )

C
C           Show coordinate 2 pixel dimension. 
C
            OUTLIN = '   Latitude pixel dimension    (deg):  #'
            CALL REPMF ( OUTLIN, '#', CO2PSZ*DPR(), NSIG, 'E', OUTLIN )
            CALL TOSTDO( OUTLIN )

C
C           Show grid center coordinates.
C
            OUTLIN = '   Pixel grid center longitude (deg):  #'
            CALL REPMF ( OUTLIN, '#', CENT1*DPR(), NSIG, 'F', OUTLIN )
            CALL TOSTDO( OUTLIN )

            OUTLIN = '   Pixel grid center latitude  (deg):  #'
            CALL REPMF ( OUTLIN, '#', CENT2*DPR(), NSIG, 'F', OUTLIN )
            CALL TOSTDO( OUTLIN )

         ELSE

            CALL SETMSG( 'This coordinate system is not supported '
     .      //           'for the equirectangular projection. '
     .      //           'Coordinate system code: #'               )
            CALL ERRINT( '#', CORSYS                               )
            CALL SIGERR( 'SPICE(NOTSUPPORTED)'                     )

         END IF


      ELSE IF (  PRJDSC( PRJXCD ) .EQ. PRJSTE  ) THEN

         IF ( CORSYS .EQ. RECSYS ) THEN
C
C           Show coordinate 1 pixel dimension.
C
            OUTLIN = '   X pixel dimension (km):             #'
            CALL REPMF ( OUTLIN, '#', CO1PSZ, NSIG, 'E', OUTLIN )
            CALL TOSTDO( OUTLIN )

C
C           Show coordinate 2 pixel dimension. 
C
            OUTLIN = '   Y pixel dimension (km):             #'
            CALL REPMF ( OUTLIN, '#', CO2PSZ, NSIG, 'E', OUTLIN )
            CALL TOSTDO( OUTLIN )

C
C           Show grid center coordinates.
C
            OUTLIN = '   Pixel grid center X (km):           #'
            CALL REPMF ( OUTLIN, '#', CENT1, NSIG, 'F', OUTLIN )
            CALL TOSTDO( OUTLIN )

            OUTLIN = '   Pixel grid center Y (km):           #'
            CALL REPMF ( OUTLIN, '#', CENT2, NSIG, 'F', OUTLIN )
            CALL TOSTDO( OUTLIN )

         ELSE

            CALL SETMSG( 'This coordinate system is not supported '
     .      //           'for the projection having code #. '
     .      //           'Coordinate system code: #'               )
            CALL ERRINT( '#', NINT( PRJDSC( PRJXCD ) )             )
            CALL ERRINT( '#', CORSYS                               )
            CALL SIGERR( 'SPICE(NOTSUPPORTED)'                     )
            CALL CHKOUT( 'SUM04'                                   )
            RETURN

         END IF

      ELSE

         CALL SETMSG ( 'Unrecognized coordinate system code: #' )
         CALL ERRINT ( '#', CORSYS                              )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                    )
         CALL CHKOUT ( 'SUM04' )
         RETURN

      END IF
       

      OUTLIN = '   Null values allowed:                #'
      IF ( NULLOK ) THEN         
         CALL REPMC  ( OUTLIN, '#', 'Yes',  OUTLIN )
      ELSE
         CALL REPMC  ( OUTLIN, '#', 'No',   OUTLIN )
      END IF
      CALL TOSTDO ( OUTLIN )
  

      IF ( NULLOK ) THEN
         OUTLIN = '   Null value parameter:               #'
         CALL REPMI  ( OUTLIN, '#', NINT(NULVAL), OUTLIN )
         CALL TOSTDO ( OUTLIN )
      END IF


C
C     Show the height scale.
C
      CALL DSKD04 ( HANDLE, DLADSC, KWHSCL, 1, 1, N, HSCALE )

      OUTLIN = '   Height units in km:                 #'
      
      CALL REPMD  ( OUTLIN, '#', HSCALE, NSIG, OUTLIN )
      CALL TOSTDO ( OUTLIN )



      OUTLIN = '   Interpolation method:               #'

      IF ( ITPDSC(ITPXCD) .EQ. ITPNO ) THEN
         ITPSTR = 'None: return raw height data'
      ELSE IF  ( ITPDSC(ITPXCD) .EQ. ITPBIL ) THEN
         ITPSTR = 'Bilinear'
      ELSE IF  ( ITPDSC(ITPXCD) .EQ. ITPNN3 ) THEN
         ITPSTR = '3 nearest neighbor linear'
      ELSE
         CALL SETMSG( 'Bad interpolation code: #' )
         CALL ERRINT( '#', NINT(ITPDSC(ITPXCD))   )
         CALL SIGERR( 'BUG'                       )
         CALL CHKOUT( 'SUM04'                     )
         RETURN
      END IF

      CALL REPMC  ( OUTLIN, '#', ITPSTR, OUTLIN )
      CALL TOSTDO ( OUTLIN )


      OUTLIN = '   Surface intercept method:           #'

      IF ( XDSC(XXCD) .EQ. XCASTP ) THEN
         XSTR = 'Constant angular step'
      ELSE IF  ( XDSC(XXCD) .EQ. XCPSTP ) THEN
         XSTR = 'Constant projected step'
      ELSE IF  ( XDSC(XXCD) .EQ. XCELBD ) THEN
         XSTR = 'Cell boundary step'
      ELSE
         CALL SETMSG( 'Bad intercept code: #' )
         CALL ERRINT( '#', NINT(XDSC(XXCD))   )
         CALL SIGERR( 'BUG'                       )
         CALL CHKOUT( 'SUM04'                     )
         RETURN
      END IF

      CALL REPMC  ( OUTLIN, '#', XSTR, OUTLIN )
      CALL TOSTDO ( OUTLIN )

      IF ( XDSC(XXCD) .EQ. XCASTP ) THEN

         OUTLIN = '      Step size                (deg):  #'
         CALL REPMF ( OUTLIN, '#', XDSC(XXSTEP)*DPR(), 
     .                NSIG,   'E', OUTLIN              )
         CALL TOSTDO( OUTLIN )

         OUTLIN = '      Convergence tolerance      (m):  #'
         CALL REPMF  ( OUTLIN, '#', XDSC(XXTOL)*1.D3,   
     .                 NSIG,   'E', OUTLIN             )
         CALL TOSTDO ( OUTLIN )

      ELSE IF ( XDSC(XXCD) .EQ. XCPSTP ) THEN

         OUTLIN = '      Step size(km):                #'
         CALL REPMF  ( OUTLIN, '#', XDSC(XXSTEP)*DPR(), 
     .                 NSIG,   'E', OUTLIN             )
         CALL TOSTDO ( OUTLIN )

         OUTLIN = '      Convergence tolerance (km):   #'
         CALL REPMF  ( OUTLIN, '#', XDSC(XXTOL),  
     .                 NSIG,   'E', OUTLIN             )
         CALL TOSTDO ( OUTLIN )

      ELSE

         CALL SETMSG( 'Bad intercept code: #' )
         CALL ERRINT( '#', NINT(XDSC(XXCD))   )
         CALL SIGERR( 'BUG'                       )
         CALL CHKOUT( 'SUM04'                     )
         RETURN

      END IF 

      
C
C     Display surface intercept acceleration information:
C
      CALL DSKD04 ( HANDLE, DLADSC, KWACCD, 1, MAXDSZ, N, ACCDSC )

      OUTLIN = '      Acceleration algorithm:          #'      

      IF ( ACCDSC(AXCD) .EQ. ACCNO ) THEN
         ACCSTR = 'Not enabled'

         CALL REPMC  ( OUTLIN, '#', ACCSTR, OUTLIN )
         CALL TOSTDO ( OUTLIN )


      ELSE IF  ( ACCDSC(AXCD) .EQ. ACCBIG ) THEN

         ACCSTR = 'Large initial step'

         CALL REPMC  ( OUTLIN, '#', ACCSTR, OUTLIN )
         CALL TOSTDO ( OUTLIN )

         OUTLIN = '         Maximum magnitude of height'

         CALL TOSTDO ( OUTLIN )
         OUTLIN = '         deltas between adjacent pixels:'
         CALL TOSTDO ( OUTLIN )

         OUTLIN = '            In COORD1 direction (km):  #'

         IF ( CORSYS .EQ. PDTSYS ) THEN
            CALL REPMC  ( OUTLIN, 'COORD1', 'longitude', OUTLIN )
         END IF
         CALL REPMD  ( OUTLIN, '#', ACCDSC(AXMXD1)*HSCALE, 
     .                 NSIG,   OUTLIN                     )
         CALL TOSTDO ( OUTLIN )

         OUTLIN = '            In COORD2 direction  (km):  #'
         IF ( CORSYS .EQ. PDTSYS ) THEN
            CALL REPMC  ( OUTLIN, 'COORD2', 'latitude', OUTLIN )
         END IF
         CALL REPMD  ( OUTLIN, '#', ACCDSC(AXMXD2)*HSCALE, 
     .                 NSIG,   OUTLIN                     )
         CALL TOSTDO ( OUTLIN )


      ELSE
         CALL SETMSG( 'Bad accleration code: #' )
         CALL ERRINT( '#', NINT(ACCDSC(AXCD))   )
         CALL SIGERR( 'BUG'                       )
         CALL CHKOUT( 'SUM04'                     )
         RETURN
      END IF



      OUTLIN = '   Normal vector calculation method:   #'

      IF  ( NVDSC(NVXCD) .EQ. NV1PAR ) THEN

         XSTR = 'Defined by partial derivatives'

      ELSE IF  ( NVDSC(NVXCD) .EQ. NVBCPT ) THEN
         XSTR = 'Bilinear compatible'
      ELSE
         CALL SETMSG( 'Bad normal vector code: #' )
         CALL ERRINT( '#', NINT(NVDSC(NVXCD))     )
         CALL SIGERR( 'BUG'                       )
         CALL CHKOUT( 'SUM04'                     )
         RETURN
      END IF

      CALL REPMC  ( OUTLIN, '#', XSTR, OUTLIN )
      CALL TOSTDO ( OUTLIN )

      CALL CHKOUT ( 'SUM04' )
      RETURN
      END



