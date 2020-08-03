C$Procedure PRCSET ( Process setup file for MKDSK---umbrella routine )

      SUBROUTINE PRCSET ( SETUP,   INPUT,   OUTPUT,  CMTFIL, 
     .                    SURFID,  CENTID,  FRAME,   FIRST,
     .                    LAST,    DCLASS,  DTYPE,   AUNITS,
     .                    DUNITS,  CORSYS,  CORPAR,  MNCOR1,   
     .                    MXCOR1,  MNCOR2,  MXCOR2,  PLTYPE,  
     .                    VOXSCL,  CGRSCL,  WRAP,    MKNCAP,
     .                    MKSCAP,  ROWMAJ,  TOPDWN,  LEFTRT,
     .                    REFVAL,  HSCALE,  NCOLS,   NROWS,
     .                    LFTCOR,  TOPCOR,  COLSTP,  ROWSTP, 
     .                    MAKVPM                            )
  
C$ Abstract
C
C     Umbrella routine for MKDSK setup file processing.
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

      INCLUDE 'mkdsk.inc'
      INCLUDE 'dskdsc.inc'

      CHARACTER*(*)         SETUP
      CHARACTER*(*)         INPUT
      CHARACTER*(*)         OUTPUT
      CHARACTER*(*)         CMTFIL
      INTEGER               SURFID
      INTEGER               CENTID
      CHARACTER*(*)         FRAME
      DOUBLE PRECISION      FIRST
      DOUBLE PRECISION      LAST
      INTEGER               DCLASS
      INTEGER               DTYPE
      CHARACTER*(*)         AUNITS
      CHARACTER*(*)         DUNITS
      INTEGER               CORSYS
      DOUBLE PRECISION      CORPAR ( * )
      DOUBLE PRECISION      MNCOR1
      DOUBLE PRECISION      MXCOR1
      DOUBLE PRECISION      MNCOR2
      DOUBLE PRECISION      MXCOR2
      INTEGER               PLTYPE
      DOUBLE PRECISION      VOXSCL
      INTEGER               CGRSCL
      LOGICAL               WRAP
      LOGICAL               MKNCAP
      LOGICAL               MKSCAP
      LOGICAL               ROWMAJ
      LOGICAL               TOPDWN
      LOGICAL               LEFTRT
      DOUBLE PRECISION      REFVAL
      DOUBLE PRECISION      HSCALE
      INTEGER               NCOLS
      INTEGER               NROWS
      DOUBLE PRECISION      LFTCOR
      DOUBLE PRECISION      TOPCOR
      DOUBLE PRECISION      COLSTP
      DOUBLE PRECISION      ROWSTP
      LOGICAL               MAKVPM

C$ Brief_I/O
C
C     VARIABLE  I/O  Entry points
C     --------  ---  --------------------------------------------------
C     SETUP      I   GETSET
C     INPUT     I/O  GETSET
C     OUTPUT    I/O  GETSET
C     CMTFIL     O   GETSET
C
C
C$ Detailed_Input
C
C     See the entry points for a description of their inputs.
C
C$ Detailed_Output
C
C     See the entry points for a description of their outputs.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     If this routine is called directly, the error SPICE(BOGUSENTRY)
C     is signaled.  See the entry points for descriptions of 
C     exceptions specific to those routines.
C
C$ Files
C
C     This suite of routines reads and returns information from an
C     MKDSK setup file. See the MKDSK setup template for a list of
C     supported setup parameters.
C
C$ Particulars
C
C     The entry points in this package are
C
C        GETSET  {Get setup file information for MKDSK}
C        GETTYP  {Get segment data type}
C        GETGEN  {Get general DSK parameters}
C        GETP02  {Get type 2 parameters}
C
C     GETSET must be called before the other entry points may be called.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     This routine is intended for use only within the program 
C     MKDSK.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     B.V. Semenov    (JPL)
C
C$ Version
C
C-    MKDSK Version 3.0.0, 08-MAR-2017 (NJB) (BVS)
C
C        Updated to support automatic voxel scale setting.
C        Updated to support plate type 5 (height grid) input
C        format. 
C
C        Some error handling bugs were corrected.
C
C        Last update 19-JAN-2016 (NJB)
C
C           Updated to support surface name-ID translation
C           and new coordinate systems. Updated header.
C
C-    MKDSK Version 2.0.0, 29-JUN-2010 (NJB)
C
C        Updated shape and DSK file keywords.
C
C-    MKDSK Version 1.0.0, 15-APR-2010 (NJB)
C
C-&

C
C     SPICELIB functions
C
      INTEGER               ESRCHC

      LOGICAL               EQSTR
      LOGICAL               EXISTS
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      CHARACTER*(*)         KWFRAM
      PARAMETER           ( KWFRAM = 'REF_FRAME_NAME'       )

      CHARACTER*(*)         KWCSYS
      PARAMETER           ( KWCSYS = 'COORDINATE_SYSTEM'    )

      CHARACTER*(*)         KWUNIT
      PARAMETER           ( KWUNIT = 'INPUT_DATA_UNITS'     )

      CHARACTER*(*)         KWSTRT
      PARAMETER           ( KWSTRT = 'START_TIME'           )

      CHARACTER*(*)         KWINP 
      PARAMETER           ( KWINP  = 'INPUT_SHAPE_FILE'     )

      CHARACTER*(*)         KWLSK 
      PARAMETER           ( KWLSK  = 'LEAPSECONDS_FILE'     )

      CHARACTER*(*)         KWOUT
      PARAMETER           ( KWOUT  = 'OUTPUT_DSK_FILE'      )

      CHARACTER*(*)         KWCMNT
      PARAMETER           ( KWCMNT = 'COMMENT_FILE'         )

      CHARACTER*(*)         KWSTOP
      PARAMETER           ( KWSTOP = 'STOP_TIME'            )

      CHARACTER*(*)         KWSURF
      PARAMETER           ( KWSURF = 'SURFACE_NAME'         )

      CHARACTER*(*)         KWCENT
      PARAMETER           ( KWCENT = 'CENTER_NAME'          )

      CHARACTER*(*)         KWDTYP
      PARAMETER           ( KWDTYP = 'DATA_TYPE'            )

      CHARACTER*(*)         KWCLSS
      PARAMETER           ( KWCLSS = 'DATA_CLASS'           )

      CHARACTER*(*)         KWMNLA
      PARAMETER           ( KWMNLA = 'MINIMUM_LATITUDE'     )

      CHARACTER*(*)         KWMXLA
      PARAMETER           ( KWMXLA = 'MAXIMUM_LATITUDE'     )

      CHARACTER*(*)         KWMNLO
      PARAMETER           ( KWMNLO = 'MINIMUM_LONGITUDE'    )

      CHARACTER*(*)         KWMXLO
      PARAMETER           ( KWMXLO = 'MAXIMUM_LONGITUDE'    )

      CHARACTER*(*)         KWMNX
      PARAMETER           ( KWMNX  = 'MINIMUM_X'    )

      CHARACTER*(*)         KWMXX
      PARAMETER           ( KWMXX  = 'MAXIMUM_X'    )

      CHARACTER*(*)         KWMNY
      PARAMETER           ( KWMNY  = 'MINIMUM_Y'    )

      CHARACTER*(*)         KWMXY
      PARAMETER           ( KWMXY  = 'MAXIMUM_Y'    )

      CHARACTER*(*)         KWRE
      PARAMETER           ( KWRE   = 'EQUATORIAL_RADIUS' )

      CHARACTER*(*)         KWRP
      PARAMETER           ( KWRP   = 'POLAR_RADIUS' )

      CHARACTER*(*)         KWPTYP
      PARAMETER           ( KWPTYP = 'PLATE_TYPE'           )

      CHARACTER*(*)         KWVSCL
      PARAMETER           ( KWVSCL = 'FINE_VOXEL_SCALE'     )

      CHARACTER*(*)         KWCSCL
      PARAMETER           ( KWCSCL = 'COARSE_VOXEL_SCALE'   )

      CHARACTER*(*)         KWWRAP
      PARAMETER           ( KWWRAP = 'WRAP_LONGITUDE' )

      CHARACTER*(*)         KWNCAP
      PARAMETER           ( KWNCAP = 'MAKE_NORTH_POLAR_CAP' )

      CHARACTER*(*)         KWSCAP
      PARAMETER           ( KWSCAP = 'MAKE_SOUTH_POLAR_CAP' )

      CHARACTER*(*)         KWRMAJ
      PARAMETER           ( KWRMAJ = 'INPUT_GRID_ORDER_ROW_MAJOR' )

      CHARACTER*(*)         KWTOPD
      PARAMETER           ( KWTOPD = 'COLUMN_VALUE_ORDER_TOP_DOWN' )

      CHARACTER*(*)         KWLFTR
      PARAMETER           ( KWLFTR = 'ROW_VALUE_ORDER_LEFT_RIGHT' )

      CHARACTER*(*)         KWHREF
      PARAMETER           ( KWHREF = 'HEIGHT_REFERENCE' )

      CHARACTER*(*)         KWHSCL
      PARAMETER           ( KWHSCL = 'HEIGHT_SCALE' )

      CHARACTER*(*)         KWNROW
      PARAMETER           ( KWNROW = 'NUMBER_OF_ROWS' )

      CHARACTER*(*)         KWNCOL
      PARAMETER           ( KWNCOL = 'NUMBER_OF_COLUMNS' )

      CHARACTER*(*)         KWLLON
      PARAMETER           ( KWLLON = 'LEFT_COLUMN_LONGITUDE' )

      CHARACTER*(*)         KWLX
      PARAMETER           ( KWLX   = 'LEFT_COLUMN_X_COORDINATE' )

      CHARACTER*(*)         KWTLAT
      PARAMETER           ( KWTLAT = 'TOP_ROW_LATITUDE' )

      CHARACTER*(*)         KWTY
      PARAMETER           ( KWTY   = 'TOP_ROW_Y_COORDINATE' )

      CHARACTER*(*)         KWCSTP
      PARAMETER           ( KWCSTP = 'COLUMN_STEP_SIZE' )

      CHARACTER*(*)         KWRSTP
      PARAMETER           ( KWRSTP = 'ROW_STEP_SIZE' )

      CHARACTER*(*)         KWVTXM
      PARAMETER           ( KWVTXM = 'MAKE_VERTEX_PLATE_MAP' )

      
      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )

      INTEGER               SHORT
      PARAMETER           ( SHORT  = 15 )

      INTEGER               TIMLEN
      PARAMETER           ( TIMLEN = 50 )

      INTEGER               NAMLEN
      PARAMETER           ( NAMLEN = 36 )

      INTEGER               NSYS
      PARAMETER           ( NSYS   = 4 )

C
C     Upper bound on coarse voxel scale. This is used to weed out
C     nonsense input values.
C
      INTEGER               MAXCGS
      PARAMETER           ( MAXCGS = 100 )

C
C     Local variables
C
      CHARACTER*(NAMLEN)    SURFNM
      CHARACTER*(LNSIZE)    CVAL
      CHARACTER*(NAMLEN)    CENTNM
      CHARACTER*(LNSIZE)    CSYNMS ( NSYS )
      CHARACTER*(LNSIZE)    SYSNAM
      CHARACTER*(1)         VTYPE
      CHARACTER*(FILSIZ)    LSK
      CHARACTER*(TIMLEN)    TIMSTR
      CHARACTER*(LNSIZE)    UNISTR ( 2 )
      CHARACTER*(LNSIZE)    WORDS  ( 2 )

      DOUBLE PRECISION      F
      DOUBLE PRECISION      RE
      DOUBLE PRECISION      RP

      INTEGER               FRAMID
      INTEGER               I
      INTEGER               N
      INTEGER               NTOKEN

      LOGICAL               CSFND
      LOGICAL               FSFND
      LOGICAL               FOUND
      LOGICAL               INIT


C
C     Saved variables
C     
      SAVE                 INIT
      SAVE                 CSYNMS

C
C     Initial values
C
      DATA                 CSYNMS / 'LATITUDINAL',
     .                              'CYLINDRICAL',
     .                              'RECTANGULAR',
     .                              'PLANETODETIC' /

      DATA                 INIT   / .FALSE.       /


      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'PRCSET' )

      CALL SIGERR ( 'SPICE(BOGUSENTRY)' )

      CALL CHKOUT ( 'PRCSET' )
      RETURN

       





C$Procedure GETSET ( Get setup file information for MKDSK )
 
      ENTRY GETSET ( SETUP, INPUT, OUTPUT, CMTFIL )
 
C$ Abstract
C
C     Get the names of the input shape file, the output DSK file, and
C     the comment file from an MKDSK setup file.
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
C
C     CHARACTER*(*)         SETUP
C     CHARACTER*(*)         INPUT
C     CHARACTER*(*)         OUTPUT
C     CHARACTER*(*)         CMTFIL
C
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     SETUP      I   Name of setup file.
C     INPUT     I/O  Name of input shape file.
C     OUTPUT    I/O  Name of output DSK file.
C     CMTFIL     O   Name of comment file.
C
C$ Detailed_Input
C
C     SETUP          is the name of an MKDSK setup file.
C
C     INPUT          is the name of a shape file to be converted
C                    to DSK format. This file conforms to the 
C                    format specification given by the MKDSK
C                    User's Guide.
C
C     OUTPUT         is the name of a DSK file to be written.
C                    OUTPUT must be a new file.
C
C$ Detailed_Output
C
C     CMTFIL         is the name of a comment file whose contents
C                    are to be added to the comment area of
C                    the DSK file created by this program.  The
C                    comment file contents precede the default
C                    comments added by MKDSK.
C
C                    If no comment file is specified, CMTFIL is
C                    returned blank.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the setup file name is blank, the error 
C         SPICE(BLANKFILENAME) is signaled.
C
C     2)  If the setup file doesn't exist, the error 
C         SPICE(FILENOTFOUND) is signaled.
C
C     3)  If the name of the surface to be represented by the DSK is
C         not specified in the setup file, the error
C         SPICE(NOSURFACENAME) is signaled.
C
C     4)  If the surface name is present but cannot be mapped to an
C         integer ID code, the error SPICE(NOTRANSLATION) is signaled.
C
C     5)  If the input shape file doesn't exist, the error
C         SPICE(FILENOTFOUND) is signaled.
C
C     6)  If the DSK file name is not specified on the command line
C         and doesn't appear in the setup file, the error
C         SPICE(NOFILESPEC) is signaled.
C
C     7)  If the output file name matches that of an existing file,
C         the error SPICE(FILEEXISTS) is signaled.
C
C     8)  If a comment file keyword is present in the setup file
C         but the associated value does not parse as a quoted string,
C         the error SPICE(TYPEMISMATCH) is signaled.
C
C     9)  If an DSK start time is present in the setup file
C         but the associated value does not parse as a quoted
C         string, the error SPICE(TYPEMISMATCH) is signaled.
C
C     10) If an DSK stop time is present in the setup file
C         but the associated value does not parse as a quoted
C         string, the error SPICE(TYPEMISMATCH) is signaled.
C
C
C$ Files
C
C     See the descriptions of INPUT and OUTPUT above.
C
C$ Particulars
C
C     None.
C
C$ Examples
C
C     None.
C
C$ Restrictions
C
C     This routine is intended for use only within the program 
C     MKDSK.
C
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman    (JPL)
C     B.V. Semenov    (JPL)
C
C$ Version
C
C-    MKDSK Version 3.0.0, 24-FEB-2017 (NJB) (BVS)
C
C        Leapseconds kernel assignment is now optional.
C     
C
C        19-JAN-2016 (NJB)
C
C        Updated to support surface name-ID translation
C        and new coordinate systems.
C
C        Alpha DSK Version 1.0.0, 15-APR-2010 (NJB)
C
C-&

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'GETSET' ) 

C
C     Check the setup file name.
C
      IF ( SETUP .EQ. ' ' ) THEN

         CALL SETMSG ( 'Setup file name may not be blank.' )
         CALL SIGERR ( 'SPICE(BLANKFILENAME)'              )
         CALL CHKOUT ( 'GETSET'                            )
         RETURN

      END IF

      IF ( .NOT. EXISTS(SETUP) ) THEN

         CALL SETMSG ( 'Setup file <#> was not found.' )
         CALL ERRCH  ( '#', SETUP                      )
         CALL SIGERR ( 'SPICE(FILENOTFOUND)'           )
         CALL CHKOUT ( 'GETSET'                        )
         RETURN
         
      END IF

C
C     Load the setup file.
C      
      CALL FURNSH ( SETUP )

      IF ( FAILED() ) THEN
      
         CALL CHKOUT ( 'GETSET' )
         RETURN

      END IF
      

C
C     Check the input file name.  If the name is available, it 
C     supersedes an input file name supplied in the setup file.
C     If the name is blank, the setup file must supply the name.
C
      IF ( INPUT .EQ. ' ' ) THEN         
C
C        Extract the input file name from the kernel pool.
C        Set the INPUT argument to the specified file name.
C
         CALL GCPOOL ( KWINP, 1, 1, N, INPUT, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'Input file was not specified on the ' //
     .                    'command line or in the setup file.'    )
            CALL SIGERR ( 'SPICE(NOFILESPEC)'                     )
            CALL CHKOUT ( 'GETSET'                                )
            RETURN

         END IF

      END IF

      IF ( .NOT. EXISTS(INPUT) ) THEN

         CALL SETMSG ( 'Input file <#> was not found.' )
         CALL ERRCH  ( '#', INPUT                      )
         CALL SIGERR ( 'SPICE(FILENOTFOUND)'           )
         CALL CHKOUT ( 'GETSET'                        )
         RETURN

      END IF
      
C
C     Check the output file name.  If the name is available, it 
C     supersedes an output file name supplied in the setup file.
C     If the name is blank, the setup file must supply the name.
C
      IF ( OUTPUT .EQ. ' ' ) THEN
C
C        Extract the output file name from the kernel pool.
C        Set the INPUT argument to the specified file name.
C
         CALL GCPOOL ( KWOUT, 1, 1, N, OUTPUT, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'Output file was not specified on the ' //
     .                    'command line or in the setup file.'     )
            CALL SIGERR ( 'SPICE(NOFILESPEC)'                      )
            CALL CHKOUT ( 'GETSET'                                 )
            RETURN

         END IF

      END IF

      IF ( EXISTS(OUTPUT)  ) THEN

         CALL SETMSG ( 'Output file <#> already exists.' )
         CALL ERRCH  ( '#', OUTPUT                       )
         CALL SIGERR ( 'SPICE(FILEEXISTS)'               )
         CALL CHKOUT ( 'GETSET'                          )
         RETURN

      END IF


C
C     Obtain the name of the leapseconds kernel and load the kernel.
C
      CALL GCPOOL ( KWLSK, 1, 1, N, LSK, FOUND )

      IF ( FOUND ) THEN

         CALL FURNSH ( LSK )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'GETSET' )
            RETURN
         END IF 
      
      END IF

C
C     See whether a comment file specification was given.
C     Capture the file name if so.
C
      CALL DTPOOL ( KWCMNT, FOUND, N, VTYPE )

      IF ( FOUND ) THEN

         IF ( VTYPE .NE. 'C' ) THEN

            CALL SETMSG ( 'Comment file name was not given a '       //
     .                    'character string value in the setup file.' )
            CALL SIGERR ( 'SPICE(TYPEMISMATCH)'                       )
            CALL CHKOUT ( 'GETSET'                                    )
            RETURN

         END IF

         CALL GCPOOL ( KWCMNT, 1, 1, N, CMTFIL, FOUND ) 
         
      END IF

      IF ( .NOT. FOUND ) THEN
         CMTFIL =  ' '
      END IF

      INIT = .TRUE.

      CALL CHKOUT ( 'GETSET' )
      RETURN






C
C     Get segment data type.
C
      ENTRY GETTYP ( DTYPE )

      CALL CHKIN( 'GETTYP' )

      CALL GIPOOL ( KWDTYP, 1, 1, N, DTYPE, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No segment data type was provided in the '   
     .   //            'setup file.'                              )
         CALL SIGERR ( 'SPICE(MISSINGDATATYPE)'                   )
         CALL CHKOUT ( 'GETTYP'                                   )
         RETURN

      END IF


      CALL CHKOUT( 'GETTYP' )
      RETURN

     


C
C     Get general DSK parameters.
C
      ENTRY GETGEN ( SURFID,  CENTID,  FRAME,   FIRST,
     .               LAST,    DCLASS,  DTYPE,   AUNITS,
     .               DUNITS,  CORSYS,  CORPAR,  MNCOR1,   
     .               MXCOR1,  MNCOR2,  MXCOR2,  MAKVPM )

C
C     Fetch the surface name; convert the name to an ID code.
C
      CALL CHKIN( 'GETGEN' )


      CALL GCPOOL ( KWSURF, 1, 1, N, SURFNM, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No surface name was provided in the '   
     .   //            'setup file. Note that the surface must '
     .   //            'be specified by a string.'               )
         CALL SIGERR ( 'SPICE(MISSINGSURFACE)'                   )
         CALL CHKOUT ( 'GETGEN'                                  )
         RETURN

      END IF

C
C     Fetch the central body name; convert the name to an ID code.
C
      CALL GCPOOL ( KWCENT, 1, 1, N, CENTNM, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No central body name was provided in the '   
     .   //            'setup file. Note that the body must '
     .   //            'be specified by a string.'               )
         CALL SIGERR ( 'SPICE(MISSINGCENTER)'                    )
         CALL CHKOUT ( 'GETGEN'                                  )
         RETURN

      END IF

C
C     Convert the central body name to an ID code.
C
      CALL BODS2C ( CENTNM, CENTID, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'The central body name # could not be mapped '   
     .   //            'to an integer ID code.'                      )
         CALL ERRCH  ( '#', CENTNM                                   )
         CALL SIGERR ( 'SPICE(NOTRANSLATION)'                        )
         CALL CHKOUT ( 'GETGEN'                                      )
         RETURN

      END IF

C
C     Convert the surface name to an ID code.
C
      CALL SRFSCC ( SURFNM, CENTID, SURFID, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'The surface name # could not be mapped '   
     .   //            'to an integer ID code.'                 )
         CALL ERRCH  ( '#', SURFNM                              )
         CALL SIGERR ( 'SPICE(NOTRANSLATION)'                   )
         CALL CHKOUT ( 'GETGEN'                                 )
         RETURN

      END IF

C
C     See whether an output DSK start time was given.
C     Capture the value if so.  
C
      CALL DTPOOL ( KWSTRT, FOUND, N, VTYPE )

      IF ( FOUND ) THEN

         IF ( VTYPE .NE. 'C' ) THEN

            CALL SETMSG ( 'DSK start time was not '                  //
     .                    'given a character string value in the '   //
     .                    'setup file.'                               )
            CALL SIGERR ( 'SPICE(TYPEMISMATCH)'                       )
            CALL CHKOUT ( 'GETGEN'                                    )
            RETURN

         END IF

         CALL GCPOOL ( KWSTRT, 1, 1, N, TIMSTR, FOUND ) 
         
         IF ( FOUND ) THEN

            CALL STR2ET ( TIMSTR, FIRST )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'GETGEN' )
               RETURN
            END IF

         ELSE

            CALL SETMSG ( 'DTPOOL says a start time was provided '
     .      //            'in the setup file, but GCPOOL '
     .      //            'can''t find it (?)'                   )
            CALL SIGERR ( 'SPICE(BUG)'                           )
            CALL CHKOUT ( 'GETGEN'                               )
            RETURN

         END IF

      ELSE

         CALL SETMSG ( 'No start time was provided in the '   
     .   //            'setup file.'                        )
         CALL SIGERR ( 'SPICE(NOSTARTTIME)'                 )
         CALL CHKOUT ( 'GETGEN'                             )
         RETURN

      END IF

C
C     See whether an output DSK stop time was given.
C     Capture the value if so. 
C
      CALL DTPOOL ( KWSTOP, FOUND, N, VTYPE )

      IF ( FOUND ) THEN

         IF ( VTYPE .NE. 'C' ) THEN

            CALL SETMSG ( 'DSK stop time was not '                   //
     .                    'given a character string value in the '   //
     .                    'setup file.'                               )
            CALL SIGERR ( 'SPICE(TYPEMISMATCH)'                       )
            CALL CHKOUT ( 'GETGEN'                                    )
            RETURN

         END IF

         CALL GCPOOL ( KWSTOP, 1, 1, N, TIMSTR, FOUND ) 
         
         IF ( FOUND ) THEN

            CALL STR2ET ( TIMSTR, LAST )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( 'GETGEN' )
               RETURN
            END IF

         ELSE

            CALL SETMSG ( 'DTPOOL says a stop time was provided '
     .      //            'in the setup file, but GCPOOL '
     .      //            'can''t find it (?)'                 )
            CALL SIGERR ( 'SPICE(BUG)'                         )
            CALL CHKOUT ( 'GETGEN'                             )
            RETURN

         END IF

      ELSE

         CALL SETMSG ( 'No stop time was provided in the '   
     .   //            'setup file.'                        )
         CALL SIGERR ( 'SPICE(NOSTOPTIME)'                  )
         CALL CHKOUT ( 'GETGEN'                             )
         RETURN

      END IF 


C
C     Get segment data class.
C
      CALL GIPOOL ( KWCLSS, 1, 1, N, DCLASS, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No segment data class was provided in the '   
     .   //            'setup file.'                              )
         CALL SIGERR ( 'SPICE(MISSINGDATACLASS)'                  )
         CALL CHKOUT ( 'GETGEN'                                   )
         RETURN

      END IF

C
C     Get segment data type.
C
      CALL GIPOOL ( KWDTYP, 1, 1, N, DTYPE, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No segment data type was provided in the '   
     .   //            'setup file.'                              )
         CALL SIGERR ( 'SPICE(MISSINGDATATYPE)'                   )
         CALL CHKOUT ( 'GETGEN'                                   )
         RETURN

      END IF

C
C     See whether a units string was provided.
C
      CALL DTPOOL ( KWUNIT, FOUND, N, VTYPE )

      IF ( FOUND ) THEN

         IF ( VTYPE .NE. 'C' ) THEN

            CALL SETMSG ( 'Units specification was not '             //
     .                    'given a character string value in the '   //
     .                    'setup file.'                               )
            CALL SIGERR ( 'SPICE(TYPEMISMATCH)'                       )
            CALL CHKOUT ( 'GETGEN'                                    )
            RETURN

         END IF

         CALL GCPOOL ( KWUNIT, 1, 2, N, UNISTR, FOUND ) 
         
         IF ( FOUND ) THEN

            IF ( N .NE. 2 ) THEN
C
C              We need both distance and angular units.
C
               CALL SETMSG ( 'Improperly formatted unit '
     .         //            'specification in setup file: number '
     .         //            'of strings found was #. Both angular '
     .         //            'units and distance units must be '
     .         //            'specified.'                          )  
               CALL ERRINT ( '#', N                                )
               CALL SIGERR ( 'SPICE(SYNTAXERROR)'                  )
               CALL CHKOUT ( 'GETGEN'                              )
               RETURN

            END IF


C
C           Parse the unit specifications.
C
            DO I = 1, N

               CALL LPARSE ( UNISTR(I), '=', 2, NTOKEN, WORDS )

               IF ( NTOKEN .NE. 2 ) THEN

                  CALL SETMSG ( 'Improperly formatted unit '
     .            //            'specification in setup file: #'    )
                  CALL ERRCH  ( '#', UNISTR(I)                      )
                  CALL SIGERR ( 'SPICE(SYNTAXERROR)'                )
                  CALL CHKOUT ( 'GETGEN'                            )
                  RETURN

               END IF

               IF ( EQSTR( WORDS(1), 'ANGLES' )  ) THEN

                  CALL LJUST ( WORDS(2), AUNITS )
                  CALL UCASE ( AUNITS,   AUNITS )
                  
               ELSE IF ( EQSTR( WORDS(1), 'DISTANCES' )  ) THEN

                  CALL LJUST ( WORDS(2), DUNITS )
                  CALL UCASE ( DUNITS,   DUNITS )
C
C                 Map "KILOMETERS" to "KM"; the latter is 
C                 recognized by CONVRT.
C
                  IF ( DUNITS .EQ. 'KILOMETERS' ) THEN
                     DUNITS = 'KM'
                  END IF

               ELSE

                  CALL SETMSG ( 'Unrecognized dimension # in unit '
     .            //            'specification in setup file: #'    )
                  CALL ERRCH  ( '#', WORDS(1)                       )
                  CALL ERRCH  ( '#', UNISTR(I)                      )
                  CALL SIGERR ( 'SPICE(SYNTAXERROR)'                )
                  CALL CHKOUT ( 'GETGEN'                            )
                  RETURN

               END IF

            END DO

         ELSE

            CALL SETMSG ( 'DTPOOL says a units assignment was '
     .      //            'provided in the setup file, '
     .      //            'but GCPOOL can''t find it (?)'        )
            CALL SIGERR ( 'SPICE(BUG)'                           )
            CALL CHKOUT ( 'GETGEN'                               )
            RETURN

         END IF

      ELSE

         CALL SETMSG ( 'No unit specification was provided in the '
     .   //            'setup file.'                               )
         CALL SIGERR ( 'SPICE(NOUNITSPEC)'                         )
         CALL CHKOUT ( 'GETGEN'                                    )
         RETURN


      END IF 


C
C     Get the reference frame name.     
C
      CALL GCPOOL ( KWFRAM, 1, 1, N, FRAME, FOUND )
 
      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No reference frame name was provided in the '   
     .   //            'setup file.'                                  )
         CALL SIGERR ( 'SPICE(MISSINGFRAME)'                          )
         CALL CHKOUT ( 'GETGEN'                                       )
         RETURN

      END IF

C
C     Verify that the frame can be mapped to an ID code.
C
      CALL NAMFRM ( FRAME, FRAMID )

      IF ( FRAMID .EQ. 0 ) THEN

         CALL SETMSG ( 'Reference frame name # could not be mapped '   
     .   //            'to a frame ID code.'                        )
         CALL ERRCH  ( '#',  FRAME                                  )
         CALL SIGERR ( 'SPICE(NOTRANSLATION)'                       )
         CALL CHKOUT ( 'GETGEN'                                     )
         RETURN
         
      END IF

C
C     Get the coordinate system name.     
C
      CALL GCPOOL ( KWCSYS, 1, 1, N, SYSNAM, FOUND )
 
      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No coordinate system name was provided in the '
     .   //            'setup file.'                                   )
         CALL SIGERR ( 'SPICE(MISSINGCOORDSYS)'                        )
         CALL CHKOUT ( 'GETGEN'                                        )
         RETURN

      END IF

C
C     Map the coordinate system name to an ID.
C
      CORSYS = ESRCHC( SYSNAM, NSYS, CSYNMS )

      IF ( CORSYS .EQ. 0 ) THEN

         CALL SETMSG ( 'Coordinate system name # was not recognized. ')
         CALL ERRCH  ( '#', SYSNAM                                    )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                          )
         CALL CHKOUT ( 'GETGEN'                                       )
         RETURN

      ELSE IF ( CORSYS .EQ. CYLSYS ) THEN

         CALL SETMSG ( 'Cylindrical coordinates are not supported '
     .   //            'by this program.'                             )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                          )
         CALL CHKOUT ( 'GETGEN'                                       )
         RETURN

      END IF

C
C     Get coordinate bounds.
C
      IF ( ( CORSYS .EQ. LATSYS ) .OR. ( CORSYS .EQ. PDTSYS )  ) THEN

         CALL GDPOOL ( KWMNLO, 1, 1, N, MNCOR1, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No minimum longitude was provided in the '   
     .      //            'setup file.'                              )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                 )
            CALL CHKOUT ( 'GETGEN'                                   )
            RETURN

         END IF

         CALL GDPOOL ( KWMXLO, 1, 1, N, MXCOR1, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No maximum longitude was provided in the '   
     .      //            'setup file.'                              )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                 )
            CALL CHKOUT ( 'GETGEN'                                   )
            RETURN

         END IF

         CALL GDPOOL ( KWMNLA, 1, 1, N, MNCOR2, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No minimum latitude was provided in the '
     .      //            'setup file.'                              )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                 )
            CALL CHKOUT ( 'GETGEN'                                   )
            RETURN

         END IF

         CALL GDPOOL ( KWMXLA, 1, 1, N, MXCOR2, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No maximum latitude was provided in the '
     .      //            'setup file.'                              )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                 )
            CALL CHKOUT ( 'GETGEN'                                   )
            RETURN

         END IF

      END IF


      IF ( CORSYS .EQ. RECSYS  ) THEN

         CALL GDPOOL ( KWMNX, 1, 1, N, MNCOR1, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No minimum X-coordinate was provided in the '
     .      //            'setup file.'                                )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                   )
            CALL CHKOUT ( 'GETGEN'                                     )
            RETURN

         END IF

         CALL GDPOOL ( KWMXX, 1, 1, N, MXCOR1, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No maximum X-coordinate was provided in the '
     .      //            'setup file.'                                )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                   )
            CALL CHKOUT ( 'GETGEN'                                     )
            RETURN

         END IF

         CALL GDPOOL ( KWMNY, 1, 1, N, MNCOR2, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No minimum Y-coordinate was provided in the '
     .      //            'setup file.'                                )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                   )
            CALL CHKOUT ( 'GETGEN'                                     )
            RETURN

         END IF

         CALL GDPOOL ( KWMXY, 1, 1, N, MXCOR2, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No maximum Y-coordinate was provided in the '
     .      //            'setup file.'                                )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                   )
            CALL CHKOUT ( 'GETGEN'                                     )
            RETURN

         END IF

      END IF


C
C     Get coordinate parameters, if necessary.
C
      CALL CLEARD ( NSYPAR, CORPAR )

      IF ( CORSYS .EQ. PDTSYS ) THEN

         CALL GDPOOL ( KWRE, 1, 1, N, RE, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No equatorial radius for the planetodetic '
     .      //            'coordinate system''s reference spheroid '
     .      //            'was provided in the setup file.'          )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                 )
            CALL CHKOUT ( 'GETGEN'                                   )
            RETURN

         END IF

         CALL GDPOOL ( KWRP, 1, 1, N, RP, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No polar radius for the planetodetic '
     .      //            'coordinate system''s reference spheroid '
     .      //            'was provided in the setup file.'          )
            CALL SIGERR ( 'SPICE(MISSINGCOORDBOUND)'                 )
            CALL CHKOUT ( 'GETGEN'                                   )
            RETURN

         END IF

         IF ( ( RE .LE. 0.D0 ) .OR. ( RP .LE. 0.D0 ) ) THEN

            CALL SETMSG ( 'In the setup file, the equatorial '
     .      //            'radius = #; the polar radius = #. '
     .      //            'Both radii must be positive.'      )
            CALL ERRDP  ( '#', RE                             )
            CALL ERRDP  ( '#', RP                             )
            CALL SIGERR ( 'SPICE(INVALIDRADII)'               )
            CALL CHKOUT ( 'GETGEN'                            )
            RETURN

         END IF

C
C        Map the radii to a flattening coefficient.
C
         CORPAR( 1 ) = RE
         F                = ( RE - RP ) / RE
         CORPAR( 2 ) = F

      END IF

C
C     Fetch the "make vertex-plate mapping" flag, if it's present.
C
      CALL GCPOOL ( KWVTXM, 1, 1, N, CVAL, FOUND )

      IF ( .NOT. FOUND ) THEN
C
C        The flag is not required to be present. By default, 
C        no map is created.
C
         MAKVPM = .FALSE.

      ELSE

         IF (  EQSTR( CVAL, 'YES' )  ) THEN

            MAKVPM = .TRUE.

         ELSE IF (  EQSTR( CVAL, 'NO' )  ) THEN

            MAKVPM = .FALSE.

         ELSE

            CALL SETMSG ( '"Make vertex-plate map" flag must '
     .      //            'be YES or NO but was #.'           )
            CALL ERRCH  ( '#', CVAL                           )
            CALL SIGERR ( 'SPICE(INVALIDFLAG)'                )
            CALL CHKOUT ( 'GETGEN'                            )
            RETURN

         END IF

      END IF



      CALL CHKOUT( 'GETGEN' )
      RETURN



C
C     Get type 2 parameters.
C
      ENTRY GETP02 ( PLTYPE, VOXSCL, CGRSCL )

      CALL CHKIN ( 'GETP02' )

C
C     Get the plate model type.     
C
      CALL GIPOOL ( KWPTYP, 1, 1, N, PLTYPE, FOUND )
 
      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No plate model type was provided in the '   
     .   //            'setup file.'                                   )
         CALL SIGERR ( 'SPICE(MISSINGPLATETYPE)'                       )
         CALL CHKOUT ( 'GETP02'                                        )
         RETURN

      END IF

C
C     Get the voxel scale.    
C
      CALL GDPOOL ( KWVSCL, 1, 1, N, VOXSCL, FSFND )

C
C     Get the coarse voxel scale.    
C
      CALL GIPOOL ( KWCSCL, 1, 1, N, CGRSCL, CSFND )

C
C     It's ok if no scales were provided; otherwise both
C     must be provided.
C 
      IF (  ( .NOT. FSFND ) .AND. ( .NOT. CSFND )  ) THEN
C
C        Return scales set to zero. The scales will be
C        determined automatically.
C
         VOXSCL = 0.D0
         CGRSCL = 0

      ELSE IF ( FSFND .AND. CSFND ) THEN
C
C        Both scales were provided; check them.
C
         IF ( VOXSCL .LE. 0.D0 ) THEN

            CALL SETMSG ( 'Fine voxel scale must be strictly '
     .      //            'positive but was #. (This scale '
     .      //            'normally should be greater than '
     .      //            'or equal to 1.0.)'                 )
            CALL ERRDP  ( '#', VOXSCL                         )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'            )
            CALL CHKOUT ( 'GETP02'                            )
            RETURN

         END IF

         IF ( CGRSCL .LT. 1 ) THEN

            CALL SETMSG ( 'Coarse voxel scale must be greater ' 
     .      //            'than or equal to 1, but was #.'     )
            CALL ERRINT ( '#', CGRSCL                          )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'             )
            CALL CHKOUT ( 'GETP02'                             )
            RETURN

         ELSE IF ( CGRSCL .GT. MAXCGS ) THEN

            CALL SETMSG ( 'Coarse voxel scale must be less ' 
     .      //            'than or equal to #, but was #. '
     .      //            '(Normally this scale should not '
     .      //            'exceed 20.)'                        )
            CALL ERRINT ( '#', MAXCGS                          )
            CALL ERRINT ( '#', CGRSCL                          )
            CALL SIGERR ( 'SPICE(VALUEOUTOFRANGE)'             )
            CALL CHKOUT ( 'GETP02'                             )
            RETURN

         END IF



      ELSE IF ( .NOT. FSFND ) THEN

         CALL SETMSG ( 'No fine voxel scale was provided in the '   
     .   //            'setup file, but a coarse voxel scale was '
     .   //            'provided. Either add an assignment for '
     .   //            'the fine voxel scale, or provide neither '
     .   //            'scale, in which case the scales will '
     .   //            'be set automatically.'                   )
         CALL SIGERR ( 'SPICE(MISSINGVOXELSCALE)'                )
         CALL CHKOUT ( 'GETP02'                                  )
         RETURN

      ELSE IF ( .NOT. CSFND ) THEN

         CALL SETMSG ( 'No coarse voxel scale was provided in the '   
     .   //            'setup file, but a fine voxel scale was '
     .   //            'provided. Either add an assignment for '
     .   //            'the coarse voxel scale, or provide neither '
     .   //            'scale, in which case the scales will '
     .   //            'be set automatically.'                   )
         CALL SIGERR ( 'SPICE(MISSINGVOXELSCALE)'                )
         CALL CHKOUT ( 'GETP02'                                  )
         RETURN

      END IF

      CALL CHKOUT( 'GETP02' )
      RETURN



C
C     Get plate type 5 (height grid) parameters.
C
      ENTRY GETG05 ( CORSYS, WRAP,   MKNCAP, MKSCAP, ROWMAJ,
     .               TOPDWN, LEFTRT, REFVAL, HSCALE, NCOLS, 
     .               NROWS,  LFTCOR, TOPCOR, COLSTP, ROWSTP )

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'GETG05' )

C
C     Initialize flags that are set conditionally.
C
      WRAP   = .FALSE.
      MKNCAP = .FALSE.
      MKSCAP = .FALSE.

C
C     Initialize REFVAL.
C
      REFVAL = 0.D0


      IF ( CORSYS .EQ. RECSYS ) THEN     

         CALL GCPOOL ( KWWRAP, 1, 1, N, CVAL, FOUND )
 
         IF ( FOUND ) THEN

            IF (  EQSTR( CVAL, 'YES' )  ) THEN

               CALL SETMSG ( 'The longitude wrap flag does not apply '
     .         //            'to rectangular coordinates. Set the '
     .         //            'flag value to ''NO'' or delete the '
     .         //            'longitude wrap assignment from '
     .         //            'the setup file.'                         )
               CALL SIGERR ( 'SPICE(SPURIOUSKEYWORD)'                  )
               CALL CHKOUT ( 'GETG05'                                  )
               RETURN

            END IF

         END IF

         CALL GCPOOL ( KWNCAP, 1, 1, N, CVAL, FOUND )
 
         IF ( FOUND ) THEN

            IF (  EQSTR( CVAL, 'YES' )  ) THEN

               CALL SETMSG ( 'The north polar cap flag does not apply '
     .         //            'to rectangular coordinates. Set the '
     .         //            'flag value to ''NO'' or delete the '
     .         //            'north polar cap flag assignment from '
     .         //            'the setup file.'                         )
               CALL SIGERR ( 'SPICE(SPURIOUSKEYWORD)'                  )
               CALL CHKOUT ( 'GETG05'                                  )
               RETURN

            END IF

         END IF

         CALL GCPOOL ( KWSCAP, 1, 1, N, CVAL, FOUND )
 
         IF ( FOUND ) THEN

            IF (  EQSTR( CVAL, 'YES' )  ) THEN
 
               CALL SETMSG ( 'The south polar cap flag does not apply '
     .         //            'to rectangular coordinates.  Set the '
     .         //            'flag value to ''NO'' or delete the '
     .         //            'south polar cap flag assignment from '
     .         //            'the setup file.'                         )
               CALL SIGERR ( 'SPICE(SPURIOUSKEYWORD)'                  )
               CALL CHKOUT ( 'GETG05'                                  )
               RETURN

            END IF

         END IF

C
C        Get the coordinate of the top row.
C
         CALL GDPOOL ( KWTY, 1, 1, N, TOPCOR, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No Y-coordinate of the top row was  '   
     .      //            'provided in the setup file.'          )
            CALL SIGERR ( 'SPICE(MISSINGTOPCOR)'                 )
            CALL CHKOUT ( 'GETG05'                               )
            RETURN

        END IF

C
C        Get the coordinate of the left column.
C
         CALL GDPOOL ( KWLX, 1, 1, N, LFTCOR, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No X-coordinate of the left column was  '   
     .      //            'provided in the setup file.'            )
            CALL SIGERR ( 'SPICE(MISSINGLEFTCOR)'                  )
            CALL CHKOUT ( 'GETG05'                                 )
            RETURN

         END IF

      ELSE
C
C        This is a lon/lat coordinate system.
C
C        Get the longitude wrap flag.
C
         CALL GCPOOL ( KWWRAP, 1, 1, N, CVAL, FOUND )
 
         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No longitude wrap flag was provided in the '
     .      //            'setup file.'                                )
            CALL SIGERR ( 'SPICE(MISSINGWRAPFLAG)'                     )
            CALL CHKOUT ( 'GETG05'                                     )
            RETURN

         END IF

         IF (  EQSTR( CVAL, 'YES' )  ) THEN

            WRAP = .TRUE.

         ELSE IF (  EQSTR( CVAL, 'NO' )  ) THEN

            WRAP = .FALSE.

         ELSE

            CALL SETMSG ( 'Longitude wrap flag must be YES or NO but '
     .      //            'was #.'                                     )
            CALL ERRCH  ( '#', CVAL                                    )
            CALL SIGERR ( 'SPICE(INVALIDFLAG)'                         )
            CALL CHKOUT ( 'GETG05'                                     )
            RETURN

         END IF


C
C        Get the north polar cap flag.
C
         CALL GCPOOL ( KWNCAP, 1, 1, N, CVAL, FOUND )
 
         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No north polar cap flag was provided in the '
     .      //            'setup file.'                                )
            CALL SIGERR ( 'SPICE(MISSINGNCAPFLAG)'                     )
            CALL CHKOUT ( 'GETG05'                                     )
            RETURN

         END IF

         IF (  EQSTR( CVAL, 'YES' )  ) THEN

            MKNCAP = .TRUE.

         ELSE IF (  EQSTR( CVAL, 'NO' )  ) THEN

            MKNCAP = .FALSE.

         ELSE

            CALL SETMSG ( 'North polar cap flag must be YES or NO but '
     .      //            'was #.'                                     )
            CALL ERRCH  ( '#', CVAL                                    )
            CALL SIGERR ( 'SPICE(INVALIDFLAG)'                         )
            CALL CHKOUT ( 'GETG05'                                     )
            RETURN

         END IF

C
C        Get the south polar cap flag.
C
         CALL GCPOOL ( KWSCAP, 1, 1, N, CVAL, FOUND )
 
         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No south polar cap flag was provided in the '
     .      //            'setup file.'                                )
            CALL SIGERR ( 'SPICE(MISSINGSCAPFLAG)'                     )
            CALL CHKOUT ( 'GETG05'                                     )
            RETURN

         END IF

         IF (  EQSTR( CVAL, 'YES' )  ) THEN

            MKSCAP = .TRUE.

         ELSE IF (  EQSTR( CVAL, 'NO' )  ) THEN

            MKSCAP = .FALSE.

         ELSE

            CALL SETMSG ( 'South polar cap flag must be YES or NO but '
     .      //            'was #.'                                     )
            CALL ERRCH  ( '#', CVAL                                    )
            CALL SIGERR ( 'SPICE(INVALIDFLAG)'                         )
            CALL CHKOUT ( 'GETG05'                                     )
            RETURN

         END IF

C
C        Get the coordinate of the top row.
C
         CALL GDPOOL ( KWTLAT, 1, 1, N, TOPCOR, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No latitude of the top row was  '   
     .      //            'provided in the setup file.'          )
            CALL SIGERR ( 'SPICE(MISSINGTOPCOR)'                 )
            CALL CHKOUT ( 'GETG05'                               )
            RETURN

         END IF

C
C        Get the coordinate of the left column.
C
         CALL GDPOOL ( KWLLON, 1, 1, N, LFTCOR, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No longitude of the left column was  '   
     .      //            'provided in the setup file.'          )
            CALL SIGERR ( 'SPICE(MISSINGLEFTCOR)'                )
            CALL CHKOUT ( 'GETG05'                               )
            RETURN

         END IF

      END IF

C
C     Get the row major grid order flag.
C
      CALL GCPOOL ( KWRMAJ, 1, 1, N, CVAL, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No row major flag was provided in the '   
     .   //            'setup file.'                             )
         CALL SIGERR ( 'SPICE(MISSINGROWMAJFLAG)'                )
         CALL CHKOUT ( 'GETG05'                                  )
         RETURN

      END IF

      IF (  EQSTR( CVAL, 'YES' )  ) THEN

         ROWMAJ = .TRUE.

      ELSE IF (  EQSTR( CVAL, 'NO' )  ) THEN

         ROWMAJ = .FALSE.

      ELSE

         CALL SETMSG ( 'Row major flag must be YES or NO but '     
     .   //            'was #.'                               )
         CALL ERRCH  ( '#', CVAL                              )
         CALL SIGERR ( 'SPICE(INVALIDFLAG)'                   )
         CALL CHKOUT ( 'GETG05'                               )
         RETURN

      END IF

C
C     Get the top-down grid order flag.
C
      CALL GCPOOL ( KWTOPD, 1, 1, N, CVAL, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No top-down flag was provided in the '   
     .   //            'setup file.'                            )
         CALL SIGERR ( 'SPICE(MISSINGTOPDOWNFLAG)'              )
         CALL CHKOUT ( 'GETG05'                                 )
         RETURN

      END IF

      IF (  EQSTR( CVAL, 'YES' )  ) THEN

         TOPDWN = .TRUE.

      ELSE IF (  EQSTR( CVAL, 'NO' )  ) THEN

         TOPDWN = .FALSE.

      ELSE

         CALL SETMSG ( 'Top-down flag must be YES or NO but '     
     .   //            'was #.'                               )
         CALL ERRCH  ( '#', CVAL                              )
         CALL SIGERR ( 'SPICE(INVALIDFLAG)'                   )
         CALL CHKOUT ( 'GETG05'                               )
         RETURN

      END IF


C
C     Get the left-right grid order flag.
C
      CALL GCPOOL ( KWLFTR, 1, 1, N, CVAL, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No left-right flag was provided in the '   
     .   //            'setup file.'                            )
         CALL SIGERR ( 'SPICE(MISSINGLEFTRTFLAG)'               )
         CALL CHKOUT ( 'GETG05'                                 )
         RETURN

      END IF

      IF (  EQSTR( CVAL, 'YES' )  ) THEN

         LEFTRT = .TRUE.

      ELSE IF (  EQSTR( CVAL, 'NO' )  ) THEN

         LEFTRT = .FALSE.

      ELSE

         CALL SETMSG ( 'Left-right flag must be YES or NO but '     
     .   //            'was #.'                               )
         CALL ERRCH  ( '#', CVAL                              )
         CALL SIGERR ( 'SPICE(INVALIDFLAG)'                   )
         CALL CHKOUT ( 'GETG05'                               )
         RETURN

      END IF


      IF ( CORSYS .NE. PDTSYS ) THEN
C
C        Get the height reference value.
C
         CALL GDPOOL ( KWHREF, 1, 1, N, REFVAL, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'No height reference was provided in the '   
     .      //            'setup file.'                            )
            CALL SIGERR ( 'SPICE(MISSINGHEIGHTREF)'                )
            CALL CHKOUT ( 'GETG05'                                 )
            RETURN

         END IF

      END IF

C
C     Get the height scale value.
C
      CALL GDPOOL ( KWHSCL, 1, 1, N, HSCALE, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No height scale factor was provided '   
     .   //            'in the setup file.'                   )
         CALL SIGERR ( 'SPICE(MISSINGHSCALE)'                 )
         CALL CHKOUT ( 'GETG05'                               )
         RETURN

      END IF

C
C     Get the column count.
C
      CALL GIPOOL ( KWNCOL, 1, 1, N, NCOLS, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No column count was provided '   
     .   //            'in the setup file.'           )
         CALL SIGERR ( 'SPICE(MISSINGNCOLS)'          )
         CALL CHKOUT ( 'GETG05'                       )
         RETURN

      END IF

C
C     Get the row count.
C
      CALL GIPOOL ( KWNROW, 1, 1, N, NROWS, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No row count was provided '   
     .   //            'in the setup file.'           )
         CALL SIGERR ( 'SPICE(MISSINGNROWS)'          )
         CALL CHKOUT ( 'GETG05'                       )
         RETURN

      END IF

C
C     Get the column step size.
C
      CALL GDPOOL ( KWCSTP, 1, 1, N, COLSTP, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No column step size was  '   
     .   //            'provided in the setup file.' )
         CALL SIGERR ( 'SPICE(MISSINGCOLSTEP)'       )
         CALL CHKOUT ( 'GETG05'                      )
         RETURN

      END IF

C
C     Get the row step size.
C
      CALL GDPOOL ( KWRSTP, 1, 1, N, ROWSTP, FOUND )

      IF ( .NOT. FOUND ) THEN

         CALL SETMSG ( 'No row step size was  '   
     .   //            'provided in the setup file.' )
         CALL SIGERR ( 'SPICE(MISSINGROWSTEP)'       )
         CALL CHKOUT ( 'GETG05'                      )
         RETURN

      END IF

      CALL CHKOUT ( 'GETG05' )
      RETURN

      END

