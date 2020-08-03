C$Procedure      OCCULT ( find occultation type at time )

      SUBROUTINE OCCULT ( TARG1,  SHAPE1, FRAME1, 
     .                    TARG2,  SHAPE2, FRAME2, 
     .                    ABCORR, OBSRVR, ET,     OCLTID )

C$ Abstract
C
C     Determines the occultation condition (not occulted, partially,
C     etc.) of one target relative to another target as seen by
C     an observer at a given time.
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
C     SPK
C     TIME
C     KERNEL
C
C$ Keywords
C
C     GEOMETRY
C     OCCULTATION
C     ELLIPSOID
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE 'gf.inc'
      INCLUDE 'occult.inc'
      INCLUDE 'dsk.inc'
      INCLUDE 'zzdsk.inc'

      CHARACTER*(*)         TARG1
      CHARACTER*(*)         SHAPE1
      CHARACTER*(*)         FRAME1
      CHARACTER*(*)         TARG2
      CHARACTER*(*)         SHAPE2
      CHARACTER*(*)         FRAME2
      CHARACTER*(*)         ABCORR
      CHARACTER*(*)         OBSRVR
      DOUBLE PRECISION      ET
      INTEGER               OCLTID


C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     TARG1      I   Name or ID of first target.
C     SHAPE1     I   Type of shape model used for first target.
C     FRAME1     I   Body-fixed, body-centered frame for first body.
C     TARG2      I   Name or ID of second target.
C     SHAPE2     I   Type of shape model used for second target.
C     FRAME2     I   Body-fixed, body-centered frame for second body.
C     ABCORR     I   Aberration correction flag.
C     OBSRVR     I   Name or ID of the observer.
C     ET         I   Time of the observation (seconds past J2000).
C     OCLTID     O   Occultation identification code.
C
C$ Detailed_Input
C
C     TARG1      is the name of the first target body. Both object
C                names and NAIF IDs are accepted. For example, both
C                'Moon' and '301' are accepted.
C
C     SHAPE1     is a string indicating the geometric model used to
C                represent the shape of the first target body. The
C                supported options are:
C
C                   'ELLIPSOID'    
C
C                       Use a triaxial ellipsoid model with radius
C                       values provided via the kernel pool. A kernel
C                       variable having a name of the form
C
C                          'BODYnnn_RADII'
C
C                       where nnn represents the NAIF integer code
C                       associated with the body, must be present in
C                       the kernel pool. This variable must be
C                       associated with three numeric values giving the
C                       lengths of the ellipsoid's X, Y, and Z
C                       semi-axes.
C
C                   'POINT'       
C
C                       Treat the body as a single point. When a point
C                       target is specified, the occultation conditions
C                       can only be total, annular, or none.
C
C                   'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                       Use topographic data provided by DSK files to
C                       model the body's shape. These data must be
C                       provided by loaded DSK files.
C
C                       The surface list specification is optional. The
C                       syntax of the list is
C
C                          <surface 1> [, <surface 2>...]
C
C                       If present, it indicates that data only for the
C                       listed surfaces are to be used; however, data
C                       need not be available for all surfaces in the
C                       list. If absent, loaded DSK data for any surface
C                       associated with the target body are used.
C
C                       The surface list may contain surface names or
C                       surface ID codes. Names containing blanks must
C                       be delimited by double quotes, for example
C
C                          SURFACES = "Mars MEGDR 128 PIXEL/DEG"
C                                         
C                       If multiple surfaces are specified, their names
C                       or IDs must be separated by commas.
C
C                       See the Particulars section below for details
C                       concerning use of DSK data.
C
C                The combinations of the shapes of the target bodies
C                TARG1 and TARG2 must be one of:
C
C                   One ELLIPSOID, one POINT
C                   Two ELLIPSOIDs
C                   One DSK, one POINT
C
C                Case and leading or trailing blanks are not
C                significant in the string SHAPE1.
C
C     FRAME1     is the name of the body-fixed, body-centered reference
C                frame associated with the first target body. Examples
C                of such names are 'IAU_SATURN' (for Saturn) and
C                'ITRF93' (for the Earth).
C
C                If the first target body is modeled as a point, FRAME1
C                should be left blank (Ex: ' ').
C
C                Case and leading or trailing blanks bracketing a
C                non-blank frame name are not significant in the string.
C
C     TARG2      is the name of the second target body. See the
C                description of TARG1 above for more details.
C
C     SHAPE2     is the shape specification for the body designated
C                by TARG2. See the description of SHAPE1 above for
C                details.
C
C     FRAME2     is the name of the body-fixed, body-centered reference
C                frame associated with the second target body. See the
C                description of FRAME1 above for more details.
C
C     ABCORR     indicates the aberration corrections to be applied to
C                the state of each target body to account for one-way
C                light time. Stellar aberration corrections are
C                ignored if specified, since these corrections don't
C                improve the accuracy of the occultation determination.
C
C                See the header of the SPICE routine SPKEZR for a
C                detailed description of the aberration correction
C                options. For convenience, the options supported by
C                this routine are listed below:
C
C                   'NONE'     Apply no correction.   
C
C                   'LT'       "Reception" case: correct for
C                              one-way light time using a Newtonian
C                              formulation.
C
C                   'CN'       "Reception" case: converged
C                              Newtonian light time correction.
C
C                   'XLT'      "Transmission" case: correct for
C                              one-way light time using a Newtonian
C                              formulation.
C
C                   'XCN'      "Transmission" case: converged
C                              Newtonian light time correction.
C
C                Case and blanks are not significant in the string
C                ABCORR.
C
C     OBSRVR     is the name of the body from which the occultation
C                is observed. See the description of TARG1 above for
C                more details.
C
C     ET         is the observation time in seconds past the J2000
C                epoch.
C
C                    
C$ Detailed_Output
C
C     OCLTID     is an integer occultation code indicating the geometric
C                relationship of the three bodies.  
C
C                The meaning of the sign of OCLTID is given below.
C
C                    Code sign          Meaning
C                    ---------          ------------------------------
C                       > 0             The second target is
C                                       partially or fully occulted
C                                       by the first.
C
C                       < 0             The first target is 
C                                       partially of fully
C                                       occulted by the second.
C
C                       = 0             No occultation.
C
C                Possible OCLTID values and meanings are given below.
C                The variable names indicate the type of occultation
C                and which target is in the back. For example, TOTAL1
C                represents a total occultation in which the first
C                target is in the back of (or is occulted by) the
C                second target.
C
C                When the target shapes are DSK and POINT, the only
C                possible occultation conditions are total, annular,
C                or none.
C
C                    Name      Code     Meaning
C                    ------    -----    ------------------------------
C                    TOTAL1     -3      Total occultation of first
C                                       target by second.
C
C                    ANNLR1     -2      Annular occultation of first
C                                       target by second. If the second
C                                       target shape is an ellipsoid,
C                                       it does not block the limb of
C                                       the first.
C
C                    PARTL1     -1      Partial occultation of first
C                                       target by second target.
C
C                    NOOCC       0      No occultation or transit: both
C                                       objects are completely visible
C                                       to the observer.
C
C                    PARTL2      1      Partial occultation of second
C                                       target by first target.
C
C                    ANNLR2      2      Annular occultation of second
C                                       target by first.
C
C                    TOTAL2      3      Total occultation of second
C                                       target by first.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the target or observer body names input by the user are
C         not recognized, the error will be diagnosed by a routine in
C         the call tree of this routine.
C
C     2)  If the input shapes are not accepted, the error will be
C         diagnosed by a routine in the call tree of this routine.
C
C     3)  If both input shapes are points, the error will be
C         diagnosed by a routine in the call tree of this routine.
C
C     4)  If the radii of a target body modeled as an ellipsoid cannot
C         be determined by searching the kernel pool for a kernel
C         variable having a name of the form
C
C            'BODYnnn_RADII' 
C
C         where nnn represents the NAIF integer code associated with
C         the body, the error will be diagnosed by a routine in the
C         call tree of this routine.
C
C     5)  If any of the target or observer bodies (TARG1, TARG2, or
C         OBSRVR) are the same, the error will be diagnosed
C         by a routine in the call tree of this routine.
C         
C     6)  If the loaded kernels provide insufficient data to 
C         compute any required state vector, the deficiency will
C         be diagnosed by a routine in the call tree of this routine.
C
C     7)  If an error occurs while reading an SPK or other kernel,
C         the error will be diagnosed by a routine in the call tree 
C         of this routine.
C
C     8)  Invalid aberration correction specifications will be
C         diagnosed by a routine in the call tree of this routine.
C
C     17) If either SHAPE1 or SHAPE2 specifies that the target surface
C         is represented by DSK data, and no DSK files are loaded for
C         the specified target, the error is signaled by a routine in
C         the call tree of this routine.
C
C     18) If either SHAPE1 or SHAPE2 specifies that the target surface
C         is represented by DSK data, but the shape specification is
C         invalid, the error is signaled by a routine in the call tree
C         of this routine.
C
C$ Files
C
C     Appropriate SPICE kernels must be loaded by the calling program
C     before this routine is called.
C
C     The following data are required:
C
C        - SPK data: the calling application must load ephemeris data
C          for the target, source and observer that cover the time
C          period specified by the window CNFINE. If aberration
C          corrections are used, the states of the target bodies and of
C          the observer relative to the solar system barycenter must be
C          calculable from the available ephemeris data. Typically
C          ephemeris data
C          are made available by loading one or more SPK files via
C          FURNSH.
C
C        - PCK data: bodies modeled as triaxial ellipsoids must have
C          semi-axis lengths provided by variables in the kernel pool.
C          Typically these data are made available by loading a text
C          PCK file via FURNSH.
C
C        - FK data: if either of the reference frames designated by
C          FRAME1 or FRAME2 are not built in to the SPICE system,
C          one or more FKs specifying these frames must be loaded.
C
C     The following data may be required:
C
C        - DSK data: if either SHAPE1 or SHAPE2 indicates that DSK
C          data are to be used, DSK files containing topographic data
C          for the target body must be loaded. If a surface list is
C          specified, data for at least one of the listed surfaces must
C          be loaded.
C
C        - Surface name-ID associations: if surface names are specified
C          in SHAPE1 or SHAPE2, the association of these names with
C          their corresponding surface ID codes must be established by
C          assignments of the kernel variables
C
C             NAIF_SURFACE_NAME
C             NAIF_SURFACE_CODE
C             NAIF_SURFACE_BODY
C
C          Normally these associations are made by loading a text
C          kernel containing the necessary assignments. An example
C          of such a set of assignments is
C
C             NAIF_SURFACE_NAME += 'Mars MEGDR 128 PIXEL/DEG'
C             NAIF_SURFACE_CODE += 1
C             NAIF_SURFACE_BODY += 499
C
C        - CK data: either of the body-fixed frames to which FRAME1 or
C          FRAME2 refer might be a CK frame. If so, at least one CK
C          file will be needed to permit transformation of vectors
C          between that frame and the J2000 frame.
C
C        - SCLK data: if a CK file is needed, an associated SCLK
C          kernel is required to enable conversion between encoded SCLK
C          (used to time-tag CK data) and barycentric dynamical time
C          (TDB).
C
C     Kernel data are normally loaded once per program run, NOT every
C     time this routine is called.
C
C$ Particulars
C
C     This routine supports the target shape combinations
C
C        POINT     - ELLIPSOID
C        POINT     - DSK
C        ELLIPSOID - ELLIPSOID
C
C     For many purposes, modeling extended bodies as triaxial
C     ellipsoids is adequate for determining whether one body is
C     occulted by another as seen from a specified observer.
C     
C
C     Using DSK data
C     ==============
C
C        DSK loading and unloading
C        -------------------------
C
C        DSK files providing data used by this routine are loaded by
C        calling FURNSH and can be unloaded by calling UNLOAD or
C        KCLEAR. See the documentation of FURNSH for limits on numbers
C        of loaded DSK files.
C
C        For run-time efficiency, it's desirable to avoid frequent
C        loading and unloading of DSK files. When there is a reason to
C        use multiple versions of data for a given target body---for
C        example, if topographic data at varying resolutions are to be
C        used---the surface list can be used to select DSK data to be
C        used for a given computation. It is not necessary to unload
C        the data that are not to be used. This recommendation presumes
C        that DSKs containing different versions of surface data for a
C        given body have different surface ID codes.
C
C
C        DSK data priority
C        -----------------
C
C        A DSK coverage overlap occurs when two segments in loaded DSK
C        files cover part or all of the same domain---for example, a
C        given longitude-latitude rectangle---and when the time
C        intervals of the segments overlap as well.
C
C        When DSK data selection is prioritized, in case of a coverage
C        overlap, if the two competing segments are in different DSK
C        files, the segment in the DSK file loaded last takes
C        precedence. If the two segments are in the same file, the
C        segment located closer to the end of the file takes
C        precedence.
C
C        When DSK data selection is unprioritized, data from competing
C        segments are combined. For example, if two competing segments
C        both represent a surface as sets of triangular plates, the
C        union of those sets of plates is considered to represent the
C        surface. 
C
C        Currently only unprioritized data selection is supported.
C        Because prioritized data selection may be the default behavior
C        in a later version of the routine, the UNPRIORITIZED keyword is
C        required in the SHAPE1 and SHAPE2 arguments.
C
C        
C        Syntax of the shape input arguments for the DSK case
C        ----------------------------------------------------
C
C        The keywords and surface list in the target shape arguments
C        SHAPE1 and SHAPE2, when DSK shape models are specified, are
C        called "clauses." The clauses may appear in any order, for
C        example
C
C           DSK/<surface list>/UNPRIORITIZED
C           DSK/UNPRIORITIZED/<surface list>
C           UNPRIORITIZED/<surface list>/DSK
C
C        The simplest form of a target argument specifying use of
C        DSK data is one that lacks a surface list, for example:
C
C           'DSK/UNPRIORITIZED'
C
C        For applications in which all loaded DSK data for the target
C        body are for a single surface, and there are no competing
C        segments, the above string suffices. This is expected to be
C        the usual case.
C
C        When, for the specified target body, there are loaded DSK
C        files providing data for multiple surfaces for that body, the
C        surfaces to be used by this routine for a given call must be
C        specified in a surface list, unless data from all of the
C        surfaces are to be used together.
C
C        The surface list consists of the string
C
C           SURFACES =
C
C        followed by a comma-separated list of one or more surface
C        identifiers. The identifiers may be names or integer codes in
C        string format. For example, suppose we have the surface
C        names and corresponding ID codes shown below:
C
C           Surface Name                              ID code
C           ------------                              -------
C           'Mars MEGDR 128 PIXEL/DEG'                1
C           'Mars MEGDR 64 PIXEL/DEG'                 2
C           'Mars_MRO_HIRISE'                         3
C
C        If data for all of the above surfaces are loaded, then
C        data for surface 1 can be specified by either
C
C           'SURFACES = 1'
C
C        or
C
C           'SURFACES = "Mars MEGDR 128 PIXEL/DEG"'
C
C        Double quotes are used to delimit the surface name because
C        it contains blank characters. 
C           
C        To use data for surfaces 2 and 3 together, any
C        of the following surface lists could be used:
C
C           'SURFACES = 2, 3'
C
C           'SURFACES = "Mars MEGDR  64 PIXEL/DEG", 3'
C
C           'SURFACES = 2, Mars_MRO_HIRISE'
C
C           'SURFACES = "Mars MEGDR 64 PIXEL/DEG", Mars_MRO_HIRISE'
C                  
C        An example of a shape argument that could be constructed
C        using one of the surface lists above is
C
C              'DSK/UNPRIORITIZED/SURFACES = '
C           // '"Mars MEGDR 64 PIXEL/DEG", 499003'
C
C
C$ Examples
C
C     1) Find whether MRO is occulted by Mars as seen by
C        the DSS-13 ground station at a few specific
C        times.
C
C        Use the meta-kernel shown below to load the required SPICE
C        kernels.
C
C
C          KPL/MK
C
C          File: occult_ex1.tm
C
C          This meta-kernel is intended to support operation of SPICE
C          example programs. The kernels shown here should not be
C          assumed to contain adequate or correct versions of data
C          required by SPICE-based user applications.
C
C          In order for an application to use this meta-kernel, the
C          kernels referenced here must be present in the user's
C          current working directory.
C
C          The names and contents of the kernels referenced
C          by this meta-kernel are as follows:
C
C             File name                        Contents
C             ---------                        --------
C             de410.bsp                        Planetary ephemeris
C             mar063.bsp                       Mars satellite ephemeris
C             pck00010.tpc                     Planet orientation and
C                                              radii
C             naif0011.tls                     Leapseconds
C             mro_psp35.bsp                    MRO ephemeris
C             megr90n000cb_plate.bds           Plate model based on
C                                              MEGDR DEM, resolution
C                                              4 pixels/degree.
C
C          \begindata
C
C             PATH_SYMBOLS    = ( 'MRO', 'GEN' )
C
C             PATH_VALUES     = (
C                                 '/ftp/pub/naif/pds/data+'
C                                 '/mro-m-spice-6-v1.0/+'
C                                 'mrosp_1000/data/spk',
C                                 '/ftp/pub/naif/generic_kernels'
C                               )
C
C             KERNELS_TO_LOAD = ( '$MRO/de410.bsp',
C                                 '$MRO/mar063.bsp',
C                                 '$MRO/mro_psp34.bsp',
C                                 '$GEN/spk/stations/+'
C                                 'earthstns_itrf93_050714.bsp',
C                                 '$GEN/pck/earth_latest_high_prec.bpc',
C                                 'pck00010.tpc',
C                                 'naif0011.tls',
C                                 'megr90n000cb_plate.bds'
C                               )
C          \begintext
C
C
C
C        Example code begins here.
C
C
C           PROGRAM OCCULT_MRO
C           IMPLICIT NONE
C
C           INCLUDE              'occult.inc'
C     C
C     C     Local parameters
C     C
C           CHARACTER*(*)         META
C           PARAMETER           ( META  = 'occult_ex1.tm' )
C
C           CHARACTER*(*)         FRMT
C           PARAMETER           ( FRMT  = '(A18, A5, A21, A5, A4, A6)' )
C
C           INTEGER               CHSIZ
C           PARAMETER           ( CHSIZ = 30 )
C
C     C
C     C     Local variables
C     C
C           CHARACTER*(CHSIZ)     ABCORR
C           CHARACTER*(CHSIZ)     FORM
C           CHARACTER*(CHSIZ)     OBSRVR
C           CHARACTER*(CHSIZ)     SHAPE1
C           CHARACTER*(CHSIZ)     SHAPE2
C           CHARACTER*(CHSIZ)     SHAPES ( 2 )
C           CHARACTER*(CHSIZ)     TARG1
C           CHARACTER*(CHSIZ)     TARG2
C           CHARACTER*(CHSIZ)     TIME
C           CHARACTER*(CHSIZ)     TSTART
C           CHARACTER*(CHSIZ)     TEND
C           CHARACTER*(CHSIZ)     OUTCH ( 4 )
C
C           DOUBLE PRECISION      ET
C           DOUBLE PRECISION      ET1
C           DOUBLE PRECISION      ETEND
C
C           INTEGER               DT
C           INTEGER               I
C           INTEGER               OCLTID
C
C     C
C     C     Saved variables
C     C
C           SAVE                  OUTCH
C           SAVE                  SHAPES
C     C
C     C     Initial values
C     C
C           DATA OUTCH ( 1 ) / 'totally occulted by'   /
C           DATA OUTCH ( 2 ) / 'transited by'          /
C           DATA OUTCH ( 3 ) / 'partially occulted by' /
C           DATA OUTCH ( 4 ) / 'not occulted by'       /
C
C           DATA SHAPES      / 'ELLIPSOID',
C          .                   'DSK/UNPRIORITIZED' /
C
C     C
C     C     Initialize the time range. Set the output time
C     C     format to PST. Set DT to an hour interval in
C     C     units of seconds.
C     C
C
C           TSTART = '2015-FEB-28 1:15:00 UTC'
C           TEND   = '2015-FEB-28 2:50:00 UTC'
C           FORM   = 'YYYY-MON-DD HR:MN ::UTC-8'
C           DT     = 1000
C
C     C
C     C     Initialize the targets, observer, and aberration
C     C     correction.
C     C
C           TARG1  = 'MRO'
C           SHAPE1 = 'POINT'
C           TARG2  = 'MARS'
C           OBSRVR = 'DSS-13'
C           ABCORR = 'CN'
C
C     C
C     C     Load kernel files via the meta-kernel.
C     C
C           CALL FURNSH ( META )
C     C
C     C     Calculate the start and stop times in ET.
C     C
C           CALL STR2ET ( TSTART, ET1   )
C           CALL STR2ET ( TEND,   ETEND )
C
C
C           DO I = 1, 2
C     C
C     C        Set the planet shape model for this pass.
C     C
C              SHAPE2 = SHAPES(I)
C
C              WRITE (*,*) ' '
C              CALL TOSTDO ( 'Mars shape: '//SHAPE2 )
C              WRITE (*,*) ' '
C
C              ET = ET1
C              DO WHILE ( ET .LT. ETEND )
C     C
C     C           Calculate the type of occultation that
C     C           corresponds to time ET.
C     C
C                 CALL OCCULT ( TARG1,  SHAPE1, ' ',
C          .                    TARG2,  SHAPE2, 'IAU_MARS',
C          .                    ABCORR, OBSRVR,  ET, OCLTID )
C     C
C     C           Output the results.
C     C
C                 CALL TIMOUT ( ET, FORM, TIME )
C
C                 IF ( OCLTID .EQ. TOTAL1 ) THEN
C                    WRITE (*,FRMT) TIME, TARG1, OUTCH(1), TARG2,
C          .                        'wrt ', OBSRVR
C
C                 ELSEIF ( OCLTID .EQ. ANNLR1 ) THEN
C                    WRITE (*,FRMT) TIME, TARG1, OUTCH(2), TARG2,
C          .                        'wrt ', OBSRVR
C
C                 ELSEIF ( OCLTID .EQ. PARTL1 ) THEN
C                    WRITE (*,FRMT) TIME, TARG1, OUTCH(3), TARG2,
C          .                        'wrt ', OBSRVR,
C          .                        'NOT POSSIBLE FOR POINT'
C
C                 ELSEIF ( OCLTID .EQ. NOOCC ) THEN
C                    WRITE (*,FRMT) TIME, TARG1, OUTCH(4), TARG2,
C          .                        'wrt ', OBSRVR
C
C                 ELSEIF ( OCLTID .EQ. PARTL2 ) THEN
C                    WRITE (*,FRMT) TIME, TARG2, OUTCH(3), TARG1,
C          .                        'wrt ', OBSRVR,
C          .                        'NOT POSSIBLE FOR POINT'
C
C                 ELSEIF ( OCLTID .EQ. ANNLR2 ) THEN
C                    WRITE (*,FRMT) TIME, TARG2, OUTCH(2), TARG1,
C          .                        'wrt ', OBSRVR
C
C                 ELSEIF ( OCLTID .EQ. TOTAL2 ) THEN
C                    WRITE (*,FRMT) TIME, TARG2, OUTCH(1), TARG1,
C          .                        'wrt ', OBSRVR
C
C                 ELSE
C                    WRITE (*,*) 'Bad occultation ID:  ', OCLTID
C
C                 END IF
C     C
C     C           Increment the time.
C     C
C                 ET = ET + DT
C
C              END DO
C
C           END DO
C           END
C
C
C        When this program was executed using gfortran on a PC Linux
C        64 bit environment, the output was:
C
C
C           Mars shape: ELLIPSOID
C
C           2015-FEB-27 17:15 MARS transited by         MRO  wrt DSS-13
C           2015-FEB-27 17:31 MRO  not occulted by      MARS wrt DSS-13
C           2015-FEB-27 17:48 MRO  totally occulted by  MARS wrt DSS-13
C           2015-FEB-27 18:04 MRO  totally occulted by  MARS wrt DSS-13
C           2015-FEB-27 18:21 MRO  not occulted by      MARS wrt DSS-13
C           2015-FEB-27 18:38 MARS transited by         MRO  wrt DSS-13
C
C           Mars shape: DSK/UNPRIORITIZED
C
C           2015-FEB-27 17:15 MARS transited by         MRO  wrt DSS-13
C           2015-FEB-27 17:31 MRO  not occulted by      MARS wrt DSS-13
C           2015-FEB-27 17:48 MRO  totally occulted by  MARS wrt DSS-13
C           2015-FEB-27 18:04 MRO  totally occulted by  MARS wrt DSS-13
C           2015-FEB-27 18:21 MRO  not occulted by      MARS wrt DSS-13
C           2015-FEB-27 18:38 MARS transited by         MRO  wrt DSS-13
C
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
C     S.C. Krening   (JPL)
C     N.J. Bachman   (JPL)
C
C$ Version 
C
C-    SPICELIB Version 2.0.0, 21-FEB-2017 (NJB)
C
C        Added FAILED tests.
C
C        01-MAR-2016 (NJB)
C
C        Upgraded to support surfaces represented by DSKs. Updated
C        example program to show use of DSKs. Updated example
C        meta-kernel. Corrected various comment typos.
C
C-    SPICELIB Version 1.0.0, 14-NOV-2013 (SCK) (NJB)
C
C
C-&

C$ Index_Entries
C
C     occultation type at a specified time
C
C-&

C
C     SPICELIB functions
C
      LOGICAL               FAILED
      LOGICAL               RETURN
C
C     External routines
C
      EXTERNAL              ZZGFOCIN
      EXTERNAL              ZZGFOCST
C
C     Local parameters
C
      INTEGER               BDNMLN
      PARAMETER           ( BDNMLN = 36 )

      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )
C
C     Local variables
C
      CHARACTER*(BDNMLN)    BACK
      CHARACTER*(BDNMLN)    BNAME
      CHARACTER*(FRNMLN)    BFRAME
      CHARACTER*(MTHLEN)    BMETHD
      CHARACTER*(MTHLEN)    BSHAPE
      CHARACTER*(MTHLEN)    FMETHD
      CHARACTER*(BDNMLN)    FNAME
      CHARACTER*(BDNMLN)    FRONT
      CHARACTER*(FRNMLN)    FFRAME
      CHARACTER*(MTHLEN)    FSHAPE
      CHARACTER*(SHPLEN)    OCCTYP ( 3 )
      CHARACTER*(CVTLEN)    PNTDEF
      CHARACTER*(SHPLEN)    PRSHP1
      CHARACTER*(SHPLEN)    PRSHP2
      CHARACTER*(MTHLEN)    SHAP1
      CHARACTER*(MTHLEN)    SHAP2
      CHARACTER*(SUBLEN)    SUBTYP
      CHARACTER*(TMTLEN)    TRMTYP

      INTEGER               I
      INTEGER               ID1
      INTEGER               ID2
      INTEGER               INDEX
      INTEGER               MLTFAC
      INTEGER               NSURF
      INTEGER               SRFLST ( MAXSRF )

      LOGICAL               ELLPS2
      LOGICAL               FOUND
      LOGICAL               OCSTAT
      LOGICAL               PRI
C
C     Saved variables
C
      SAVE                  OCCTYP

C
C     The variable OCCTYP associates the string of an occultation
C     type (from gf.inc) with its positive integer code (from
C     occult.inc).  The variable OCCTYP is set up so each string is
C     stored at the index relating to that configuration's positive
C     integer code.  The positive integer codes assume the first
C     target is occulting (in front of) the second target.
C
C                 Ex:  PARTL2 = 1               (from occult.inc)
C                      OCCTYP ( 1 ) = 'PARTIAL' (from gf.inc)
C
C     The table below shows the relation between each index of OCCTYP,
C     the occultation condition, which target is in front and back, the
C     multiplication factor, and the output integer occultation code.
C     Note that the output integer occultation code is the integer index
C     of OCCTYP multiplied by the multiplication factor.  
C
C                 OCLTID = Index * MLTFAC
C
C     MLTFAC is 1 if TARG1 is in front, and -1 if TARG1 is in back.
C     The setup of OCCTYP could be changed, but it is important to keep
C     the output integer occultation codes consistent with the values
C     from occult.inc.
C
C         Index   Occult. Condition   Front   Back   MLTFAC  OCLTID
C         -----   -----------------   -----   -----  ------  ------
C           1          Partial        TARG1   TARG2    1       1
C           1          Partial        TARG2   TARG1   -1      -1
C           2          Annular        TARG1   TARG2    1       2
C           2          Annular        TARG2   TARG1   -1      -2
C           3          Total          TARG1   TARG2    1       3
C           3          Total          TARG2   TARG1   -1      -3
C

      DATA OCCTYP ( PARTL2 ) / PARTL  /
      DATA OCCTYP ( ANNLR2 ) / ANNULR /
      DATA OCCTYP ( TOTAL2 ) / FULL   /

C
C     Standard SPICE error handling.
C
      IF ( RETURN() ) THEN
         RETURN
      END IF
      
      CALL CHKIN ( 'OCCULT' )

C
C     Left justify the shapes and target names and make them upper case.
C
      CALL LJUST ( SHAPE1, SHAP1 )
      CALL UCASE ( SHAP1,  SHAP1 )

      CALL LJUST ( SHAPE2, SHAP2 )
      CALL UCASE ( SHAP2,  SHAP2 )

      CALL LJUST ( TARG1,  FNAME )
      CALL UCASE ( FNAME,  FNAME )

      CALL LJUST ( TARG2,  BNAME )
      CALL UCASE ( BNAME,  BNAME )

C
C     The variable ELLPS2 is a flag that indicates whether both targets
C     are represented as ellipsoids.
C
      ELLPS2 =   ( SHAP1 .EQ. EDSHAP )  .AND.
     .           ( SHAP2 .EQ. EDSHAP ) 

C
C     Parse the input shapes. We need the target ID codes
C     for this.
C
      IF ( SHAP1 .EQ. 'POINT' ) THEN

         PRSHP1 = SHAP1(:SHPLEN)

      ELSE

         CALL BODS2C ( FNAME, ID1, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'First target name # could not be '
     .      //            'mapped to an ID code.'            )
            CALL ERRCH  ( '#', FNAME                         )
            CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'            )
            CALL CHKOUT ( 'OCCULT'                           )
            RETURN

         END IF

         CALL ZZPRSMET ( ID1, SHAP1, MAXSRF, PRSHP1, SUBTYP,
     .                   PRI, NSURF, SRFLST, PNTDEF, TRMTYP )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'OCCULT' )
            RETURN
         END IF

      END IF

      IF ( SHAP2 .EQ. 'POINT' ) THEN

         PRSHP2 = SHAP2(:SHPLEN)

      ELSE

         CALL BODS2C ( BNAME, ID2, FOUND )

         IF ( .NOT. FOUND ) THEN

            CALL SETMSG ( 'Second target name # could not be '
     .      //            'mapped to an ID code.'            )
            CALL ERRCH  ( '#', BNAME                         )
            CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'            )
            CALL CHKOUT ( 'OCCULT'                           )
            RETURN

         END IF

         CALL ZZPRSMET ( ID2, SHAP2, MAXSRF, PRSHP2, SUBTYP,
     .                   PRI, NSURF, SRFLST, PNTDEF, TRMTYP )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'OCCULT' )
            RETURN
         END IF

      END IF

C
C     Test two main cases:
C
C     1) The first target is the front body.
C     2) The second target is the front body.
C
C     First, initialize the occultation code to reflect no occultation.
C
      OCLTID = NOOCC

      DO I = 1, 2 
C
C        The first time through, make the first target the
C        front. On the second time, make the second target the front.
C        For details on the variable MLTFAC, please see the detailed
C        explanation of the OCCTYP variable near the start of the code.
C
         IF ( I .EQ. 1 ) THEN

            FRONT  = FNAME
            FMETHD = SHAP1
            BMETHD = SHAP2
            FSHAPE = PRSHP1
            FFRAME = FRAME1
            BACK   = BNAME
            BSHAPE = PRSHP2
            BFRAME = FRAME2
            MLTFAC = 1
            
         ELSE 
            
            FRONT  = BNAME
            FMETHD = SHAP2
            BMETHD = SHAP1
            FSHAPE = PRSHP2
            FFRAME = FRAME2
            BACK   = FNAME
            BSHAPE = PRSHP1
            BFRAME = FRAME1
            MLTFAC = -1
            
         END IF
    
C
C        Check if there is any occultation with the current front/back
C        configuration. ZZGFOCIN performs initializations. ZZGFOCST
C        returns a true/false logical indicating if there is an
C        occultation.
C
         CALL ZZGFOCIN ( ANY, FRONT,  FMETHD, FFRAME, 
     .                        BACK,   BMETHD, BFRAME,  
     .                        OBSRVR, ABCORR        )
         
         CALL ZZGFOCST ( ET, OCSTAT )
         
         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'OCCULT' )
            RETURN
         END IF
         
C
C        If there was an occultation, and both targets are represented
C        as ellipsoids, test the three types of occultations: partial,
C        annular, and full. Note: If the integer parameters within
C        occult.inc are changed, the following DO loop will need to be
C        altered.
C
         IF ( OCSTAT ) THEN
C
C           An occultation exists.
C
            IF ( ELLPS2 ) THEN
C
C              Both shapes are ellipsoids.
C              
               DO  INDEX = PARTL2, TOTAL2
               
                  CALL ZZGFOCIN ( OCCTYP(INDEX), 
     .                            FRONT,  FSHAPE, FFRAME,
     .                            BACK,   BSHAPE, BFRAME, 
     .                            OBSRVR, ABCORR        )
                              
                  CALL ZZGFOCST ( ET, OCSTAT )
               
                  IF ( FAILED() ) THEN
                     CALL CHKOUT ( 'OCCULT' )
                     RETURN
                  END IF
C
C                 If the occultation condition is true, return the
C                 integer occultation ID code.
C
                  IF ( OCSTAT ) THEN
                  
                     OCLTID = MLTFAC * INDEX

                     CALL CHKOUT ( 'OCCULT' )
                     RETURN

                  END IF      
C
C                 End the DO loop that checks the occultation type.
C         
               END DO


            ELSE IF (      ( FSHAPE .EQ. EDSHAP ) 
     .                .OR. ( FSHAPE .EQ. DSSHAP )  ) THEN
C
C              The front target is an ellipsoid or DSK shape: this
C              is a total occultation. (Other target is a point).
C           
               OCLTID = MLTFAC * TOTAL2

               CALL CHKOUT ( 'OCCULT' )
               RETURN


            ELSE IF (      ( BSHAPE .EQ. EDSHAP ) 
     .                .OR. ( BSHAPE .EQ. DSSHAP )  ) THEN
C
C              The back target is an ellipsoid or DSK shape: this is an
C              annular occultation. (Other target is a point).
C   
               OCLTID = MLTFAC * ANNLR2

               CALL CHKOUT ( 'OCCULT' )
               RETURN

            END IF

         END IF
C
C        End the DO loop that checks the front/back orientation of
C        the input targets.
C
      END DO

C
C     If the occultation searches show that there was no occultation
C     at the given time, return an occultation code that indicates
C     no occultation. If this part of the code has been reached and 
C     the occultation code indicates an occultation was found, an error
C     has occurred.  
C     
      IF ( OCLTID .NE. NOOCC ) THEN
         
         CALL SETMSG ( 'This error should never be reached; '
     .        //       'the occultation code result # is invalid.')
         CALL ERRINT ( '#', OCLTID  )
         CALL SIGERR ( 'SPICE(BUG)' )
         CALL CHKOUT ( 'OCCULT'      )
         RETURN
         
      END IF
      
      CALL CHKOUT ( 'OCCULT' )
      RETURN
      END

