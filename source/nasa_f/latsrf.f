C$Procedure LATSRF ( Latitudinal grid to surface points )

      SUBROUTINE LATSRF ( METHOD, TARGET, ET, 
     .                    FIXREF, NPTS,   LONLAT, SRFPTS )

C$ Abstract
C
C     Map array of planetocentric longitude/latitude coordinate pairs
C     to surface points on a specified target body.
C
C     The surface of the target body may be represented by a triaxial
C     ellipsoid or by topographic data provided by DSK files.
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
C     FRAMES
C     PCK
C     SPK
C     TIME
C
C$ Keywords
C
C     COORDINATES
C     DSK
C     GEOMETRY
C     SURFACE
C
C$ Declarations
  
      IMPLICIT NONE

      INCLUDE 'dsk.inc'
      INCLUDE 'gf.inc'
      INCLUDE 'zzctr.inc'
      INCLUDE 'zzdsk.inc'

      CHARACTER*(*)         METHOD
      CHARACTER*(*)         TARGET
      DOUBLE PRECISION      ET
      CHARACTER*(*)         FIXREF
      INTEGER               NPTS
      DOUBLE PRECISION      LONLAT ( 2, * )
      DOUBLE PRECISION      SRFPTS ( 3, * )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     METHOD     I   Computation method.
C     TARGET     I   Name of target body.
C     ET         I   Epoch in TDB seconds past J2000 TDB.
C     FIXREF     I   Body-fixed, body-centered target body frame.
C     NPTS       I   Number of coordinate pairs in input array.
C     LONLAT     I   Array of longitude/latitude coordinate pairs.
C     SRFPTS     O   Array of surface points.
C
C$ Detailed_Input
C
C
C     METHOD      is a short string providing parameters defining
C                 the computation method to be used. In the syntax
C                 descriptions below, items delimited by brackets
C                 are optional.
C                
C                 METHOD may be assigned the following values:   
C
C                    'ELLIPSOID'
C 
C                       The surface point computation uses a triaxial
C                       ellipsoid to model the surface of the target
C                       body. The ellipsoid's radii must be available
C                       in the kernel pool.
C
C                    'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                       The surface point computation uses topographic
C                       data to model the surface of the target body.
C                       These data must be provided by loaded DSK
C                       files.
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
C
C                 Neither case nor white space are significant in
C                 METHOD, except within double-quoted strings. For
C                 example, the string ' eLLipsoid ' is valid.
C
C                 Within double-quoted strings, blank characters are
C                 significant, but multiple consecutive blanks are
C                 considered equivalent to a single blank. Case is 
C                 not significant. So
C
C                    "Mars MEGDR 128 PIXEL/DEG"
C
C                 is equivalent to 
C
C                    " mars megdr  128  pixel/deg "
C
C                 but not to
C
C                    "MARS MEGDR128PIXEL/DEG"
C 
C                 
C     TARGET      is the name of the target body. TARGET is
C                 case-insensitive, and leading and trailing blanks in
C                 TARGET are not significant. Optionally, you may
C                 supply a string containing the integer ID code for
C                 the object. For example both 'MOON' and '301' are
C                 legitimate strings that indicate the Moon is the
C                 target body.
C
C                 When the target body's surface is represented by a
C                 tri-axial ellipsoid, this routine assumes that a
C                 kernel variable representing the ellipsoid's radii is
C                 present in the kernel pool. Normally the kernel
C                 variable would be defined by loading a PCK file.
C
C
C     ET          is the epoch for which target surface data will be
C                 selected, if the surface is modeled using DSK data.
C                 In this case, only segments having time coverage that
C                 includes the epoch ET will be used.
C
C                 ET is ignored if the target is modeled as an
C                 ellipsoid.
C
C                 ET is expressed as TDB seconds past J2000 TDB.
C
C                                  
C     FIXREF      is the name of a body-fixed reference frame centered
C                 on the target body.
C
C                 If the target shape is given by DSK data, FIXREF may
C                 designate any such frame supported by the SPICE
C                 system, including built-in frames (documented in the
C                 Frames Required Reading) and frames defined by a
C                 loaded frame kernel (FK).
C
C                 When the target surface is modeled as an ellipsoid,
C                 the reference frame designated by FIXREF (described
C                 below) must have its coordinate axes aligned with the
C                 respective principal axes of the reference ellipsoid.
C
C                 The string FIXREF is case-insensitive, and leading
C                 and trailing blanks in FIXREF are not significant.
C
C                 The output surface points in the array SRFPTS will be
C                 expressed relative to this reference frame.
C
C
C     NPTS        is the number of coordinate pairs in the array LONLAT.
C
C
C     LONLAT      is an array of pairs of planetocentric longitudes and
C                 latitudes of surface points. Elements
C
C                    LONLAT(1,I)
C                    LONLAT(2,I)
C
C                 are, respectively, the planetocentric longitude and
C                 latitude of the Ith surface point.
C
C                 The units of longitude and latitude are radians.
C
C
C$ Detailed_Output
C
C     SRFPTS      is an array of target body surface points
C                 corresponding to the pairs of coordinates in the
C                 input LONLAT array. Elements
C
C                    SRFPTS(1,I)
C                    SRFPTS(2,I)
C                    SRFPTS(3,I)
C
C                 are the Cartesian coordinates, expressed in the
C                 reference frame designated by FIXREF, of the surface
C                 point corresponding to the Ith pair of input
C                 coordinates.
C
C                 If there are multiple solutions for a given input
C                 coordinate pair, this routine will return the point
C                 at those coordinates having the greatest distance
C                 from the origin of the coordinate system.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C
C     1)  If the target body name
C         input string cannot be converted to an integer ID code, the
C         error SPICE(IDCODENOTFOUND) is signaled.
C
C     2)  If the input target body-fixed frame FIXREF is not
C         recognized, the error SPICE(NOFRAME) is signaled. A frame
C         name may fail to be recognized because a required frame
C         specification kernel has not been loaded; another cause is a
C         misspelling of the frame name.
C
C     3)  If the input frame FIXREF is not centered at the target body,
C         the error SPICE(INVALIDFRAME) is signaled.
C
C     4)  If data are not available to convert between the frame
C         FIXREF and the frame of a DSK segment of interest, the error
C         will be signaled by a routine in the call tree of this
C         routine.
C
C     5)  If the input argument METHOD cannot be parsed, the error
C         will be signaled either by this routine or by a routine in
C         the call tree of this routine.
C
C     6)  If the computation method specifies an ellipsoidal target
C         model, and if triaxial radii of the target body have not been
C         loaded into the kernel pool prior to calling LATSRF, the
C         error will be diagnosed and signaled by a routine in the call
C         tree of this routine.
C
C     7)  The target must be an extended body: if the computation
C         method specifies an ellipsoidal target model, and if any of
C         the radii of the target body are non-positive, the error will
C         be signaled by routines in the call tree of this routine.
C
C     8)  If METHOD specifies that the target surface is represented by
C         DSK data, and no DSK files are loaded for the specified
C         target, the error is signaled by a routine in the call tree
C         of this routine.
C         
C     9)  If METHOD specifies that the target surface is represented
C         by DSK data, and data representing the portion of the surface
C         corresponding to the coordinates provided in LONLAT are not
C         available, an error will be signaled by a routine in the call
C         tree of this routine.
C
C    10)  If a surface point cannot be computed because the ray
C         corresponding to a longitude/latitude pair fails to intersect
C         the target surface as defined by the plate model, the error
C         SPICE(NOINTERCEPT) is signaled.
C
C    11)  If the surface point corresponding to a longitude/latitude
C         pair in LONLAT does not have matching longitude and latitude
C         (because it is on the opposite side of the origin), the error
C         SPICE(SHAPENOTSUPPORTED) is signaled.
C
C    12)  If the target shape is "ellipsoid" and not all radii of the
C         ellipsoid are strictly positive, the error
C         SPICE(BADAXISLENGTH) will be signaled. If the radii are not
C         available in the kernel pool, the error will be signaled
C         by a routine in the call tree of this routine.
C
C$ Files
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - Shape data for the target body:
C                
C            PCK data: 
C
C               If the target shape is modeled as an ellipsoid,
C               triaxial radii for the target body must be loaded into
C               the kernel pool. Typically this is done by loading a
C               text PCK file via FURNSH.
C
C            DSK data: 
C
C               If the target shape is modeled by DSK data, DSK files
C               containing topographic data for the target body must be
C               loaded. If a surface list is specified, data for at
C               least one of the listed surfaces must be loaded.
C
C        - Target body orientation data: these may be provided in a
C          text or binary PCK file. In some cases, target body
C          orientation may be provided by one more more CK files. In
C          either case, data are made available by loading the files
C          via FURNSH.
C
C     The following data may be required:
C
C        - Frame data: if a frame definition is required to convert 
C          between the body-fixed frame of the target and the frame of
C          a DSK segment providing topographic data, that definition
C          must be available in the kernel pool. Typically the
C          definition is supplied by loading a frame kernel via FURNSH.
C
C        - Surface name-ID associations: if surface names are specified
C          in METHOD, the association of these names with their
C          corresponding surface ID codes must be established by 
C          assignments of the kernel variables
C
C             NAIF_SURFACE_NAME
C             NAIF_SURFACE_CODE
C             NAIF_SURFACE_BODY
C
C          Normally these associations are made by loading a text
C          kernel containing the necessary assignments. An example of
C          such a set of assignments is
C
C             NAIF_SURFACE_NAME += 'Mars MEGDR 128 PIXEL/DEG'
C             NAIF_SURFACE_CODE += 1                    
C             NAIF_SURFACE_BODY += 499
C
C        - SCLK data: if the target body's orientation is provided by
C          CK files, an associated SCLK kernel must be loaded.
C
C     In all cases, kernel data are normally loaded once per program
C     run, NOT every time this routine is called. 
C
C
C$ Particulars
C
C     This routine is intended to be used for target body surfaces that
C     have a unique radius for each pair of planetocentric longitude
C     and latitude coordinates.
C
C     If the target surface is represented by topographic data, it is
C     possible for there to be multiple surface points at a given
C     planetocentric longitude and latitude. For example, this can
C     occur if the surface has features such as cliffs, caves, or
C     arches.
C
C     For more complex surfaces, the routine
C
C        DSKSXV {DSK, ray-surface intercept, vectorized}
C
C     may be more suitable. That routine works with rays having vertices
C     anywhere outside of the target body.
C
C
C     Planetocentric coordinates
C     ==========================
C
C     Planetocentric longitude and latitude are defined as follows:
C
C        Longitude of a point P is the angle between the prime meridian
C        and the meridian containing P. The direction of increasing
C        longitude is from the +X axis towards the +Y axis.
C 
C        Latitude of a point P is the angle from the XY plane of the
C        ray from the origin through the point.
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
C        required in the METHOD argument.
C
C        
C        Syntax of the METHOD input argument
C        -----------------------------------
C
C        The keywords and surface list in the METHOD argument
C        are called "clauses." The clauses may appear in any
C        order, for example
C
C           DSK/<surface list>/UNPRIORITIZED
C           DSK/UNPRIORITIZED/<surface list>
C           UNPRIORITIZED/<surface list>/DSK
C
C        The simplest form of the METHOD argument specifying use of
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
C        An example of a METHOD argument that could be constructed
C        using one of the surface lists above is
C
C           'DSK/UNPRIORITIZED/SURFACES = "Mars MEGDR 64 PIXEL/DEG", 3'
C
C 
C$ Examples
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as input,
C     the compiler and supporting libraries, and the machine specific
C     arithmetic implementation.
C
C
C     1)  In the following example program, the file
C
C            phobos512.bds
C
C         is a DSK file containing a type 2 segment that provides a
C         plate model representation of the surface of Phobos.
C
C         Find the surface points on a target body corresponding to a
C         given planetocentric longitude/latitude grid. In order to
C         duplicate the example output, the kernel name
C
C            phobos512.bds
C
C         should be supplied at the prompt.
C                
C
C         Example code begins here.
C
C
C              PROGRAM EX1
C              IMPLICIT NONE
C        C
C        C     SPICELIB functions
C        C
C              DOUBLE PRECISION      DPR
C              DOUBLE PRECISION      RPD
C        C
C        C     Local parameters
C        C
C              CHARACTER*(*)         FMT1
C              PARAMETER           ( FMT1   = '(1X,A,F11.6)' )
C
C              INTEGER               BDNMLN
C              PARAMETER           ( BDNMLN = 36 )
C
C              INTEGER               FILSIZ
C              PARAMETER           ( FILSIZ = 255 )
C
C              INTEGER               FRNMLN
C              PARAMETER           ( FRNMLN = 32 )
C
C              INTEGER               LNSIZE
C              PARAMETER           ( LNSIZE = 79 )
C
C              INTEGER               MAXN
C              PARAMETER           ( MAXN   = 100 )
C
C              INTEGER               MTHLEN
C              PARAMETER           ( MTHLEN = 80 )
C        C
C        C     Local variables
C        C
C              CHARACTER*(FILSIZ)    DSK
C              CHARACTER*(FRNMLN)    FIXREF
C              CHARACTER*(MTHLEN)    METHOD
C              CHARACTER*(LNSIZE)    OUTLIN
C              CHARACTER*(BDNMLN)    TARGET
C
C              DOUBLE PRECISION      DLAT
C              DOUBLE PRECISION      DLON
C              DOUBLE PRECISION      ET
C              DOUBLE PRECISION      GRID   ( 2, MAXN )
C              DOUBLE PRECISION      LAT
C              DOUBLE PRECISION      LAT0
C              DOUBLE PRECISION      LON
C              DOUBLE PRECISION      LON0
C              DOUBLE PRECISION      SRFPTS ( 3, MAXN )
C              DOUBLE PRECISION      XLAT
C              DOUBLE PRECISION      XLON
C              DOUBLE PRECISION      XR
C
C              INTEGER               I
C              INTEGER               J
C              INTEGER               N
C              INTEGER               NLAT
C              INTEGER               NLON
C
C        C
C        C     Set target, reference frame, and epoch.
C        C
C              TARGET = 'PHOBOS'
C              FIXREF = 'IAU_PHOBOS'
C              ET     = 0.D0
C        C
C        C     Use DSK data to represent the surface.
C        C
C              METHOD = 'DSK/UNPRIORITIZED'
C        C
C        C     Set the grid dimensions.
C        C
C              NLON   = 6
C              NLAT   = 3
C        C
C        C     Derive evenly spaced grid separations and starting
C        C     values in the longitude and latitude dimensions.
C        C     Units are degrees.
C        C
C              LAT0 = 90.D0
C              LON0 =  0.D0
C
C              DLAT = 180.D0 / (NLAT + 1)
C              DLON = 360.D0 /  NLON
C        C
C        C     Prompt for the name of the DSK to read.
C        C
C              CALL PROMPT ( 'Enter DSK name    > ', DSK )
C        C
C        C     Load the DSK file.
C        C
C              CALL FURNSH ( DSK )
C        C
C        C     Now generate the grid points.  We generate
C        C     points along latitude bands, working from
C        C     north to south.  The latitude range is selected
C        C     to range from +45 to -45 degrees.  Longitude
C        C     ranges from 0 to 300 degrees.  The increment
C        C     is 45 degrees for latitude and 60 degrees for
C        C     longitude.
C        C
C              N = 0
C
C              DO I = 1, NLAT
C
C                 LAT = RPD() * ( LAT0 - I*DLAT )
C
C                 DO J = 1, NLON
C
C                    N   = N + 1
C                    LON = RPD() * ( LON0 + (J-1)*DLON )
C
C                    GRID(1,N) = LON
C                    GRID(2,N) = LAT
C
C                 END DO
C
C              END DO
C        C
C        C     Find the surface points corresponding to the grid points.
C        C
C              CALL LATSRF ( METHOD, TARGET, ET,
C             .              FIXREF, N,      GRID, SRFPTS )
C        C
C        C     Print out the surface points in latitudinal
C        C     coordinates and compare the derived lon/lat values
C        C     to those of the input grid.
C        C
C              DO I = 1, N
C        C
C        C        Use RECRAD rather than RECLAT to produce
C        C        non-negative longitudes.
C        C
C                 CALL RECRAD ( SRFPTS(1,I), XR, XLON, XLAT )
C
C                 WRITE (*,*) ' '
C
C                 OUTLIN = 'Intercept for grid point #:'
C                 CALL REPMI ( OUTLIN, '#', I, OUTLIN )
C
C                 WRITE(*,*)  OUTLIN
C                 OUTLIN = '  Cartesian coordinates: (#, #, #)'
C
C                 DO J = 1, 3
C                    CALL REPMF( OUTLIN, '#', SRFPTS(J,I),
C             .                  5,      'E', OUTLIN      )
C                 END DO
C
C                 WRITE (*,*) OUTLIN
C
C                 WRITE (*,*)    '  Latitudinal Coordinates:'
C                 WRITE (*,FMT1) '   Longitude (deg): ', XLON*DPR()
C                 WRITE (*,FMT1) '   Latitude  (deg): ', XLAT*DPR()
C                 WRITE (*,FMT1) '   Radius     (km): ', XR
C                 WRITE (*,*)    ' '
C                 WRITE (*,*)    '  Original Grid Coordinates:'
C                 WRITE (*,FMT1) '   Longitude (deg): ', GRID(1,I)*DPR()
C                 WRITE (*,FMT1) '   Latitude  (deg): ', GRID(2,I)*DPR()
C
C              END DO
C              WRITE (*,*) ' '
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit
C     platform, the output for the first 3 points and the last 3 points
C     (the rest of the output is not shown due to its large volume)
C     was:
C
C
C      Enter DSK name    > phobos512.bds
C
C       Intercept for grid point 1:
C         Cartesian coordinates: (7.1817E+00, 0.0000E+00, 7.1817E+00)
C         Latitudinal Coordinates:
C          Longitude (deg):    0.000000
C          Latitude  (deg):   45.000000
C          Radius     (km):   10.156402
C
C         Original Grid Coordinates:
C          Longitude (deg):    0.000000
C          Latitude  (deg):   45.000000
C
C       Intercept for grid point 2:
C         Cartesian coordinates: (3.5820E+00, 6.2042E+00, 7.1640E+00)
C         Latitudinal Coordinates:
C          Longitude (deg):   60.000000
C          Latitude  (deg):   45.000000
C          Radius     (km):   10.131412
C
C         Original Grid Coordinates:
C          Longitude (deg):   60.000000
C          Latitude  (deg):   45.000000
C
C       Intercept for grid point 3:
C         Cartesian coordinates: (-3.6854E+00, 6.3832E+00, 7.3707E+00)
C         Latitudinal Coordinates:
C          Longitude (deg):  120.000000
C          Latitude  (deg):   45.000000
C          Radius     (km):   10.423766
C
C         Original Grid Coordinates:
C          Longitude (deg):  120.000000
C          Latitude  (deg):   45.000000
C
C           ...
C
C       Intercept for grid point 16:
C         Cartesian coordinates: (-8.2374E+00, 1.5723E-15, -8.2374E+00)
C         Latitudinal Coordinates:
C          Longitude (deg):  180.000000
C          Latitude  (deg):  -45.000000
C          Radius     (km):   11.649512
C
C         Original Grid Coordinates:
C          Longitude (deg):  180.000000
C          Latitude  (deg):  -45.000000
C
C       Intercept for grid point 17:
C         Cartesian coordinates: (-3.6277E+00, -6.2833E+00, -7.2553E+00)
C         Latitudinal Coordinates:
C          Longitude (deg):  240.000000
C          Latitude  (deg):  -45.000000
C          Radius     (km):   10.260572
C
C         Original Grid Coordinates:
C          Longitude (deg):  240.000000
C          Latitude  (deg):  -45.000000
C
C       Intercept for grid point 18:
C         Cartesian coordinates: (3.2881E+00, -5.6952E+00, -6.5762E+00)
C         Latitudinal Coordinates:
C          Longitude (deg):  300.000000
C          Latitude  (deg):  -45.000000
C          Radius     (km):    9.300154
C
C         Original Grid Coordinates:
C          Longitude (deg):  300.000000
C          Latitude  (deg):  -45.000000
C
C
C$ Restrictions
C
C     1)  This routine assumes that the origin of the body-fixed
C         reference frame associated with the target body is located in
C         the interior of that body.
C
C     2)  The results returned by this routine may not be meaningful
C         if the target surface has multiple surface points associated
C         with some (longitude, latitude) coordinates.
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
C-    SPICELIB Version 1.0.0, 21-FEB-2017 (NJB) 
C
C        Original version 01-JUL-2016 (NJB) 
C
C-&

C$ Index_Entries
C
C     map latitudinal coordinates to Cartesian surface points
C     map latitudinal coordinates to DSK surface points
C
C-&


C
C     SPICELIB functions
C
      DOUBLE PRECISION      DPR
      DOUBLE PRECISION      VDOT

      LOGICAL               EQSTR
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      CHARACTER*(*)         RNAME
      PARAMETER           ( RNAME = 'LATSRF' )

      INTEGER               BDNMLN
      PARAMETER           ( BDNMLN = 36 )

      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )

C
C     Local variables
C
      CHARACTER*(CVTLEN)    LMBTYP
      CHARACTER*(MTHLEN)    PRVMTH
      CHARACTER*(SHPLEN)    SHPSTR
      CHARACTER*(SUBLEN)    SUBTYP
      CHARACTER*(TMTLEN)    TRMTYP

      DOUBLE PRECISION      R
      DOUBLE PRECISION      RAYDIR ( 3 )
      DOUBLE PRECISION      X      ( 3 )

      INTEGER               FIXFID
      INTEGER               FXCENT
      INTEGER               FXCLSS
      INTEGER               FXCLID
      INTEGER               I
      INTEGER               N
      INTEGER               NSURF
      INTEGER               SHAPE
      INTEGER               SRFCTR ( CTRSIZ )
      INTEGER               SRFLST ( MAXSRF )
      INTEGER               TRGCDE

      LOGICAL               FIRST
      LOGICAL               FND
      LOGICAL               PRI
      LOGICAL               SURFUP

C
C     Saved name/ID item declarations.
C
      INTEGER               SVCTR1 ( CTRSIZ )
      CHARACTER*(BDNMLN)    SVTARG
      INTEGER               SVTCDE
      LOGICAL               SVFND1

C
C     Saved frame name/ID item declarations.
C
      INTEGER               SVCTR2 ( CTRSIZ )
      CHARACTER*(FRNMLN)    SVFREF
      INTEGER               SVFXFC

C
C     Saved target radius values.
C
      INTEGER               SVCTR3 ( CTRSIZ )
      INTEGER               SVPRVT
      DOUBLE PRECISION      SVRADI ( 3 )

C
C     Saved name/ID items.
C
      SAVE                  SVCTR1
      SAVE                  SVTARG
      SAVE                  SVTCDE
      SAVE                  SVFND1
C
C     Saved frame name/ID items.
C
      SAVE                  SVCTR2
      SAVE                  SVFREF
      SAVE                  SVFXFC

C
C     Saved target radius values.
C
      SAVE                  SVCTR3 
      SAVE                  SVPRVT
      SAVE                  SVRADI

C
C     Saved values
C
      SAVE                  NSURF
      SAVE                  FIRST
      SAVE                  PRI
      SAVE                  PRVMTH
      SAVE                  SHAPE
      SAVE                  SRFCTR
      SAVE                  SRFLST

C
C     Initial values
C
      DATA                  FIRST  / .TRUE. /
      DATA                  PRVMTH / ' '    /
      DATA                  SVPRVT / 0      /



      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( RNAME )

      IF ( FIRST ) THEN
C
C        Initialize local surface counter.
C
         CALL ZZCTRUIN ( SRFCTR )
C
C        Initialize target, frame, and radius counters.
C
         CALL ZZCTRUIN ( SVCTR1 )
         CALL ZZCTRUIN ( SVCTR2 )
         CALL ZZCTRUIN ( SVCTR3 )

      END IF

C
C     Obtain integer code for the target.
C 
      CALL ZZBODS2C ( SVCTR1, SVTARG, SVTCDE, SVFND1,
     .                TARGET, TRGCDE, FND    )
            
      IF ( .NOT. FND ) THEN
      
         CALL SETMSG ( 'The target, '
     .   //            '''#'', is not a recognized name for an '
     .   //            'ephemeris object. The cause of this '
     .   //            'problem may be that you need an updated '
     .   //            'version of the SPICE Toolkit, or that you '
     .   //            'failed to load a kernel containing a '
     .   //            'name-ID mapping for this body.'           )
         CALL ERRCH  ( '#', TARGET                                )
         CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                    )
         CALL CHKOUT ( RNAME                                      )
         RETURN
      
      END IF

C
C     Determine the attributes of the frame designated by FIXREF.
C
      CALL ZZNAMFRM ( SVCTR2, SVFREF, SVFXFC, FIXREF, FIXFID )

      CALL FRINFO ( FIXFID, FXCENT, FXCLSS, FXCLID, FND )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

      IF ( .NOT. FND ) THEN

         CALL SETMSG ( 'Reference frame # is not recognized by ' //
     .                 'the SPICE frame subsystem. Possibly '    //
     .                 'a required frame definition kernel has ' //
     .                 'not been loaded.'                        )
         CALL ERRCH  ( '#',  FIXREF                              )
         CALL SIGERR ( 'SPICE(NOFRAME)'                          )
         CALL CHKOUT ( RNAME                                     )
         RETURN

      END IF

C
C     Make sure that FIXREF is centered at the target body's center.
C
      IF ( FXCENT .NE. TRGCDE ) THEN

         CALL SETMSG ( 'Reference frame # is not centered at the ' 
     .   //            'target body #. The ID code of the frame '
     .   //            'center is #.'                             )
         CALL ERRCH  ( '#',  FIXREF                               )
         CALL ERRCH  ( '#',  TARGET                               )
         CALL ERRINT ( '#',  FXCENT                               )
         CALL SIGERR ( 'SPICE(INVALIDFRAME)'                      )
         CALL CHKOUT ( RNAME                                      )
         RETURN

      END IF

C
C     Check whether the surface name/ID mapping has been updated.
C
      CALL ZZSRFTRK ( SRFCTR, SURFUP )

C
C     Initialize the SINCPT utility package for the next computation.
C     The choice of initialization routine depends on the target
C     surface type.
C     
      IF ( FIRST  .OR.  SURFUP  .OR.  ( METHOD .NE. PRVMTH ) ) THEN
C
C        Set the previous method string to an invalid value, so it
C        cannot match any future, valid input. This will force this
C        routine to parse the input method on the next call if any
C        failure occurs in this branch. Once success is assured, we can
C        record the current method in the previous method string.
C
         PRVMTH = ' '

C
C        Parse the method string. If the string is valid, the
C        outputs SHAPE and SUBTYP will always be be set. However,
C        SUBTYP is not used in this routine. 
C
C        For DSK shapes, the surface list array and count will be set
C        if the method string contains a surface list.
C
         CALL ZZPRSMET ( TRGCDE, METHOD, MAXSRF, SHPSTR, SUBTYP,
     .                   PRI,    NSURF,  SRFLST, LMBTYP, TRMTYP )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         IF (  EQSTR( SHPSTR, 'ELLIPSOID' )  ) THEN

            SHAPE = ELLSHP

         ELSE IF (  EQSTR( SHPSTR, 'DSK' )  ) THEN

            SHAPE = DSKSHP

         ELSE
C
C           This is a backstop check.
C
            CALL SETMSG ( '[1] Returned shape value from method '
     .      //            'string was <#>.'                      )
            CALL ERRCH  ( '#', SHPSTR                            )
            CALL SIGERR ( 'SPICE(BUG)'                           )
            CALL CHKOUT ( RNAME                                  )
            RETURN

         END IF

C
C        There should be no subtype specification in the method 
C        string.
C
         IF ( SUBTYP .NE. ' ' ) THEN
            
            CALL SETMSG ( 'Spurious sub-observer point type <#> '
     .      //            'was present in the method string #. ' 
     .      //            'The sub-observer type is valid in '
     .      //            'the method strings for SUBPNT and '
     .      //            'SUBSLR, but is not applicable for LATSRF.' )
            CALL ERRCH  ( '#', SUBTYP                                 )
            CALL ERRCH  ( '#', METHOD                                 )
            CALL SIGERR ( 'SPICE(INVALIDMETHOD)'                      )
            CALL CHKOUT ( RNAME                                       )
            RETURN

         END IF

         PRVMTH = METHOD

      END IF

C
C     At this point, the first pass actions were successful.
C
      FIRST = .FALSE.

C
C     Check the target body shape.
C
      IF ( SHAPE .EQ. ELLSHP ) THEN
 
         IF ( TRGCDE .NE. SVPRVT ) THEN
C
C           Reset counter to force lookup.
C           
            CALL ZZCTRUIN ( SVCTR3 )

         END IF
C
C        Look up target radii using counter.
C
         CALL ZZBODVCD ( TRGCDE, 'RADII', 3, SVCTR3, N, SVRADI )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         IF (  MIN( SVRADI(1), SVRADI(2), SVRADI(3) ) .LE. 0.D0 ) THEN

            CALL SETMSG ( 'Body # radii should be positive but '
     .      //            'were # # #.'                         )
            CALL ERRCH  ( '#',  TARGET                          )
            CALL ERRDP  ( '#',  SVRADI(1)                       )
            CALL ERRDP  ( '#',  SVRADI(2)                       )
            CALL ERRDP  ( '#',  SVRADI(3)                       )
            CALL SIGERR ( 'SPICE(BADAXISLENGTH)'                )
            CALL CHKOUT ( RNAME                                 )
            RETURN
            
         END IF

C
C        The radii are valid. Update the previous target ID.
C
         SVPRVT = TRGCDE

C
C        Generate surface points.
C
         DO I = 1, NPTS
C
C           Let X be a point having norm 1 and located at the Ith input
C           longitude and latitude.
C
            CALL LATREC ( 1.D0, LONLAT(1,I), LONLAT(2,I), X )
C
C           Scale X to place it on the ellipsoid surface. 
C
            CALL EDPNT ( X,         SVRADI(1), 
     .                   SVRADI(2), SVRADI(3), SRFPTS(1,I) )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF

         END DO



      ELSE IF  ( SHAPE .EQ. DSKSHP ) THEN
C
C        Initialize the DSK ray-surface intercept algorithm to use a
C        DSK model for the surface of the target body.
C
         CALL ZZSUDSKI ( TRGCDE, NSURF, SRFLST, FIXFID )

C
C        Get the radius of an outer bounding sphere for the body. Scale
C        up to avoid getting too close to the target surface.
C
         CALL ZZMAXRAD ( R )
         R = 2 * R
         
         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF


C
C        Generate surface points.
C
         DO I = 1, NPTS

            CALL LATREC ( R, LONLAT(1,I), LONLAT(2,I), X )

            CALL VMINUS ( X, RAYDIR )
C
C           Find the ray-surface intercept for the ray emanating
C           from X and pointing in the -X direction, where the
C           surface is represented by DSK data for the specified
C           body and surface list (the surface list was supplied
C           to ZZSUDSKI).
C            
            CALL ZZRAYSFX ( X, RAYDIR, ET, SRFPTS(1,I), FND )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF

            IF ( .NOT. FND ) THEN
               
               CALL SETMSG ( 'No surface point was found on body # '
     .         //            'at planetocentric longitude # (# deg), '
     .         //            'latitude # (# deg). This problem may '
     .         //            'be due to insufficient DSK data having '
     .         //            'been loaded for the body. It also '
     .         //            'could be due to the body having a '
     .         //            'shape not suitable for this computation '
     .         //            ', for example, a torus.'                 )
               CALL ERRCH  ( '#', TARGET                               )
               CALL ERRDP  ( '#', LONLAT(1,I)                          )
               CALL ERRDP  ( '#', LONLAT(1,I) * DPR()                  )
               CALL ERRDP  ( '#', LONLAT(2,I)                          )
               CALL ERRDP  ( '#', LONLAT(2,I) * DPR()                  )
               CALL SIGERR ( 'SPICE(POINTNOTFOUND)'                    )
               CALL CHKOUT ( RNAME                                     )
               RETURN
            
            END IF

C
C           Make sure the intercept is on the correct side of the
C           object.
C
            IF (  VDOT( X, SRFPTS(1,I) )  .LT. 0.D0  ) THEN

               CALL SETMSG ( 'A surface point was found on body # '
     .         //            'for the input planetocentric '
     .         //            'longitude # (# deg), '
     .         //            'latitude # (# deg), but this point '
     .         //            'is on the opposite side of the body. '
     .         //            'This likely indicates the '
     .         //            'the body does not contain the origin '
     .         //            'of the coordinate system. LATSRF does '
     .         //            'not work with such surfaces. Consider '
     .         //            'using DSKSXV for this computation.'    )
               CALL ERRCH  ( '#', TARGET                             )
               CALL ERRDP  ( '#', LONLAT(1,I)                        )
               CALL ERRDP  ( '#', LONLAT(1,I) * DPR()                )
               CALL ERRDP  ( '#', LONLAT(2,I)                        )
               CALL ERRDP  ( '#', LONLAT(2,I) * DPR()                )
               CALL SIGERR ( 'SPICE(SHAPENOTSUPPORTED)'              )
               CALL CHKOUT ( RNAME                                   )
               RETURN

            END IF

         END DO

         
      ELSE

         CALL SETMSG ( 'Input method <#> does not specify '
     .   //            'the target shape as either ELLIPSOID '
     .   //            'or DSK.'                               )
         CALL ERRCH  ( '#', METHOD                             )
         CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                   )
         CALL CHKOUT ( RNAME                                   )
         RETURN

      END IF
   
      CALL CHKOUT ( RNAME )
      RETURN
      END
