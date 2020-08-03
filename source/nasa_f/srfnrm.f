C$Procedure SRFNRM ( Map surface points to outward normal vectors )

      SUBROUTINE SRFNRM ( METHOD, TARGET, ET,    FIXREF,
     .                    NPTS,   SRFPTS, NORMLS        )

C$ Abstract
C
C     Map array of surface points on a specified target body to
C     the corresponding unit length outward surface normal vectors.
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
      INCLUDE 'dsktol.inc'
      INCLUDE 'gf.inc'
      INCLUDE 'zzctr.inc'
      INCLUDE 'zzdsk.inc'
      
      CHARACTER*(*)         METHOD
      CHARACTER*(*)         TARGET
      DOUBLE PRECISION      ET
      CHARACTER*(*)         FIXREF
      INTEGER               NPTS
      DOUBLE PRECISION      SRFPTS ( 3, * )
      DOUBLE PRECISION      NORMLS ( 3, * )

C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     METHOD     I   Computation method.
C     TARGET     I   Name of target body.
C     ET         I   Epoch in TDB seconds past J2000 TDB.
C     FIXREF     I   Body-fixed, body-centered target body frame.
C     NPTS       I   Number of surface points in input array.
C     SRFPTS     I   Array of surface points.
C     NORMLS     O   Array of outward, unit length normal vectors.
C     PTMEMM     P   Default point-surface membership margin.     
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
C                       The normal vector computation uses a triaxial
C                       ellipsoid to model the surface of the target
C                       body. The ellipsoid's radii must be available
C                       in the kernel pool.
C
C
C                    'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                       The normal vector computation uses topographic
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
C                 on the target body. FIXREF may be any such frame
C                 supported by the SPICE system, including built-in
C                 frames (documented in the Frames Required Reading)
C                 and frames defined by a loaded frame kernel (FK). The
C                 string FIXREF is case-insensitive, and leading and
C                 trailing blanks in FIXREF are not significant.
C
C                 The input surface points in the array SRFPTS are
C                 expressed relative to this reference frame, as are
C                 the normal vectors computed by this routine.
C
C
C     NPTS        is the number of surface points in the array SRFPTS.
C
C
C     SRFPTS      is an array of target body surface points. Elements
C
C                    SRFPTS(1,I)
C                    SRFPTS(2,I)
C                    SRFPTS(3,I)
C
C                 are the Cartesian coordinates, expressed in the
C                 reference frame designated by FIXREF, of the Ith
C                 surface point in the array. Each surface point 
C                 represents an offset from the center of that frame.
C
C                 All surface points must actually be "on" the surface,
C                 that is, the distance of each point from the surface
C                 must be less than a small margin. See the Parameters
C                 section below for a description of this margin.
C
C$ Detailed_Output
C
C     NORMLS      is an array of unit length, outward normal vectors
C                 corresponding to the points in SRFPTS. Elements
C
C                    NORMLS(1,I)
C                    NORMLS(2,I)
C                    NORMLS(3,I)
C
C                 are the Cartesian coordinates, expressed in the
C                 reference frame designated by FIXREF, of the Ith
C                 normal vector in the array.
C
C$ Parameters
C
C     PTMEMM      is the default point-surface membership margin. This
C                 margin limits the distance an input point can be from
C                 a surface and still be considered to lie on that
C                 surface.
C
C                 The details of the application of PTMEMM are
C                 implementation-dependent. In the DSK case, roughly
C                 speaking, a point-surface distance limit within a DSK
C                 segment is set to
C
C                    PTMEMM * MAXR
C
C                 where MAXR is the radius of an outer bounding sphere
C                 for the segment.
C
C                 For shapes modeled as ellipsoids, the expression
C                 above is applied to the maximum radius of the
C                 ellipsoid.
C
C                 See the include file
C
C                    dsktol.inc
C
C                 for the declaration of PTMEMM.
C
C                 This margin can be overridden. See dsktol.inc
C                 and DSKSTL for details.
C
C$ Exceptions
C
C     1)  If the target body name specified in the input string cannot
C         be converted to an integer ID code, the error
C         SPICE(IDCODENOTFOUND) is signaled.
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
C         loaded into the kernel pool prior to calling SRFNRM, the
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
C     9)  If METHOD specifies that the target surface is represented by
C         DSK data, and data representing the portion of the surface
C         corresponding to the surface points provided in SRFPTS are
C         not available, an error will be signaled by a routine in the
C         call tree of this routine.
C
C    10)  If an input surface point is not within a small tolerance
C         of the specified surface, the error SPICE(POINTNOTONSURFACE)
C         is signaled. See the Parameters section for details.
C
C    11)  If the target shape is "ellipsoid" and not all radii of the
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
C     1) Compute outward normal vectors at surface points on a target
C        body, where the points correspond to a given planetocentric
C        longitude/latitude grid. Use both ellipsoid and DSK shape
C        models.
C
C        Use the meta-kernel shown below to load the required SPICE
C        kernels.
C
C
C           KPL/MK
C
C           File: srfnrm_ex1.tm
C
C           This meta-kernel is intended to support operation of SPICE
C           example programs. The kernels shown here should not be
C           assumed to contain adequate or correct versions of data
C           required by SPICE-based user applications.
C
C           In order for an application to use this meta-kernel, the
C           kernels referenced here must be present in the user's
C           current working directory.
C
C           The names and contents of the kernels referenced
C           by this meta-kernel are as follows:
C
C              File name                        Contents
C              ---------                        --------
C              pck00010.tpc                     Planet orientation and
C                                               radii
C              phobos512.bds                    DSK based on
C                                               Gaskell ICQ Q=512
C                                               plate model
C           \begindata
C
C              PATH_SYMBOLS    = 'GEN'
C              PATH_VALUES     = '/ftp/pub/naif/generic_kernels'
C
C              KERNELS_TO_LOAD = ( '$GEN/pck/pck00010.tpc',
C                                  '$GEN/dsk/phobos/phobos512.bds' )
C           \begintext
C
C
C        Example code begins here.
C 
C
C             PROGRAM EX1
C             IMPLICIT NONE
C       C
C       C     SPICELIB functions
C       C
C             DOUBLE PRECISION      DPR
C             DOUBLE PRECISION      RPD
C       C
C       C     Local parameters
C       C
C             CHARACTER*(*)         FMT1
C             PARAMETER           ( FMT1   = '(1X,A,F11.6)' )
C
C             CHARACTER*(*)         META
C             PARAMETER           ( META   = 'srfnrm_ex1.tm' )
C
C
C             INTEGER               BDNMLN
C             PARAMETER           ( BDNMLN = 36 )
C
C             INTEGER               FILSIZ
C             PARAMETER           ( FILSIZ = 255 )
C
C             INTEGER               FRNMLN
C             PARAMETER           ( FRNMLN = 32 )
C
C             INTEGER               LNSIZE
C             PARAMETER           ( LNSIZE = 79 )
C
C             INTEGER               MAXN
C             PARAMETER           ( MAXN   = 100000 )
C
C             INTEGER               MTHLEN
C             PARAMETER           ( MTHLEN = 80 )
C       C
C       C     Local variables
C       C
C             CHARACTER*(FRNMLN)    FIXREF
C             CHARACTER*(MTHLEN)    METHOD ( 2 )
C             CHARACTER*(LNSIZE)    OUTLIN
C             CHARACTER*(BDNMLN)    TARGET
C
C             DOUBLE PRECISION      DLAT
C             DOUBLE PRECISION      DLON
C             DOUBLE PRECISION      ET
C             DOUBLE PRECISION      GRID   ( 2, MAXN )
C             DOUBLE PRECISION      LAT
C             DOUBLE PRECISION      LAT0
C             DOUBLE PRECISION      LON
C             DOUBLE PRECISION      LON0
C             DOUBLE PRECISION      NORMLS ( 3, MAXN, 2 )
C             DOUBLE PRECISION      NRMLAT
C             DOUBLE PRECISION      NRMLON
C             DOUBLE PRECISION      NRMRAD
C             DOUBLE PRECISION      SRFPTS ( 3, MAXN, 2 )
C             DOUBLE PRECISION      XLAT
C             DOUBLE PRECISION      XLON
C             DOUBLE PRECISION      XR
C
C             INTEGER               I
C             INTEGER               J
C             INTEGER               N
C             INTEGER               NLAT
C             INTEGER               NLON
C
C       C
C       C     Saved variables
C       C
C             SAVE                  GRID
C             SAVE                  NORMLS
C             SAVE                  SRFPTS
C
C       C
C       C     Load kernels.
C       C
C             CALL FURNSH ( META )
C       C
C       C     Set target, reference frame, and epoch.
C       C
C             TARGET = 'PHOBOS'
C             FIXREF = 'IAU_PHOBOS'
C             ET     = 0.D0
C       C
C       C     Use both a reference ellipsoid and DSK data
C       C     to represent the surface.
C       C
C             METHOD(1) = 'ELLIPSOID'
C             METHOD(2) = 'DSK/UNPRIORITIZED'
C       C
C       C     Set the grid dimensions.
C       C
C             NLON   = 6
C             NLAT   = 3
C       C
C       C     Derive evenly spaced grid separations and starting
C       C     values in the longitude and latitude dimensions.
C       C     Units are degrees.
C       C
C             LAT0 = 90.D0
C             LON0 =  0.D0
C
C             DLAT = 180.D0 / (NLAT + 1)
C             DLON = 360.D0 /  NLON
C
C       C
C       C     Now generate the grid points.  We generate
C       C     points along latitude bands, working from
C       C     north to south.  The latitude range is selected
C       C     to range from +45 to -45 degrees.  Longitude
C       C     ranges from 0 to 300 degrees.  The increment
C       C     is 45 degrees for latitude and 60 degrees for
C       C     longitude.
C       C
C             N = 0
C
C             DO I = 1, NLAT
C
C                LAT = RPD() * ( LAT0 - I*DLAT )
C
C                DO J = 1, NLON
C
C                   N   = N + 1
C                   LON = RPD() * ( LON0 + (J-1)*DLON )
C
C                   GRID(1,N) = LON
C                   GRID(2,N) = LAT
C
C                END DO
C
C             END DO
C
C       C
C       C     Find the surface points corresponding to the grid points.
C       C
C       C
C       C     Compute outward normal vectors at the surface points,
C       C     using both surface representations.
C       C
C             DO I = 1, 2
C
C                CALL LATSRF ( METHOD(I), TARGET, ET,
C            .                 FIXREF,    N,      GRID,  SRFPTS(1,1,I) )
C
C                CALL SRFNRM ( METHOD(I),
C            .                 TARGET,    ET,            FIXREF,
C            .                 N,         SRFPTS(1,1,I), NORMLS(1,1,I) )
C             END DO
C
C       C
C       C     Print out the surface points in latitudinal
C       C     coordinates and compare the derived lon/lat values
C       C     to those of the input grid.
C       C
C             DO I = 1, N
C       C
C       C        Use RECRAD rather than RECLAT to produce
C       C        non-negative longitudes.
C       C
C                CALL RECRAD ( SRFPTS(1,I,1), XR, XLON, XLAT )
C
C                WRITE (*,*) ' '
C
C                OUTLIN = ' Surface point for grid point #:'
C                CALL REPMI  ( OUTLIN, '#', I, OUTLIN )
C                CALL TOSTDO ( OUTLIN )
C
C                WRITE (*,*)    '  Latitudinal Coordinates:'
C                WRITE (*,FMT1) '    Longitude           (deg): ',
C            .                  XLON*DPR()
C                WRITE (*,FMT1) '    Latitude            (deg): ',
C            .                  XLAT*DPR()
C                WRITE (*,FMT1) '    Ellipsoid Radius     (km): ',
C            .                  XR
C
C                CALL RECRAD ( SRFPTS(1,I,2), XR, XLON, XLAT )
C
C                WRITE (*,FMT1) '    DSK Radius           (km): ',
C            .                  XR
C       C
C       C        Convert the Ith normal vector to latitudinal
C       C        coordinates.
C       C
C                CALL RECRAD ( NORMLS(1,I,1), NRMRAD, NRMLON, NRMLAT )
C
C                WRITE (*,*)    '  Ellipsoid normal vector direction:'
C                WRITE (*,FMT1) '    Longitude           (deg): ',
C            .                       NRMLON*DPR()
C                WRITE (*,FMT1) '    Latitude            (deg): ',
C            .                      NRMLAT*DPR()
C
C                CALL RECRAD ( NORMLS(1,I,2), NRMRAD, NRMLON, NRMLAT )
C
C                WRITE (*,*)    '  DSK normal vector direction:'
C                WRITE (*,FMT1) '    Longitude           (deg): ',
C            .                  NRMLON*DPR()
C                WRITE (*,FMT1) '    Latitude            (deg): ',
C            .                  NRMLAT*DPR()
C
C             END DO
C             WRITE (*,*) ' '
C             END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit
C     platform, the output for the first 3 points 
C     (the rest of the output is not shown due to its large volume)
C     was: 
C
C
C        Surface point for grid point 1:
C          Latitudinal Coordinates:
C            Longitude           (deg):    0.000000
C            Latitude            (deg):   45.000000
C            Ellipsoid Radius     (km):   10.542977
C            DSK Radius           (km):   10.156402
C          Ellipsoid normal vector direction:
C            Longitude           (deg):    0.000000
C            Latitude            (deg):   63.895146
C          DSK normal vector direction:
C            Longitude           (deg):  341.337568
C            Latitude            (deg):   62.610726
C
C        Surface point for grid point 2:
C          Latitudinal Coordinates:
C            Longitude           (deg):   60.000000
C            Latitude            (deg):   45.000000
C            Ellipsoid Radius     (km):   10.172847
C            DSK Radius           (km):   10.131412
C          Ellipsoid normal vector direction:
C            Longitude           (deg):   66.059787
C            Latitude            (deg):   58.877649
C          DSK normal vector direction:
C            Longitude           (deg):   48.859884
C            Latitude            (deg):   56.924717
C
C        Surface point for grid point 3:
C          Latitudinal Coordinates:
C            Longitude           (deg):  120.000000
C            Latitude            (deg):   45.000000
C            Ellipsoid Radius     (km):   10.172847
C            DSK Radius           (km):   10.423766
C          Ellipsoid normal vector direction:
C            Longitude           (deg):  113.940213
C            Latitude            (deg):   58.877649
C          DSK normal vector direction:
C            Longitude           (deg):  118.553200
C            Latitude            (deg):   55.906774
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
C     N.J. Bachman    (JPL)
C
C$ Version
C
C-    SPICELIB Version 1.0.0, 22-FEB-2017 (NJB) 
C
C        Added FAILED call.
C
C        01-JUL-2016 (NJB) 
C
C-&

C$ Index_Entries
C
C     map Cartesian surface points to normal vectors
C     compute normal vectors on topographic surface
C     compute normal vectors on dsk surface
C
C-&


C
C     SPICELIB functions
C
      LOGICAL               EQSTR
      LOGICAL               FAILED
      LOGICAL               RETURN

C
C     Local parameters
C
      CHARACTER*(*)         RNAME
      PARAMETER           ( RNAME  = 'SRFNRM' )

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
      
      DOUBLE PRECISION      A2
      DOUBLE PRECISION      B2
      DOUBLE PRECISION      C2
      DOUBLE PRECISION      LEVEL
      DOUBLE PRECISION      LIMIT
      DOUBLE PRECISION      MAXR
      DOUBLE PRECISION      PTSRFM

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
      SAVE                  MAXR
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
            
      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF


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
     .      //            'SUBSLR, but is not applicable for SRFNRM.' )
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
C        Compute the point-surface distance limit.
C
         MAXR   = MAX ( SVRADI(1), SVRADI(2), SVRADI(3) )

         CALL DSKGTL ( KEYPTM, PTSRFM )

         LIMIT  = PTSRFM * MAXR

C
C        Generate normal vectors.
C
         DO I = 1, NPTS
            
            A2    = SVRADI(1)*SVRADI(1)
            B2    = SVRADI(2)*SVRADI(2)
            C2    = SVRADI(3)*SVRADI(3)

            LEVEL = (   ( SRFPTS(1,I)*SRFPTS(1,I) / A2 )
     .                + ( SRFPTS(2,I)*SRFPTS(2,I) / B2 )
     .                + ( SRFPTS(3,I)*SRFPTS(3,I) / C2 )  ) ** 0.5D0

C
C           The test below is a distance test if the target shape
C           is a sphere. For other ellipsoids, it's an approximation.
C
            IF ( ABS( LEVEL - 1.D0 ) .GE. LIMIT ) THEN

               CALL SETMSG ( 'Input point at index # is not on '
     .         //            'the target body surface. The '
     .         //            'level surface parameter (x/a)**2 '
     .         //            '+ (y/b)**2 + (z/c)**2 for this '
     .         //            'point is #.'                      )
               CALL ERRINT ( '#', I                             )
               CALL ERRDP  ( '#', LEVEL                         )
               CALL SIGERR ( 'SPICE(POINTNOTONSURFACE)'         )
               CALL CHKOUT ( RNAME                              )
               RETURN

            END IF

            CALL SURFNM ( SVRADI(1),   SVRADI(2),  SVRADI(3), 
     .                    SRFPTS(1,I), NORMLS(1,I)            )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF


         END DO


      ELSE IF  ( SHAPE .EQ. DSKSHP ) THEN
C
C        Generate normal vectors.
C
         DO I = 1, NPTS
C
C           Use the DSK API segment buffering system to efficiently
C           select relevant segments and compute normals.
C 
            CALL ZZSBFNRM ( TRGCDE, NSURF,       SRFLST,     ET,    
     .                      FIXFID, SRFPTS(1,I), NORMLS(1,I)     )
           
            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF
            
C
C           Make sure normals have unit length.
C           
            CALL VHATIP ( NORMLS(1,I) )

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
