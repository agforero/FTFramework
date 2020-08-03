C$Procedure DSKXSI (DSK, ray-surface intercept with source information)
 
      SUBROUTINE DSKXSI ( PRI,    TARGET, NSURF,  SRFLST, ET,           
     .                    FIXREF, VERTEX, RAYDIR, MAXD,   MAXI,
     .                    XPT,    HANDLE, DLADSC, DSKDSC, DC,     
     .                    IC,     FOUND                        )

C$ Abstract
C
C     Compute a ray-surface intercept using data provided by 
C     multiple loaded DSK segments. Return information about 
C     the source of the data defining the surface on which the
C     intercept was found: DSK handle, DLA and DSK descriptors,
C     and DSK data type-dependent parameters.
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
C     CK
C     DSK
C     FRAMES
C     PCK
C     SPK
C     TIME
C
C$ Keywords
C
C     GEOMETRY
C     INTERCEPT
C     SURFACE
C     TOPOGRAPHY
C
C$ Declarations
 
      IMPLICIT NONE

      INCLUDE 'dsk.inc'
      INCLUDE 'dsktol.inc'    
      INCLUDE 'dsksrc.inc'
      INCLUDE 'srftrn.inc'
      INCLUDE 'zzctr.inc'
        
      LOGICAL               PRI
      CHARACTER*(*)         TARGET
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      DOUBLE PRECISION      ET
      CHARACTER*(*)         FIXREF
      DOUBLE PRECISION      VERTEX ( 3 )
      DOUBLE PRECISION      RAYDIR ( 3 )
      INTEGER               MAXD
      INTEGER               MAXI
      DOUBLE PRECISION      XPT    ( 3 )
      INTEGER               HANDLE
      INTEGER               DLADSC ( * )
      DOUBLE PRECISION      DSKDSC ( * )
      DOUBLE PRECISION      DC     ( * )
      INTEGER               IC     ( * )
      LOGICAL               FOUND
 
C$ Brief_I/O
C
C     VARIABLE  I/O  DESCRIPTION
C     --------  ---  --------------------------------------------------
C     PRI        I   Data prioritization flag. 
C     TARGET     I   Target body name.
C     NSURF      I   Number of surface IDs in list.
C     SRFLST     I   Surface ID list.
C     ET         I   Epoch, expressed as seconds past J2000 TDB.
C     FIXREF     I   Name of target body-fixed reference frame.
C     VERTEX     I   Vertex of ray.
C     RAYDIR     I   Direction vector of ray.
C     MAXD       I   Size of DC array.
C     MAXI       I   Size of IC array.
C     XPT        O   Intercept point.
C     HANDLE     O   Handle of segment contributing surface data.
C     DLADSC     O   DLA descriptor of segment.
C     DSKDSC     O   DSK descriptor of segment.
C     DC         O   Double precision component of source info.
C     IC         O   Integer component of source info.
C     FOUND      O   Found flag. 
C     DCSIZE     P   Required size of DC array.
C     ICSIZE     P   Required size of IC array.
C
C$ Detailed_Input
C
C     PRI        is a logical flag indicating whether to perform
C                a prioritized or unprioritized DSK segment search.
C                In an unprioritized search, no segment masks another:
C                data from all specified segments are used to  
C                define the surface of interest.
C
C                The search is unprioritized if and only if PRI 
C                is set to .FALSE. In the N0066 SPICE Toolkit, this
C                is the only allowed value.
C
C
C     TARGET     is the name of the target body on which a surface
C                intercept is sought.
C
C
C     NSURF,
C     SRFLST     are, respectively, a count of surface ID codes in a
C                list and an array containing the list. Only DSK
C                segments for the body designated by TARGET and having
C                surface IDs in this list will be considered in the
C                intercept computation. If the list is empty, all DSK
C                segments for TARGET will be considered.
C
C
C     ET         is the epoch of the intersection computation,
C                expressed as seconds past J2000 TDB. This epoch is
C                used only for DSK segment selection. Segments used in
C                the intercept computation must include ET in their
C                time coverage intervals.
C
C
C     FIXREF     is the name of a body-fixed, body-centered reference
C                frame associated with the target. The input ray vectors
C                are specified in this frame, as is the output intercept
C                point.
C      
C                The frame designated by FIXREF must have a fixed
C                orientation relative to the frame of any DSK segment
C                used in the computation.
C
C
C     VERTEX,
C     RAYDIR     are, respectively, the vertex and direction vector of
C                the ray to be used in the intercept computation. 
C
C                Both the vertex and ray's direction vector must be
C                represented in the reference frame designated by
C                FIXREF. The vertex is considered to be an offset from
C                the target body.
C     
C
C     MAXD,
C     MAXI       are, respectively, the declared sizes of the arrays
C                DC and IC. MAXD must be at least DCSIZE, while
C                MAXI must be at least ICSIZE. See the Parameters
C                section for details.
C
C$ Detailed_Output
C
C
C     XPT        is the intercept of the input ray on the surface 
C                specified by the inputs
C
C                   PRI
C                   TARGET
C                   NSURF
C                   SRFLST
C                   ET
C
C                if such an intercept exists. If the ray intersects the
C                surface at multiple points, the one closest to the
C                ray's vertex is selected.
C
C                XPT is defined if and only if FOUND is .TRUE.
C
C                Units are km.
C
C
C     HANDLE,
C     DLADSC,
C     DSKDSK     are, respectively, the DSK file handle, DLA descriptor,
C                and DSK descriptor of the DSK file and segment that 
C                contributed the surface data on which the intercept 
C                was found.
C
C                These outputs are defined if and only if FOUND is
C                .TRUE.
C
C     DC,
C     IC         are, respectively, double precision and integer arrays
C                that may contain additional information associated
C                with the segment contributing the surface data on
C                which the intercept was found. The information is 
C                DSK data type-dependent.
C
C                   For DSK type 2 segments
C
C                      IC(1) is the intercept plate ID.
C                      DC is unused.
C                  
C                These outputs are defined if and only if FOUND is
C                .TRUE.
C
C                The declared length of DC must be at least DSIZE;
C                the declared length of IC must be at least ISIZE.
C                See the Parameter section for details.
C
C
C     FOUND      is a logical flag that is set to .TRUE. if and only if
C                and intercept was found.
C
C
C$ Parameters
C
C     See the include file 
C
C        dsksrc.inc
C
C     for declarations of size parameters
C
C        DCSIZE
C        ICSIZE
C
C     for the output arguments
C
C        DC
C        IC
C
C     See the include files
C
C        dla.inc
C        dskdsc.inc
C
C     for declarations of DLA and DSK descriptor sizes and
C     documentation of the contents of these descriptors.
C
C     See the include file
C
C        dsktol.inc
C
C     for the values of tolerance parameters used by default by the
C     ray-surface intercept algorithm. These are discussed in in the
C     Particulars section below.
C
C$ Exceptions
C
C     1)  If the input prioritization flag PRI is set to .TRUE.,
C         the error SPICE(BADPRIORITYSPEC) is signaled.
C     
C     2)  If the input body name TARGET cannot be mapped to an
C         ID code, the error SPICE(IDCODENOTFOUND) is signaled.
C
C     3)  If the input frame name FIXREF cannot be mapped to an
C         ID code, the error SPICE(IDCODENOTFOUND) is signaled.
C
C     4)  If the frame center associated with FIXREF cannot be
C         retrieved, the error SPICE(NOFRAMEINFO) is signaled.
C
C     5)  If the frame center associated with FIXREF is not
C         the target body, the error SPICE(INVALIDFRAME) is signaled.
C
C     6)  If MAXD is less than DCSIZE or MAXI is less than ICSIZE,
C         the error SPICE(ARRAYTOOSMALL) will be signaled.
C
C     7)  If NSURF is less than 0, the error SPICE(INVALIDCOUNT)
C         is signaled.
C
C     8)  Any errors that occur during the intercept computation
C         will be signaled by routines in the call tree of this
C         routine.
C$ Files
C    
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - SPK data: ephemeris data for the positions of the centers
C          of DSK reference frames relative to the target body are
C          required if those frames are not centered at the target
C          body center.
C
C          Typically ephemeris data are made available by loading one
C          or more SPK files via FURNSH.
C
C        - DSK data: DSK files containing topographic data for the
C          target body must be loaded. If a surface list is specified,
C          data for at least one of the listed surfaces must be loaded.
C
C        - Frame data: if a frame definition is required to convert
C          DSK segment data to the body-fixed frame designated by
C          FIXREF, the target, that definition must be available in the
C          kernel pool. Typically the definitions of frames not already
C          built-in to SPICE are supplied by loading a frame kernel.
C
C        - CK data: if the frame to which FIXREF refers is a CK frame,
C          and if any DSK segments used in the computation have a
C          different frame, at least one CK file will be needed to
C          permit transformation of vectors between that frame and both
C          the J2000 and the target body-fixed frames.
C
C        - SCLK data: if a CK file is needed, an associated SCLK
C          kernel is required to enable conversion between encoded SCLK
C          (used to time-tag CK data) and barycentric dynamical time
C          (TDB).
C
C     In all cases, kernel data are normally loaded once per program
C     run, NOT every time this routine is called. 
C
C
C$ Particulars
C
C
C     This is the lowest-level public interface for computing
C     ray-surface intercepts, where the surface is modeled using
C     topographic data provided by DSK files. The highest-level
C     interface for this purpose is SINCPT.
C
C     In cases where the data source information returned by this
C     routine are not needed, the routine DSKXV may be more suitable.
C
C     This routine works with multiple DSK files. It places no
C     restrictions on the data types or coordinate systems of the DSK
C     segments used in the computation. DSK segments using different
C     reference frames may be used in a single computation. The only
C     restriction is that any pair of reference frames used directly or
C     indirectly are related by a constant rotation.
C
C     This routine enables calling applications to identify the source
C     of the data defining the surface on which an intercept was found.
C     The file, segment, and segment-specific information such as a DSK
C     type 2 plate ID are returned.
C
C     This routine can be used for improved efficiency in situations 
C     in which multiple ray-surface intercepts are to be performed
C     using a constant ray vertex.
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
C        in a later version of the routine, the presence of the PRI
C        argument is required.
C
C
C        Round-off errors and mitigating algorithms
C        ------------------------------------------
C
C        When topographic data are used to represent the surface of a
C        target body, round-off errors can produce some results that
C        may seem surprising.
C
C        Note that, since the surface in question might have mountains,
C        valleys, and cliffs, the points of intersection found for
C        nearly identical sets of inputs may be quite far apart from
C        each other: for example, a ray that hits a mountain side in a
C        nearly tangent fashion may, on a different host computer, be
C        found to miss the mountain and hit a valley floor much farther
C        from the observer, or even miss the target altogether.
C        
C        Round-off errors can affect segment selection: for example, a
C        ray that is expected to intersect the target body's surface
C        near the boundary between two segments might hit either
C        segment, or neither of them; the result may be
C        platform-dependent.
C
C        A similar situation exists when a surface is modeled by a set
C        of triangular plates, and the ray is expected to intersect the
C        surface near a plate boundary.
C        
C        To avoid having the routine fail to find an intersection when
C        one clearly should exist, this routine uses two "greedy"
C        algorithms:
C       
C           1) If the ray passes sufficiently close to any of the 
C              boundary surfaces of a segment (for example, surfaces of
C              maximum and minimum longitude or latitude), that segment
C              is tested for an intersection of the ray with the
C              surface represented by the segment's data.
C
C              This choice prevents all of the segments from being
C              missed when at least one should be hit, but it could, on
C              rare occasions, cause an intersection to be found in a
C              segment other than the one that would be found if higher
C              precision arithmetic were used.
C              
C           2) For type 2 segments, which represent surfaces as 
C              sets of triangular plates, each plate is expanded very
C              slightly before a ray-plate intersection test is
C              performed. The default plate expansion factor is 
C
C                 1 + XFRACT
C
C              where XFRACT is declared in 
C
C                 dsktol.inc
C
C              For example, given a value for XFRACT of 1.e-10, the
C              sides of the plate are lengthened by 1/10 of a micron
C              per km. The expansion keeps the centroid of the plate
C              fixed.
C
C              Plate expansion prevents all plates from being missed
C              in cases where clearly at least one should be hit.
C
C              As with the greedy segment selection algorithm, plate
C              expansion can occasionally cause an intercept to be
C              found on a different plate than would be found if higher
C              precision arithmetic were used. It also can occasionally
C              cause an intersection to be found when the ray misses
C              the target by a very small distance. 
C
C$ Examples
C
C     The numerical results shown for these examples may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine 
C     specific arithmetic implementation. 
C      
C     1)  Compute surface intercepts of rays emanating from a set of
C         vertices distributed on a longitude-latitude grid. All
C         vertices are outside the target body, and all rays point
C         toward the target's center.
C
C         Check intercepts against expected values. Indicate the
C         number of errors, the number of computations, and the
C         number of intercepts found.
C
C        Use the meta-kernel below to load example SPICE kernels. 
C
C
C           KPL/MK
C
C           File: dskxsi_ex1.tm
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
C              phobos512.bds                    DSK based on
C                                               Gaskell ICQ Q=512
C                                               plate model
C           \begindata
C
C              PATH_SYMBOLS    = 'GEN'
C              PATH_VALUES     = '/ftp/pub/naif/generic_kernels'
C
C              KERNELS_TO_LOAD = ( '$GEN/dsk/phobos/phobos512.bds' )
C
C           \begintext
C
C
C     Example code begins here.
C
C
C              PROGRAM SPEAR
C              IMPLICIT NONE
C        C
C        C     Multi-segment spear program.
C        C
C        C     This program expects all loaded DSKs
C        C     to represent the same body and surface.
C        C
C        C        Syntax: spear <meta-kernel>
C        C
C        C
C              INCLUDE 'dla.inc'
C              INCLUDE 'dsk.inc'
C              INCLUDE 'dskdsc.inc'
C              INCLUDE 'dsksrc.inc'
C              INCLUDE 'srftrn.inc'
C        C
C        C     SPICELIB functions
C        C
C              DOUBLE PRECISION      RPD
C              DOUBLE PRECISION      VDIST
C
C              LOGICAL               FAILED
C        C
C        C     Local parameters
C        C
C              DOUBLE PRECISION      DTOL
C              PARAMETER           ( DTOL   = 1.D-14 )
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
C              INTEGER               TYPLEN
C              PARAMETER           ( TYPLEN = 4 )
C
C        C
C        C     Local variables
C        C
C              CHARACTER*(FILSIZ)    DSK1
C              CHARACTER*(TYPLEN)    FILTYP
C              CHARACTER*(FRNMLN)    FIXREF
C              CHARACTER*(FILSIZ)    META
C              CHARACTER*(FILSIZ)    SOURCE
C              CHARACTER*(BDNMLN)    TARGET
C
C              DOUBLE PRECISION      D
C              DOUBLE PRECISION      DC     ( DCSIZE )
C              DOUBLE PRECISION      DSKDSC ( DSKDSZ )
C              DOUBLE PRECISION      ET
C              DOUBLE PRECISION      RAYDIR ( 3 )
C              DOUBLE PRECISION      LAT
C              DOUBLE PRECISION      LATCRD ( 3 )
C              DOUBLE PRECISION      LATSTP
C              DOUBLE PRECISION      LON
C              DOUBLE PRECISION      LONSTP
C              DOUBLE PRECISION      POLMRG
C              DOUBLE PRECISION      R
C              DOUBLE PRECISION      RADIUS
C              DOUBLE PRECISION      VERTEX ( 3 )
C              DOUBLE PRECISION      XPT    ( 3 )
C              DOUBLE PRECISION      XYZHIT ( 3 )
C
C              INTEGER               BODYID
C              INTEGER               DLADSC ( DLADSZ )
C              INTEGER               DTYPE
C              INTEGER               FRAMID
C              INTEGER               HANDLE
C              INTEGER               IC     ( ICSIZE )
C              INTEGER               NCASES
C              INTEGER               NDERR
C              INTEGER               NHITS
C              INTEGER               NLSTEP
C              INTEGER               NSURF
C              INTEGER               PLID
C              INTEGER               SRFLST ( MAXSRF )
C              INTEGER               SURFID
C
C              LOGICAL               FOUND
C
C
C              CALL CHKIN ( 'SPEAR' )
C        C
C        C     Load kernel.
C        C
C              CALL GETCML ( META )
C
C              IF ( META .EQ. ' ' ) THEN
C                 CALL TOSTDO( 'Syntax: spear <meta-kernel>' )
C                 CALL BYEBYE( 'SUCCESS' )
C              END IF
C        C
C        C     Load kernels.
C        C
C              CALL FURNSH ( META )
C        C
C        C     Get a handle for one of the loaded DSKs,
C        C     then find the first segment and extract
C        C     the body and surface IDs.
C        C
C              CALL KDATA  ( 1,      'DSK',  DSK1, FILTYP,
C             .              SOURCE, HANDLE, FOUND )
C
C              CALL DLABFS ( HANDLE, DLADSC, FOUND )
C
C              IF ( .NOT. FOUND ) THEN
C                 CALL SIGERR ( 'SPICE(NOSEGMENT)' )
C              END IF
C
C              CALL DSKGD ( HANDLE, DLADSC, DSKDSC )
C
C              BODYID = NINT( DSKDSC(CTRIDX) )
C              SURFID = NINT( DSKDSC(SRFIDX) )
C              FRAMID = NINT( DSKDSC(FRMIDX) )
C
C              CALL BODC2N ( BODYID, TARGET, FOUND )
C
C              IF ( .NOT. FOUND ) THEN
C                 CALL SETMSG ( 'Cannot map body ID # to a name.' )
C                 CALL ERRINT ( '#',  BODYID                      )
C                 CALL SIGERR ( 'SPICE(BODYNAMENOTFOUND)'         )
C              END IF
C
C              CALL FRMNAM ( FRAMID, FIXREF )
C
C              IF ( FIXREF .EQ. ' ' ) THEN
C                 CALL SETMSG ( 'Cannot map frame ID # to a name.' )
C                 CALL ERRINT ( '#',  FRAMID                       )
C                 CALL SIGERR ( 'SPICE(FRAMENAMENOTFOUND)'         )
C              END IF
C
C        C
C        C     Set the magnitude of the ray vertices. Use a large
C        C     number to ensure the vertices are outside of
C        C     any realistic target.
C        C
C              R = 1.D10
C
C        C
C        C     Spear the target with rays pointing toward
C        C     the origin.  Use a grid of ray vertices
C        C     located on a sphere enclosing the target.
C        C
C        C     The variable POLMRG ("pole margin") can
C        C     be set to a small positive value to reduce
C        C     the number of intercepts done at the poles.
C        C     This may speed up the computation for
C        C     the multi-segment case, since rays parallel
C        C     to the Z axis will cause all segments converging
C        C     at the pole of interest to be tested for an
C        C     intersection.
C        C
C              POLMRG = 5.D-1
C              LATSTP = 1.D0
C              LONSTP = 2.D0
C
C              NCASES = 0
C              NHITS  = 0
C              NDERR  = 0
C
C              LON    = -180.D0
C              LAT    = 90.D0
C              NLSTEP = 0
C
C
C              WRITE (*,*) 'Computing intercepts...'
C
C              DO WHILE ( LON .LT. 180.D0 )
C
C                 DO WHILE ( NLSTEP .LE. 180  )
C
C                    IF ( LON .EQ. -180.D0 ) THEN
C
C                       LAT = 90.D0 - NLSTEP*LATSTP
C
C                    ELSE
C
C                       IF ( NLSTEP .EQ. 0 ) THEN
C
C                          LAT =  90.D0 - POLMRG
C
C                       ELSE IF ( NLSTEP .EQ. 180 ) THEN
C
C                          LAT = -90.D0 + POLMRG
C
C                       ELSE
C
C                          LAT = 90.D0 - NLSTEP*LATSTP
C
C                       END IF
C
C                    END IF
C
C                    NCASES = NCASES + 1
C
C                    CALL LATREC ( R, LON*RPD(), LAT*RPD(), VERTEX )
C                    CALL VMINUS ( VERTEX, RAYDIR )
C
C                    NSURF     = 1
C                    SRFLST(1) = SURFID
C
C                    CALL DSKXSI ( .FALSE., TARGET, NSURF,  SRFLST,
C             .                    ET,      FIXREF, VERTEX, RAYDIR,
C             .                    DCSIZE,  ICSIZE, XPT,    HANDLE,
C             .                    DLADSC, DSKDSC,  DC,     IC,
C             .                    FOUND                           )
C
C                    IF ( .NOT. FAILED() .AND. FOUND ) THEN
C        C
C        C              Record that a new intercept was found.
C        C
C                       NHITS = NHITS + 1
C        C
C        C              Compute the latitude and longitude of
C        C              the intercept. Make sure these agree
C        C              well with those of the vertex.
C        C
C                       CALL RECLAT ( XPT,       LATCRD(1),
C             .                       LATCRD(2), LATCRD(3) )
C
C                       RADIUS = LATCRD(1)
C
C                       CALL LATREC ( RADIUS,   LON*RPD(),
C             .                       LAT*RPD(), XYZHIT    )
C
C                       D = VDIST ( XPT, XYZHIT )
C
C                       IF ( D/R .GT. DTOL ) THEN
C        C
C        C                 Get the intercept segment's plate ID if
C        C                 applicable.
C        C
C                          DTYPE = NINT( DSKDSC(TYPIDX) )
C
C                          WRITE (*,*) '======================'
C                          WRITE (*,*) 'LON, LAT       = ', LON, LAT
C                          WRITE (*,*) 'Bad intercept'
C                          WRITE (*,*) 'Distance error = ', D
C                          WRITE (*,*) 'XPT            = ', XPT
C                          WRITE (*,*) 'XYZHIT         = ', XYZHIT
C
C                          IF ( DTYPE .EQ. 2 ) THEN
C                             PLID = IC(1)
C                             WRITE (*,*) 'Plate ID      = ', PLID
C                          END IF
C
C                          NDERR = NDERR + 1
C
C                       END IF
C
C                    ELSE
C        C
C        C              Missing the target entirely is a fatal error.
C        C
C                       WRITE (*,*) '======================'
C                       WRITE (*,*) 'LON, LAT = ', LON, LAT
C                       WRITE (*,*) 'No intercept'
C                       WRITE (*,*) 'NCASES = ', NCASES
C                       STOP
C
C                    END IF
C
C                    NLSTEP = NLSTEP + 1
C
C                 END DO
C
C                 LON    = LON + LONSTP
C                 LAT    = 90.D0
C                 NLSTEP = 0
C
C              END DO
C
C              WRITE (*,*) 'Done.'
C              WRITE (*,*) ' '
C              WRITE (*,*) 'NCASES = ', NCASES
C              WRITE (*,*) 'NHITS  = ', NHITS
C              WRITE (*,*) 'NDERR  = ', NDERR
C              WRITE (*,*) ' '
C              END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was:
C
C
C        Computing intercepts...
C        Done.
C
C        NCASES =        32580
C        NHITS  =        32580
C        NDERR  =            0
C
C
C$ Restrictions
C
C     1)  The frame designated by FIXREF must have a fixed
C         orientation relative to the frame of any DSK segment
C         used in the computation. This routine has no 
C         practical way of ensuring that this condition is met;
C         so this responsibility is delegated to the calling
C         application.
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
C-    SPICELIB Version 1.0.0, 04-APR-2017 (NJB) 
C
C        Original 26-FEB-2016 (NJB) 
C
C
C-&
 
C$ Index_Entries
C
C     dsk ray-surface intercept with source information
C     dsk ray-surface intercept with handle and descriptors
C     dsk ray-surface intercept with handle and descriptors
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
      INTEGER               BDNMLN
      PARAMETER           ( BDNMLN = 36 )

      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )

C
C     Local variables
C
      CHARACTER*(FRNMLN)    PRVFRM
      CHARACTER*(BDNMLN)    SVTNAM

      INTEGER               FIXFID
      INTEGER               FRMCTR ( CTRSIZ )
      INTEGER               FXCENT
      INTEGER               FXCLSS
      INTEGER               FXTPID
      INTEGER               PRVTCD
      INTEGER               SVTCDE
      INTEGER               TRGCDE
      INTEGER               TRGCTR ( CTRSIZ )

      LOGICAL               FIRST
      LOGICAL               FRMFND
      LOGICAL               NEWFRM
      LOGICAL               NEWTRG
      LOGICAL               SVTFND
      LOGICAL               TRGFND
      LOGICAL               UPDATE

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  FIXFID
      SAVE                  FRMCTR
      SAVE                  PRVFRM
      SAVE                  PRVTCD
      SAVE                  SVTCDE
      SAVE                  SVTFND
      SAVE                  SVTNAM
      SAVE                  TRGCDE
      SAVE                  TRGCTR
C
C     Initial values
C
      DATA                  FIRST  / .TRUE. /
      DATA                  PRVFRM / ' '    /
      DATA                  PRVTCD / 0      /
      

      IF ( RETURN() ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'DSKXSI' )


      IF ( FIRST ) THEN
C
C        Initialize counters.
C
         CALL ZZCTRUIN( TRGCTR )
         CALL ZZCTRUIN( FRMCTR )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'DSKXSI' )
            RETURN
         END IF

      END IF

C
C     Reject PRI if not set properly.
C
      IF ( PRI ) THEN

         CALL SETMSG ( 'In the N0066 SPICE Toolkit, PRI must be '
     .   //            'set to .FALSE., indicating that an '
     .   //            'unprioritized search is to be performed.' )
         CALL SIGERR ( 'SPICE(BADPRIORITYSPEC)'                   )
         CALL CHKOUT ( 'DSKXSI'                                   )
         RETURN

      END IF

C
C     Reject NSURF if not set properly. Zero is a valid value.
C
      IF ( NSURF .LT. 0 ) THEN

         CALL SETMSG ( 'The surface count NSURF must be non-negative '
     .   //            'but was #.'                                   )
         CALL ERRINT ( '#',  NSURF                                    )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                          )
         CALL CHKOUT ( 'DSKXSI'                                       )
         RETURN

      END IF

C
C     Check output array sizes.
C
      IF (  ( MAXD .LT. DCSIZE ) .OR. ( MAXI .LT. ICSIZE )  ) THEN

         CALL SETMSG ( 'Output array size MAXD must be at least #; '
     .   //            'output array size MAXI must be at least #. '
     .   //            'Actual sizes were # and # respectively.'    )
         CALL ERRINT ( '#',  DCSIZE                                 )
         CALL ERRINT ( '#',  ICSIZE                                 )
         CALL ERRINT ( '#',  MAXD                                   )
         CALL ERRINT ( '#',  MAXI                                   )
         CALL SIGERR ( 'SPICE(ARRAYTOOSMALL)'                       )
         CALL CHKOUT ( 'DSKXSI'                                     )
         RETURN

      END IF

C
C     Obtain integer codes for the target and reference frame.
C 
      CALL ZZBODS2C ( TRGCTR, SVTNAM, SVTCDE, SVTFND,
     .                TARGET, TRGCDE, TRGFND         )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'DSKXSI' )
         RETURN
      END IF
            
      IF ( .NOT. TRGFND ) THEN
      
         CALL SETMSG ( 'The target, '
     .   //            '''#'', is not a recognized name for an '
     .   //            'ephemeris object. The cause of this '
     .   //            'problem may be that you need an updated '
     .   //            'version of the SPICE Toolkit, or that you '
     .   //            'failed to load a kernel containing a '
     .   //            'name-ID mapping for this body.'           )
         CALL ERRCH  ( '#', TARGET                                )
         CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                    )
         CALL CHKOUT ( 'DSKXSI'                                   )
         RETURN
      
      END IF


      NEWFRM = ( FIXREF .NE. PRVFRM ) .OR. FIRST
      NEWTRG = ( TRGCDE .NE. PRVTCD ) .OR. FIRST

C
C     Get the frame ID if the pool state has changed. The
C     first call to ZZPCKTRCK will indicate an update.
C
      CALL ZZPCTRCK ( FRMCTR, UPDATE )

      IF ( UPDATE .OR. NEWFRM .OR. NEWTRG ) THEN

         CALL NAMFRM ( FIXREF, FIXFID )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'DSKXSI' )
            RETURN
         END IF

         IF ( FIXFID .EQ. 0 ) THEN

            CALL SETMSG ( 'Reference frame # is not recognized by ' 
     .      //            'the SPICE frame subsystem. Possibly '    
     .      //            'a required frame definition kernel has ' 
     .      //            'not been loaded.'                        )
            CALL ERRCH  ( '#',  FIXREF                              )
            CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                   )
            CALL CHKOUT ( 'DSKXSI'                                  )
            RETURN

         END IF

C
C        Determine the attributes of the frame designated by FIXREF.
C
         CALL FRINFO ( FIXFID, FXCENT, FXCLSS, FXTPID, FRMFND )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'DSKXSI' )
            RETURN
         END IF

         IF ( .NOT. FRMFND ) THEN

            CALL SETMSG ( 'Attributes for reference frame # '
     .      //            'could not be obtained from ' 
     .      //            'the SPICE frame subsystem. Possibly '    
     .      //            'a required frame definition kernel has ' 
     .      //            'not been loaded.'                        )
            CALL ERRCH  ( '#',  FIXREF                              )
            CALL SIGERR ( 'SPICE(NOFRAMEINFO)'                      )
            CALL CHKOUT ( 'DSKXSI'                                  )
            RETURN

         END IF
C
C        Make sure that FIXREF is centered at the target body's center.
C
         IF ( FXCENT .NE. TRGCDE ) THEN

            CALL SETMSG ( 'Reference frame # is not centered at the ' 
     .      //            'target body #. The ID code of the frame '
     .      //            'center is #.'                             )
            CALL ERRCH  ( '#',  FIXREF                               )
            CALL ERRCH  ( '#',  TARGET                               )
            CALL ERRINT ( '#',  FXCENT                               )
            CALL SIGERR ( 'SPICE(INVALIDFRAME)'                      )
            CALL CHKOUT ( 'DSKXSI'                                   )
            RETURN

         END IF

C
C        We have a valid frame at this point. Save the name.
C
         FIRST  = .FALSE.
         PRVFRM = FIXREF

C
C        Update the previous target ID code as well.
C
         PRVTCD = TRGCDE

      END IF

C
C     TRGCDE and FIXFID are set.
C
C
C     Perform the intercept computation.
C
      CALL ZZSBFXRI ( TRGCDE, NSURF,  SRFLST, ET,
     .                FIXFID, VERTEX, RAYDIR, XPT, HANDLE,
     .                DLADSC, DSKDSC, DC,     IC,  FOUND   )


      CALL CHKOUT ( 'DSKXSI' )
      RETURN
      END



