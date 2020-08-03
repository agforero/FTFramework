C$Procedure DSKXV ( DSK, ray-surface intercept, vectorized )
 
      SUBROUTINE DSKXV ( PRI,    TARGET, NSURF,    
     .                   SRFLST, ET,     FIXREF, NRAYS,
     .                   VTXARR, DIRARR, XPTARR, FNDARR )

C$ Abstract
C
C     Compute ray-surface intercepts for a set of rays, using data
C     provided by multiple loaded DSK segments.
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

      INCLUDE 'dsktol.inc'    
      INCLUDE 'zzctr.inc'
        
      LOGICAL               PRI
      CHARACTER*(*)         TARGET
      INTEGER               NSURF
      INTEGER               SRFLST ( * )
      DOUBLE PRECISION      ET
      CHARACTER*(*)         FIXREF
      INTEGER               NRAYS
      DOUBLE PRECISION      VTXARR ( 3, * )
      DOUBLE PRECISION      DIRARR ( 3, * )
      DOUBLE PRECISION      XPTARR ( 3, * )
      LOGICAL               FNDARR ( * )
 
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
C     NRAYS      I   Number of rays.
C     VTXARR     I   Array of vertices of rays.
C     DIRARR     I   Array of direction vectors of rays.
C     XPTARR     O   Intercept point array.
C     FNDARR     O   Found flag array.
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
C                list and an containing the list. Only DSK segments for
C                for the body designated by TARGET and having surface
C                IDs in this list will considered in the intercept
C                computation. If the list is empty, all DSK segments
C                for TARGET will be considered.
C
C
C     ET         is the epoch of the intersection computation,
C                expressed as seconds past J2000 TDB. This epoch is
C                used only for DSK segment selection. Segments used
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
C     NRAYS,
C     VTXARR, 
C     DIRARR     are, respectively, a count of rays, an array containing
C                the vertices of rays, and an array containing the
C                direction vectors of the rays.
C
C                The ray's vertices are considered to represent offsets
C                from the center of the target body.
C
C                The rays' vertices and direction vectors are
C                represented in the reference frame designated by
C                FIXREF.
C     
C$ Detailed_Output
C
C
C     XPTARR     is an array containing the intercepts of the input
C                rays on the surface specified by the inputs
C
C                   PRI
C                   TARGET
C                   NSURF
C                   SRFLST
C                   ET
C
C                The Ith element of XPTARR is the intercept
C                corresponding to the Ith ray, if such an intercept
C                exists. If a ray intersects the surface at multiple
C                points, the intercept closest to the ray's vertex is
C                selected.
C
C                The Ith element of XPTARR is defined if and only if the
C                Ith element of FNDARR is .TRUE.
C
C                Units are km.
C
C
C     FNDARR     is an array of logical flags indicating whether the
C                input rays intersect the surface. The Ith element of
C                FNDARR is set to .TRUE. if and only if an intercept
C                was found for the Ith ray.
C
C
C$ Parameters
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
C     6)  If NRAYS is less than 1, the error SPICE(INVALIDCOUNT)
C         is signaled.
C
C     7)  If NSURF is less than 0, the error SPICE(INVALIDCOUNT)
C         is signaled.
C
C     8)  Any errors that occur during the intercept computation
C         will be signaled by routines in the call tree of this
C         routine.
C
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
C     This routine is suitable for efficient ray-surface intercept
C     computations in which the relative observer-target geometry is
C     constant but the rays vary.
C
C     For cases in which it is necessary to know the source of the
C     data defining the surface on which an intercept was found,
C     use the SPICELIB routine DSKXSI. 
C
C     For cases in which a ray's vertex is not explicitly known but is
C     defined by relative observer-target geometry, the SPICELIB
C     ray-surface intercept routine SINCPT should be used.
C
C     This routine works with multiple DSK files. It places no
C     restrictions on the data types or coordinate systems of the DSK
C     segments used in the computation. DSK segments using different
C     reference frames may be used in a single computation. The only
C     restriction is that any pair of reference frames used directly or
C     indirectly are related by a constant rotation.
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
C     1) Compute surface intercepts of rays emanating from a set of
C        vertices distributed on a longitude-latitude grid. All
C        vertices are outside the target body, and all rays point
C        toward the target's center.
C
C        Check intercepts against expected values. Indicate the
C        number of errors, the number of computations, and the
C        number of intercepts found.
C
C        Use the meta-kernel below to load example SPICE kernels. 
C
C
C           KPL/MK
C
C           File: dskxv_ex1.tm
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
C              PROGRAM VSPEAR
C              IMPLICIT NONE
C        C
C        C     Multi-segment, vectorized spear program.
C        C
C        C     This program expects all loaded DSKs
C        C     to represent the same body and surface.
C        C
C        C        Syntax: vspear <meta-kernel>
C        C
C              INCLUDE 'dla.inc'
C              INCLUDE 'dsk.inc'
C              INCLUDE 'dskdsc.inc'
C        C
C        C     SPICELIB functions
C        C
C              DOUBLE PRECISION      RPD
C              DOUBLE PRECISION      VDIST
C        C
C        C     Local parameters
C        C
C              DOUBLE PRECISION      DTOL
C              PARAMETER           ( DTOL   = 1.D-14 )
C
C              INTEGER               BDNMLN
C              PARAMETER           ( BDNMLN = 36 )
C
C              INTEGER               CMDLEN
C              PARAMETER           ( CMDLEN = 1000 )
C
C              INTEGER               FILSIZ
C              PARAMETER           ( FILSIZ = 255 )
C
C              INTEGER               LNSIZE
C              PARAMETER           ( LNSIZE = 80 )
C
C              INTEGER               FRNMLN
C              PARAMETER           ( FRNMLN = 32 )
C
C              INTEGER               MAXN
C              PARAMETER           ( MAXN   = 100000 )
C
C              INTEGER               TYPLEN
C              PARAMETER           ( TYPLEN = 4 )
C
C        C
C        C     Local variables
C        C
C              CHARACTER*(CMDLEN)    CMD
C              CHARACTER*(FILSIZ)    DSK1
C              CHARACTER*(TYPLEN)    FILTYP
C              CHARACTER*(FRNMLN)    FIXREF
C              CHARACTER*(FILSIZ)    FNAME
C              CHARACTER*(LNSIZE)    IDCH
C              CHARACTER*(FILSIZ)    SOURCE
C              CHARACTER*(BDNMLN)    TARGET
C
C              DOUBLE PRECISION      D
C              DOUBLE PRECISION      DSKDSC ( DSKDSZ )
C              DOUBLE PRECISION      DIRARR ( 3, MAXN )
C              DOUBLE PRECISION      ET
C              DOUBLE PRECISION      LAT
C              DOUBLE PRECISION      LATCRD ( 3 )
C              DOUBLE PRECISION      LATSTP
C              DOUBLE PRECISION      LON
C              DOUBLE PRECISION      LONSTP
C              DOUBLE PRECISION      POLMRG
C              DOUBLE PRECISION      R
C              DOUBLE PRECISION      RADIUS
C              DOUBLE PRECISION      VLAT
C              DOUBLE PRECISION      VLON
C              DOUBLE PRECISION      VRAD
C              DOUBLE PRECISION      VTXARR ( 3, MAXN )
C              DOUBLE PRECISION      XPTARR ( 3, MAXN )
C              DOUBLE PRECISION      XYZHIT ( 3 )
C
C              INTEGER               BODYID
C              INTEGER               DLADSC ( DLADSZ )
C              INTEGER               FRAMID
C              INTEGER               HANDLE
C              INTEGER               I
C              INTEGER               J
C              INTEGER               NDERR
C              INTEGER               NHITS
C              INTEGER               NLSTEP
C              INTEGER               NRAYS
C              INTEGER               NSURF
C              INTEGER               SRFLST ( MAXSRF )
C              INTEGER               SURFID
C
C              LOGICAL               FNDARR ( MAXN )
C              LOGICAL               FOUND
C
C        C
C        C     Saved variables
C        C
C        C     Save large arrays to avoid stack problems.
C        C
C              SAVE                  DIRARR
C              SAVE                  FNDARR
C              SAVE                  XPTARR
C
C
C              CALL CHKIN ( 'SPEAR' )
C        C
C        C     Load kernel.
C        C
C              CALL GETCML ( CMD )
C
C              IF ( CMD .EQ. ' ' ) THEN
C                 CALL TOSTDO( 'Syntax: spear <meta-kernel>' )
C                 CALL BYEBYE( 'SUCCESS' )
C              END IF
C
C        C
C        C     Pick the meta-kernel name from the command.
C        C
C              CALL NEXTWD ( CMD, FNAME, IDCH )
C        C
C        C     Load DSKs.
C        C
C              CALL FURNSH ( FNAME )
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
C        C     number to to ensure the vertices are outside of
C        C     any realistic target.
C        C
C              R = 1.D10
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
C              NHITS  = 0
C              NDERR  = 0
C
C              LON    = -180.D0
C              LAT    = 90.D0
C              NLSTEP = 0
C              NRAYS  = 0
C
C        C
C        C     Generate rays.
C        C
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
C                          LAT =  90.D0 - NLSTEP*LATSTP
C
C                       END IF
C
C                    END IF
C
C                    NRAYS  = NRAYS  + 1
C
C                    CALL LATREC ( R,               LON*RPD(),
C             .                    LAT*RPD(),       VTXARR(1,NRAYS) )
C                    CALL VMINUS ( VTXARR(1,NRAYS), DIRARR(1,NRAYS) )
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
C        C
C        C     Assign surface ID list.
C        C
C        C     Note that, if we knew that all files had the desired
C        C     surface ID, we could set `nsurf' to 0 and omit the
C        C     initialization of the surface ID list.
C        C
C              NSURF     = 1
C              SRFLST(1) = SURFID
C
C
C              WRITE (*,*) 'Computing intercepts...'
C
C              CALL DSKXV ( .FALSE., TARGET, NSURF, SRFLST,
C             .             ET,      FIXREF, NRAYS, VTXARR,
C             .             DIRARR,  XPTARR, FNDARR        )
C
C              WRITE (*,*) 'Done.'
C        C
C        C     Check results.
C        C
C              DO I = 1, NRAYS
C
C                 IF ( FNDARR(I) ) THEN
C        C
C        C           Record that a new intercept was found.
C        C
C                    NHITS = NHITS + 1
C        C
C        C           Compute the latitude and longitude of
C        C           the intercept. Make sure these agree
C        C           well with those of the vertex.
C        C
C                    CALL RECLAT ( XPTARR(1,I), LATCRD(1),
C             .                    LATCRD(2),   LATCRD(3) )
C
C                    RADIUS = LATCRD(1)
C
C        C
C        C           Recover the vertex longitude and latitude.
C        C
C                    CALL RECLAT ( VTXARR(1,I), VRAD, VLON, VLAT )
C
C                    CALL LATREC ( RADIUS,  VLON,
C             .                    VLAT,    XYZHIT )
C
C                    D = VDIST ( XPTARR(1,I), XYZHIT )
C
C                    IF ( D/R .GT. DTOL ) THEN
C        C
C        C              Get the intercept segment's plate ID if
C        C              applicable.
C        C
C                       WRITE (*,*) '======================'
C                       WRITE (*,*) 'LON, LAT       = ', LON, LAT
C                       WRITE (*,*) 'Bad intercept'
C                       WRITE (*,*) 'Distance error = ', D
C                       WRITE (*,*) 'XPT            = ',
C             .                     ( XPTARR(J,I), J = 1, 3 )
C                       WRITE (*,*) 'XYZHIT         = ', XYZHIT
C
C                       NDERR = NDERR + 1
C
C                    END IF
C
C                 ELSE
C        C
C        C           Missing the target entirely is a fatal error.
C        C
C        C           This is true only for this program, not in
C        C           general. For example, if the target shape is
C        C           a torus, many rays would miss the target.
C        C
C                    WRITE (*,*) '======================'
C                    WRITE (*,*) 'LON, LAT = ', LON, LAT
C                    WRITE (*,*) 'No intercept'
C                    WRITE (*,*) 'I        = ', I
C                    STOP
C
C                 END IF
C
C              END DO
C
C              WRITE (*,*) 'NRAYS  = ', NRAYS
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
C        NRAYS  =        32580
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
C-    SPICELIB Version 1.0.0, 21-FEB-2017 (NJB) 
C
C        Original 25-FEB-2016 (NJB) 
C
C
C-&
 
C$ Index_Entries
C
C     vectorized ray-surface intercept 
C     vectorized ray-dsk intercept      
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
      INTEGER               I
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

      CALL CHKIN ( 'DSKXV' )


      IF ( FIRST ) THEN
C
C        Initialize counters.
C
         CALL ZZCTRUIN( TRGCTR )
         CALL ZZCTRUIN( FRMCTR )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'DSKXV' )
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
         CALL CHKOUT ( 'DSKXV'                                    )
         RETURN

      END IF

C
C     Reject NRAYS if not set properly.
C
      IF ( NRAYS .LT. 1 ) THEN

         CALL SETMSG ( 'The ray count NRAYS must be at least 1 but '
     .   //            'was #.'                                     )
         CALL ERRINT ( '#',  NRAYS                                  )
         CALL SIGERR ( 'SPICE(INVALIDCOUNT)'                        )
         CALL CHKOUT ( 'DSKXV'                                      )
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
         CALL CHKOUT ( 'DSKXV'                                        )
         RETURN

      END IF

C
C     Obtain integer codes for the target and reference frame.
C 
      CALL ZZBODS2C ( TRGCTR, SVTNAM, SVTCDE, SVTFND,
     .                TARGET, TRGCDE, TRGFND         )
            
      IF ( FAILED() ) THEN
         CALL CHKOUT ( 'DSKXV' )
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
         CALL CHKOUT ( 'DSKXV'                                   )
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
            CALL CHKOUT ( 'DSKXV' )
            RETURN
         END IF

         IF ( FIXFID .EQ. 0 ) THEN

            CALL SETMSG ( 'Reference frame # is not recognized by ' 
     .      //            'the SPICE frame subsystem. Possibly '    
     .      //            'a required frame definition kernel has ' 
     .      //            'not been loaded.'                        )
            CALL ERRCH  ( '#',  FIXREF                              )
            CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                   )
            CALL CHKOUT ( 'DSKXV'                                   )
            RETURN

         END IF

C
C        Determine the attributes of the frame designated by FIXREF.
C
         CALL FRINFO ( FIXFID, FXCENT, FXCLSS, FXTPID, FRMFND )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( 'DSKXV' )
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
            CALL CHKOUT ( 'DSKXV'                                   )
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
            CALL CHKOUT ( 'DSKXV'                                    )
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
C     Perform the intercept computations.
C
      DO I = 1, NRAYS

         CALL ZZSBFXR ( TRGCDE,      NSURF,       SRFLST,
     .                  ET,          FIXFID,      VTXARR(1,I),
     .                  DIRARR(1,I), XPTARR(1,I), FNDARR(I)   )

         IF ( FAILED() ) THEN
            CALL CHKOUT( 'DSKXV' )
            RETURN
         END IF

      END DO

      CALL CHKOUT ( 'DSKXV' )
      RETURN
      END



