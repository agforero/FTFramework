C$Procedure SINCPT ( Surface intercept )
 
      SUBROUTINE SINCPT ( METHOD,  TARGET,  ET,      FIXREF, 
     .                    ABCORR,  OBSRVR,  DREF,    DVEC,    
     .                    SPOINT,  TRGEPC,  SRFVEC,  FOUND  )
      
C$ Abstract
C
C     Given an observer and a direction vector defining a ray, compute
C     the surface intercept of the ray on a target body at a specified
C     epoch, optionally corrected for light time and stellar
C     aberration.
C
C     The surface of the target body may be represented by a triaxial
C     ellipsoid or by topographic data provided by DSK files.
C
C     This routine supersedes SRFXPT.
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
C     NAIF_IDS
C     PCK
C     SCLK
C     SPK
C     TIME
C
C$ Keywords
C
C     GEOMETRY
C
C$ Declarations

      IMPLICIT NONE

      INCLUDE               'dsk.inc'      
      INCLUDE               'frmtyp.inc'
      INCLUDE               'gf.inc'
      INCLUDE               'zzabcorr.inc'
      INCLUDE               'zzctr.inc'
      INCLUDE               'zzdsk.inc'

      
      CHARACTER*(*)         METHOD
      CHARACTER*(*)         TARGET
      DOUBLE PRECISION      ET
      CHARACTER*(*)         FIXREF
      CHARACTER*(*)         ABCORR
      CHARACTER*(*)         OBSRVR
      CHARACTER*(*)         DREF
      DOUBLE PRECISION      DVEC   ( 3 )
      DOUBLE PRECISION      SPOINT ( 3 )
      DOUBLE PRECISION      TRGEPC
      DOUBLE PRECISION      SRFVEC ( 3 )
      LOGICAL               FOUND
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     METHOD     I   Computation method.
C     TARGET     I   Name of target body.
C     ET         I   Epoch in ephemeris seconds past J2000 TDB.
C     FIXREF     I   Body-fixed, body-centered target body frame.
C     ABCORR     I   Aberration correction.
C     OBSRVR     I   Name of observing body.
C     DREF       I   Reference frame of ray's direction vector.
C     DVEC       I   Ray's direction vector.
C     SPOINT     O   Surface intercept point on the target body.
C     TRGEPC     O   Intercept epoch.
C     SRFVEC     O   Vector from observer to intercept point.
C     FOUND      O   Flag indicating whether intercept was found.
C
C$ Detailed_Input
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
C                       The intercept computation uses a triaxial
C                       ellipsoid to model the surface of the target
C                       body. The ellipsoid's radii must be available
C                       in the kernel pool.
C
C
C                    'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                       The intercept computation uses topographic data
C                       to model the surface of the target body. These
C                       data must be provided by loaded DSK files.
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
C                 supply a string containing the integer ID code 
C                 for the object. For example both 'MOON' and '301'
C                 are legitimate strings that indicate the Moon is the
C                 target body.
C
C                 When the target body's surface is represented by a
C                 tri-axial ellipsoid, this routine assumes that a
C                 kernel variable representing the ellipsoid's radii is
C                 present in the kernel pool. Normally the kernel
C                 variable would be defined by loading a PCK file.
C
C
C     ET          is the epoch of participation of the observer,
C                 expressed as ephemeris seconds past J2000 TDB: ET is
C                 the epoch at which the observer's state is computed.
C
C                 When aberration corrections are not used, ET is also
C                 the epoch at which the state and orientation of the
C                 target body are computed.
C
C                 When aberration corrections are used, the position
C                 and orientation of the target body are computed at
C                 ET-LT or ET+LT, where LT is the one-way light time
C                 between the intercept point and the observer, and the
C                 sign applied to LT depends on the selected
C                 correction. See the description of ABCORR below for
C                 details.
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
C                 The output intercept point SPOINT and the observer-to-
C                 intercept vector SRFVEC will be expressed relative to
C                 this reference frame.
C
C
C     ABCORR      indicates the aberration corrections to be applied
C                 when computing the target's position and orientation.
C
C                 For remote sensing applications, where the apparent
C                 surface intercept point seen by the observer is
C                 desired, normally the correction
C              
C                    'CN+S'
C     
C                 should be used. This and the other supported options
C                 are described below. ABCORR may be any of the 
C                 following:
C
C                    'NONE'     Apply no correction. Return the 
C                               geometric surface intercept point on the
C                               target body.
C
C                 Let LT represent the one-way light time between the
C                 observer and the surface intercept point (note: NOT
C                 between the observer and the target body's center).
C                 The following values of ABCORR apply to the
C                 "reception" case in which photons depart from the
C                 intercept point's location at the light-time
C                 corrected epoch ET-LT and *arrive* at the observer's
C                 location at ET:
C
C                    'LT'       Correct for one-way light time (also
C                               called "planetary aberration") using a
C                               Newtonian formulation. This correction
C                               yields the location of the surface
C                               intercept point at the moment it
C                               emitted photons arriving at the
C                               observer at ET.
C 
C                               The light time correction uses an
C                               iterative solution of the light time
C                               equation. The solution invoked by the
C                               'LT' option uses one iteration.
C
C                               Both the target position as seen by the
C                               observer, and rotation of the target
C                               body, are corrected for light time.
C
C                    'LT+S'     Correct for one-way light time and
C                               stellar aberration using a Newtonian
C                               formulation. This option modifies the
C                               surface intercept obtained with the
C                               'LT' option to account for the
C                               observer's velocity relative to the
C                               solar system barycenter. These
C                               computations yield the apparent surface
C                               intercept point.
C
C                    'CN'       Converged Newtonian light time
C                               correction. In solving the light time
C                               equation, the 'CN' correction iterates
C                               until the solution converges. Both the
C                               position and rotation of the target
C                               body are corrected for light time.
C
C                    'CN+S'     Converged Newtonian light time and
C                               stellar aberration corrections. This
C                               option produces a solution that is at
C                               least as accurate at that obtainable
C                               with the 'LT+S' option. Whether the
C                               'CN+S' solution is substantially more
C                               accurate depends on the geometry of the
C                               participating objects and on the
C                               accuracy of the input data. In all
C                               cases this routine will execute more
C                               slowly when a converged solution is
C                               computed.
C
C                               For reception-case applications
C                               involving intercepts near the target
C                               body limb, this option should be used.
C
C                 The following values of ABCORR apply to the
C                 "transmission" case in which photons *depart* from
C                 the observer's location at ET and arrive at the
C                 intercept point at the light-time corrected epoch
C                 ET+LT:
C
C                    'XLT'      "Transmission" case: correct for
C                               one-way light time using a Newtonian
C                               formulation. This correction yields the
C                               intercept location at the moment it
C                               receives photons emitted from the
C                               observer's location at ET. 
C
C                               The light time correction uses an
C                               iterative solution of the light time
C                               equation. The solution invoked by the
C                               'XLT' option uses one iteration.
C
C                               Both the target position as seen by the
C                               observer, and rotation of the target
C                               body, are corrected for light time.
C
C                    'XLT+S'    "Transmission" case: correct for
C                               one-way light time and stellar
C                               aberration using a Newtonian
C                               formulation  This option modifies the
C                               intercept obtained with the 'XLT'
C                               option to account for the observer's
C                               velocity relative to the solar system
C                               barycenter.
C
C                    'XCN'      Converged Newtonian light time
C                               correction. This is the same as XLT
C                               correction but with further iterations
C                               to a converged Newtonian light time
C                               solution. 
C
C                    'XCN+S'    "Transmission" case: converged
C                               Newtonian light time and stellar
C                               aberration corrections. This option
C                               produces a solution that is at least as
C                               accurate at that obtainable with the
C                               'XLT+S' option. Whether the 'XCN+S'
C                               solution is substantially more accurate
C                               depends on the geometry of the
C                               participating objects and on the
C                               accuracy of the input data. In all
C                               cases this routine will execute more
C                               slowly when a converged solution is
C                               computed.
C
C                               For transmission-case applications
C                               involving intercepts near the target
C                               body limb, this option should be used.
C
C                Case and embedded blanks are not significant in
C                ABCORR. For example, the string
C
C                   'Cn + s'
C
C                 is valid.
C
C
C     OBSRVR      is the name of the observing body. This is typically
C                 a spacecraft, the earth, or a surface point on the
C                 earth. OBSRVR is case-insensitive, and leading and
C                 trailing blanks in OBSRVR are not significant.
C                 Optionally, you may supply a string containing the
C                 integer ID code for the object. For example both
C                 'MOON' and '301' are legitimate strings that indicate
C                 the Moon is the observer.
C
C
C     DREF        is the name of the reference frame relative to which
C                 the ray's direction vector is expressed. This may be
C                 any frame supported by the SPICE system, including
C                 built-in frames (documented in the Frames Required
C                 Reading) and frames defined by a loaded frame kernel
C                 (FK). The string DREF is case-insensitive, and
C                 leading and trailing blanks in DREF are not
C                 significant.
C
C                 When DREF designates a non-inertial frame, the
C                 orientation of the frame is evaluated at an epoch
C                 dependent on the frame's center and, if the center is
C                 not the observer, on the selected aberration
C                 correction. See the description of the direction
C                 vector DVEC for details.
C
C
C     DVEC        Ray direction vector emanating from the observer. The
C                 intercept with the target body's surface of the ray
C                 defined by the observer and DVEC is sought.
C
C                 DVEC is specified relative to the reference frame
C                 designated by DREF.
C
C                 Non-inertial reference frames are treated as follows:
C                 if the center of the frame is at the observer's
C                 location, the frame is evaluated at ET. If the
C                 frame's center is located elsewhere, then letting
C                 LTCENT be the one-way light time between the observer
C                 and the central body associated with the frame, the
C                 orientation of the frame is evaluated at ET-LTCENT,
C                 ET+LTCENT, or ET depending on whether the requested
C                 aberration correction is, respectively, for received
C                 radiation, transmitted radiation, or is omitted.
C                 LTCENT is computed using the method indicated by
C                 ABCORR.
C                 
C
C$ Detailed_Output
C
C
C     SPOINT      is the surface intercept point on the target body of
C                 the ray defined by the observer and the direction
C                 vector. If the ray intersects the target body in
C                 multiple points, the selected intersection point is
C                 the one closest to the observer. The output argument
C                 FOUND (see below) indicates whether an intercept was
C                 found.
C
C                 SPOINT is expressed in Cartesian coordinates,
C                 relative to the target body-fixed frame designated by
C                 FIXREF. The body-fixed target frame is evaluated at
C                 the intercept epoch TRGEPC (see description below).
C
C                 When light time correction is used, the duration of
C                 light travel between SPOINT to the observer is
C                 considered to be the one way light time. When both
C                 light time and stellar aberration corrections are
C                 used, SPOINT is compute such that, when the vector
C                 from the observer to SPOINT is corrected for light
C                 time and stellar aberration, the resulting vector
C                 lies on the ray defined by the observer's location
C                 and DVEC.
C
C                 The components of SPOINT are given in units of km.
C
C
C     TRGEPC      is the "intercept epoch." TRGEPC is defined as
C                 follows: letting LT be the one-way light time between
C                 the observer and the intercept point, TRGEPC is the
C                 epoch ET-LT, ET+LT, or ET depending on whether the
C                 requested aberration correction is, respectively, for
C                 received radiation, transmitted radiation, or
C                 omitted. LT is computed using the method indicated by
C                 ABCORR.
C
C                 TRGEPC is expressed as seconds past J2000 TDB.
C
C
C     SRFVEC      is the vector from the observer's position at ET to
C                 the aberration-corrected (or optionally, geometric)
C                 position of SPOINT, where the aberration corrections
C                 are specified by ABCORR. SRFVEC is expressed in the
C                 target body-fixed reference frame designated by
C                 FIXREF, evaluated at TRGEPC.
C  
C                 The components of SRFVEC are given in units of km.
C
C                 One can use the SPICELIB function VNORM to obtain the
C                 distance between the observer and SPOINT:
C
C                    DIST = VNORM ( SRFVEC )
C
C                 The observer's position OBSPOS, relative to the
C                 target body's center, where the center's position is
C                 corrected for aberration effects as indicated by
C                 ABCORR, can be computed via the call:
C
C                    CALL VSUB ( SPOINT, SRFVEC, OBSPOS )
C
C                 To transform the vector SRFVEC from a reference frame
C                 FIXREF at time TRGEPC to a time-dependent reference
C                 frame REF at time ET, the routine PXFRM2 should be
C                 called. Let XFORM be the 3x3 matrix representing the
C                 rotation from the reference frame FIXREF at time
C                 TRGEPC to the reference frame REF at time ET. Then
C                 SRFVEC can be transformed to the result REFVEC as
C                 follows:
C
C                     CALL PXFRM2 ( FIXREF, REF,    TRGEPC, ET, XFORM )
C                     CALL MXV    ( XFORM,  SRFVEC, REFVEC )
C
C                 The second example in the Examples header section
C                 below presents a complete program that demonstrates
C                 this procedure.
C     
C
C     FOUND       A logical flag indicating whether or not the ray
C                 intersects the target. If an intersection exists
C                 FOUND will be returned as .TRUE. If the ray misses
C                 the target, FOUND will be returned as .FALSE.
C
C$ Parameters
C
C     None.
C
C$ Exceptions
C 
C
C     1)  If the specified aberration correction is unrecognized, the
C         error will be diagnosed and signaled by a routine in the call
C         tree of this routine.
C
C     2)  If either the target or observer input strings cannot be
C         converted to an integer ID code, the error
C         SPICE(IDCODENOTFOUND) is signaled.
C
C     3)  If OBSRVR and TARGET map to the same NAIF integer ID code,
C         the error SPICE(BODIESNOTDISTINCT) is signaled.
C
C     4)  If the input target body-fixed frame FIXREF is not
C         recognized, the error SPICE(NOFRAME) is signaled. A frame
C         name may fail to be recognized because a required frame
C         specification kernel has not been loaded; another cause is a
C         misspelling of the frame name.
C
C     5)  If the input frame FIXREF is not centered at the target body,
C         the error SPICE(INVALIDFRAME) is signaled.
C
C     6)  If the input argument METHOD cannot be parsed, the error
C         is signaled either by this routine or by a routine in the
C         call tree of this routine.
C
C     7)  If the target and observer have distinct identities but are
C         at the same location (for example, the target is Mars and the
C         observer is the Mars barycenter), the error
C         SPICE(NOSEPARATION) is signaled.
C
C     8)  If insufficient ephemeris data have been loaded prior to
C         calling SINCPT, the error will be diagnosed and signaled by a
C         routine in the call tree of this routine. Note that when
C         light time correction is used, sufficient ephemeris data must
C         be available to propagate the states of both observer and
C         target to the solar system barycenter.
C
C     9)  If the computation method specifies an ellipsoidal target
C         shape and triaxial radii of the target body have not been
C         loaded into the kernel pool prior to calling SINCPT, the
C         error will be diagnosed and signaled by a routine in the call
C         tree of this routine.
C
C     10) The target must be an extended body: if any of the radii of
C         the target body are non-positive, the error will be
C         diagnosed and signaled by routines in the call tree of this
C         routine.
C
C     11) If PCK data specifying the target body-fixed frame
C         orientation have not been loaded prior to calling SINCPT,
C         the error will be diagnosed and signaled by a routine in the
C         call tree of this routine.
C
C     12) If the reference frame designated by DREF is not recognized
C         by the SPICE frame subsystem, the error SPICE(NOFRAME)
C         will be signaled.
C
C     13) If the direction vector DVEC is the zero vector, the error
C         SPICE(ZEROVECTOR) will be signaled.
C
C     14) If METHOD specifies that the target surface is represented by
C         DSK data, and no DSK files are loaded for the specified
C         target, the error is signaled by a routine in the call tree
C         of this routine.
C         
C     15) If METHOD specifies that the target surface is represented
C         by DSK data, and DSK data are not available for a portion of 
C         the target body's surface, an intercept might not be found.
C         This routine does not revert to using an ellipsoidal surface
C         in this case.
C
C$ Files
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - SPK data: ephemeris data for target and observer must be
C          loaded. If aberration corrections are used, the states of
C          target and observer relative to the solar system barycenter
C          must be calculable from the available ephemeris data.
C          Typically ephemeris data are made available by loading one
C          or more SPK files via FURNSH.
C
C        - PCK data: if the computation method is specified as
C          "Ellipsoid," triaxial radii for the target body must be 
C          loaded into the kernel pool. Typically this is done by
C          loading a text PCK file via FURNSH.
C
C        - Further PCK data: rotation data for the target body must
C          be loaded. These may be provided in a text or binary PCK
C          file. 
C
C     The following data may be required:
C
C        - DSK data: if METHOD indicates that DSK data are to be used,
C          DSK files containing topographic data for the target body
C          must be loaded. If a surface list is specified, data for
C          at least one of the listed surfaces must be loaded.
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
C          kernel containing the necessary assignments. An example
C          of such an assignment is
C
C             NAIF_SURFACE_NAME += 'Mars MEGDR 128 PIXEL/DEG'
C             NAIF_SURFACE_CODE += 1                    
C             NAIF_SURFACE_BODY += 499
C
C        - Frame data: if a frame definition is required to convert
C          the observer and target states to the body-fixed frame of
C          the target, that definition must be available in the kernel
C          pool. Similarly, the frame definition required to map
C          between the frame designated by DREF and the target
C          body-fixed frame must be available. Typically the
C          definitions of frames not already built-in to SPICE are
C          supplied by loading a frame kernel.
C
C        - CK data: if the frame to which DREF refers is fixed to a
C          spacecraft instrument or structure, at least one CK file
C          will be needed to permit transformation of vectors between
C          that frame and both the J2000 and the target body-fixed
C          frames.
C
C        - SCLK data: if a CK file is needed, an associated SCLK
C          kernel is required to enable conversion between encoded SCLK
C          (used to time-tag CK data) and barycentric dynamical time
C          (TDB).
C
C     In all cases, kernel data are normally loaded once per program
C     run, NOT every time this routine is called.
C
C$ Particulars
C
C     Given a ray defined by a direction vector and the location of an
C     observer, SINCPT computes the surface intercept point of the ray
C     on a specified target body. SINCPT also determines the vector
C     from the observer to the surface intercept point. If the ray
C     intersects the target in multiple locations, the intercept
C     closest to the observer is selected.
C
C     When aberration corrections are used, this routine finds the
C     value of SPOINT such that, if SPOINT is regarded as an ephemeris
C     object, after the selected aberration corrections are applied to
C     the vector from the observer to SPOINT, the resulting vector is
C     parallel to the direction vector DVEC.
C
C     This routine computes light time corrections using light time
C     between the observer and the surface intercept point, as opposed
C     to the center of the target. Similarly, stellar aberration
C     corrections done by this routine are based on the direction of
C     the vector from the observer to the light-time corrected
C     intercept point, not to the target center. This technique avoids
C     errors due to the differential between aberration corrections
C     across the target body. Therefore it's valid to use aberration
C     corrections with this routine even when the observer is very
C     close to the intercept point, in particular when the
C     observer-intercept point distance is much less than the
C     observer-target center distance. It's also valid to use stellar
C     aberration corrections even when the intercept point is near or
C     on the limb (as may occur in occultation computations using a
C     point target).
C
C     When comparing surface intercept point computations with results
C     from sources other than SPICE, it's essential to make sure the
C     same geometric definitions are used.
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
C                 1 + 1.E-10
C
C              In other words, the sides of the plate are lengthened by
C              1/10 of a micron per km. The expansion keeps the centroid
C              of the plate fixed.
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
C         
C        Aberration corrections
C        ----------------------
C
C        For irregularly shaped target bodies, the distance between the
C        observer and the nearest surface intercept need not be a
C        continuous function of time; hence the one-way light time
C        between the intercept and the observer may be discontinuous as
C        well. In such cases, the computed light time, which is found
C        using an iterative algorithm, may converge slowly or not at
C        all. In all cases, the light time computation will terminate,
C        but the result may be less accurate than expected.
C     
C
C$ Examples
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as
C     input, the compiler and supporting libraries, and the machine 
C     specific arithmetic implementation. 
C
C
C     1) The following program computes surface intercept points on Mars
C        for the boresight and FOV boundary vectors of the MGS MOC
C        narrow angle camera. The intercepts are computed for a single
C        observation epoch. Light time and stellar aberration
C        corrections are used. For simplicity, camera distortion is
C        ignored.
C
C        Intercepts are computed using both triaxial ellipsoid and 
C        topographic surface models. 
C
C        The topographic model is based on data from the MGS MOLA DEM
C        megr90n000cb, which has a resolution of 4 pixels/degree. A
C        triangular plate model was produced by computing a 720 x 1440
C        grid of interpolated heights from this DEM, then tessellating
C        the height grid. The plate model is stored in a type 2 segment
C        in the referenced DSK file.
C
C        Use the meta-kernel shown below to load the required SPICE
C        kernels.
C           
C
C           KPL/MK
C
C           File: sincpt_ex1.tm
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
C              de430.bsp                        Planetary ephemeris
C              mar097.bsp                       Mars satellite ephemeris
C              pck00010.tpc                     Planet orientation and
C                                               radii
C              naif0011.tls                     Leapseconds 
C              mgs_moc_v20.ti                   MGS MOC instrument
C                                               parameters
C              mgs_sclkscet_00061.tsc           MGS SCLK coefficients
C              mgs_sc_ext12.bc                  MGS s/c bus attitude
C              mgs_ext12_ipng_mgs95j.bsp        MGS ephemeris
C              megr90n000cb_plate.bds           Plate model based on
C                                               MEGDR DEM, resolution
C                                               4 pixels/degree.
C
C           \begindata
C
C              KERNELS_TO_LOAD = ( 'de430.bsp',
C                                  'mar097.bsp',
C                                  'pck00010.tpc',
C                                  'naif0011.tls',
C                                  'mgs_moc_v20.ti',
C                                  'mgs_sclkscet_00061.tsc',
C                                  'mgs_sc_ext12.bc',
C                                  'mgs_ext12_ipng_mgs95j.bsp',
C                                  'megr90n000cb_plate.bds'      )
C           \begintext
C
C
C        Example code begins here.
C
C
C           PROGRAM EX1
C           IMPLICIT NONE
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      VNORM
C
C     C
C     C     Local parameters
C     C
C           CHARACTER*(*)         META
C           PARAMETER           ( META   = 'sincpt_ex1.tm' )
C
C           INTEGER               ABCLEN
C           PARAMETER           ( ABCLEN = 20 )
C
C           INTEGER               LNSIZE
C           PARAMETER           ( LNSIZE = 78 )
C
C           INTEGER               METLEN
C           PARAMETER           ( METLEN = 40 )
C
C           INTEGER               NAMLEN
C           PARAMETER           ( NAMLEN = 32 )
C
C           INTEGER               TIMLEN
C           PARAMETER           ( TIMLEN = 50 )
C
C           INTEGER               SHPLEN
C           PARAMETER           ( SHPLEN = 80 )
C
C           INTEGER               NCORNR
C           PARAMETER           ( NCORNR = 4 )
C
C           INTEGER               NMETH
C           PARAMETER           ( NMETH  = 2 )
C
C     C
C     C     Local variables
C     C
C           CHARACTER*(ABCLEN)    ABCORR
C           CHARACTER*(NAMLEN)    CAMERA
C           CHARACTER*(NAMLEN)    DREF
C           CHARACTER*(NAMLEN)    FIXREF
C           CHARACTER*(METLEN)    METHDS ( NMETH )
C           CHARACTER*(METLEN)    METHOD
C           CHARACTER*(NAMLEN)    OBSRVR
C           CHARACTER*(SHPLEN)    SHAPE
C           CHARACTER*(NAMLEN)    SRFTYP ( NMETH )
C           CHARACTER*(NAMLEN)    TARGET
C           CHARACTER*(LNSIZE)    TITLE
C           CHARACTER*(TIMLEN)    UTC
C
C           DOUBLE PRECISION      BOUNDS ( 3, NCORNR )
C           DOUBLE PRECISION      BSIGHT ( 3 )
C           DOUBLE PRECISION      DIST
C           DOUBLE PRECISION      DPR
C           DOUBLE PRECISION      DVEC   ( 3 )
C           DOUBLE PRECISION      ET
C           DOUBLE PRECISION      LAT
C           DOUBLE PRECISION      LON
C           DOUBLE PRECISION      RADIUS
C           DOUBLE PRECISION      SPOINT ( 3 )
C           DOUBLE PRECISION      SRFVEC ( 3 )
C           DOUBLE PRECISION      TRGEPC
C
C           INTEGER               CAMID
C           INTEGER               I
C           INTEGER               J
C           INTEGER               K
C           INTEGER               N
C
C           LOGICAL               FOUND
C
C           DATA                  ABCORR / 'CN+S'              /
C           DATA                  CAMERA / 'MGS_MOC_NA'        /
C           DATA                  FIXREF / 'IAU_MARS'          /
C           DATA                  METHDS / 'ELLIPSOID',
C          .                               'DSK/UNPRIORITIZED' /
C           DATA                  OBSRVR / 'MGS'               /
C           DATA                  SRFTYP / 'Ellipsoid',
C          .                      'MGS/MOLA topography, 4 pixel/deg'  /
C           DATA                  TARGET / 'Mars'                     /
C           DATA                  UTC    / '2003 OCT 13 06:00:00 UTC' /
C
C     C
C     C     Load kernel files:
C     C
C           CALL FURNSH ( META )
C
C     C
C     C     Convert the UTC request time to ET (seconds past
C     C     J2000, TDB).
C     C
C           CALL STR2ET ( UTC, ET )
C
C     C
C     C     Get the MGS MOC Narrow angle camera (MGS_MOC_NA)
C     C     ID code. Then look up the field of view (FOV)
C     C     parameters by calling GETFOV.
C     C
C           CALL BODN2C ( CAMERA, CAMID, FOUND )
C
C           IF ( .NOT. FOUND ) THEN
C              CALL SETMSG ( 'Could not find ID code for ' //
C          .                 'instrument #.'               )
C              CALL ERRCH  ( '#', CAMERA                   )
C              CALL SIGERR ( 'SPICE(NOTRANSLATION)'        )
C           END IF
C
C     C
C     C     GETFOV will return the name of the camera-fixed frame
C     C     in the string DREF, the camera boresight vector in
C     C     the array BSIGHT, and the FOV corner vectors in the
C     C     array BOUNDS.
C     C
C           CALL GETFOV ( CAMID,  NCORNR, SHAPE,  DREF,
C          .              BSIGHT, N,      BOUNDS       )
C
C
C           WRITE (*,*) ' '
C           WRITE (*,*) 'Surface Intercept Locations for Camera'
C           WRITE (*,*) 'FOV Boundary and Boresight Vectors'
C           WRITE (*,*) ' '
C           WRITE (*,*) '   Instrument:            ', CAMERA
C           WRITE (*,*) '   Epoch:                 ', UTC
C           WRITE (*,*) '   Aberration correction: ', ABCORR
C           WRITE (*,*) ' '
C
C     C
C     C     Now compute and display the surface intercepts for the
C     C     boresight and all of the FOV boundary vectors.
C     C
C           DO I = 1, NCORNR+1
C
C              IF ( I .LE. NCORNR ) THEN
C
C                 TITLE = 'Corner vector #'
C                 CALL REPMI ( TITLE, '#', I, TITLE )
C
C                 CALL VEQU ( BOUNDS(1,I), DVEC )
C
C              ELSE
C
C                 TITLE = 'Boresight vector'
C                 CALL VEQU ( BSIGHT, DVEC )
C
C              END IF
C
C              WRITE (*,*) ' '
C              WRITE (*,*) TITLE
C
C              TITLE = '  Vector in # frame = '
C              CALL REPMC ( TITLE, '#', DREF, TITLE )
C
C              WRITE (*,*) ' '
C              WRITE (*,*) TITLE
C
C              IF ( I .LE. NCORNR ) THEN
C                 WRITE (*, '(1X,3F20.14)') ( BOUNDS(J,I), J = 1, 3 )
C              ELSE
C                 WRITE (*, '(1X,3F20.14)') BSIGHT
C              END IF
C
C              WRITE (*,*) ' '
C              WRITE (*,*) '  Intercept:'
C
C     C
C     C        Compute the surface intercept point using
C     C        the specified aberration corrections. Loop
C     C        over the set of computation methods.
C     C
C              DO K = 1, NMETH
C
C                 METHOD = METHDS(K)
C
C                 CALL SINCPT ( METHOD, TARGET, ET,
C          .                    FIXREF, ABCORR, OBSRVR,
C          .                    DREF,   DVEC,   SPOINT,
C          .                    TRGEPC, SRFVEC, FOUND   )
C
C                 IF ( FOUND ) THEN
C     C
C     C              Compute range from observer to apparent
C     C              intercept.
C     C
C                    DIST = VNORM ( SRFVEC )
C     C
C     C              Convert rectangular coordinates to
C     C              planetocentric latitude and longitude.
C     C              Convert radians to degrees.
C     C
C                    CALL RECLAT ( SPOINT, RADIUS, LON, LAT )
C
C                    LON = LON * DPR ()
C                    LAT = LAT * DPR ()
C     C
C     C              Display the results.
C     C
C
C                    WRITE (*,*) ' '
C                    CALL TOSTDO ( '     Surface representation: '
C          .         //            SRFTYP(K)                      )
C                    WRITE (*,*) ' '
C                    WRITE (*,*)
C          .         '     Radius                   (km)  = ', RADIUS
C                    WRITE (*,*)
C          .         '     Planetocentric Latitude  (deg) = ', LAT
C                    WRITE (*,*)
C          .         '     Planetocentric Longitude (deg) = ', LON
C                    WRITE (*,*)
C          .         '     Range                    (km)  = ', DIST
C
C                 ELSE
C
C                    CALL TOSTDO ( '   Surface representation: '
C          .         //            SRFTYP(K)                     )
C                    WRITE (*,*) '     Intercept not found.'
C                    WRITE (*,*) ' '
C
C                 END IF
C
C              END DO
C
C              WRITE (*,*) ' '
C
C           END DO
C
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was:
C
C
C        Surface Intercept Locations for Camera
C        FOV Boundary and Boresight Vectors
C
C           Instrument:            MGS_MOC_NA
C           Epoch:                 2003 OCT 13 06:00:00 UTC
C           Aberration correction: CN+S
C
C
C        Corner vector 1
C
C          Vector in MGS_MOC_NA frame =
C            0.00000185713838   -0.00380156226592    0.99999277403434
C
C          Intercept:
C
C            Surface representation: Ellipsoid
C
C             Radius                   (km)  =    3384.9411357607282
C             Planetocentric Latitude  (deg) =   -48.477482367206768
C             Planetocentric Longitude (deg) =   -123.47407481971256
C             Range                    (km)  =    388.98308225698992
C
C            Surface representation: MGS/MOLA topography, 4 pixel/deg
C
C             Radius                   (km)  =    3387.6408267726060
C             Planetocentric Latitude  (deg) =   -48.492259559975274
C             Planetocentric Longitude (deg) =   -123.47541193495911
C             Range                    (km)  =    386.14510040407890
C
C
C        Corner vector 2
C
C          Vector in MGS_MOC_NA frame =
C            0.00000185713838    0.00380156226592    0.99999277403434
C
C          Intercept:
C
C            Surface representation: Ellipsoid
C
C             Radius                   (km)  =    3384.9396985743228
C             Planetocentric Latitude  (deg) =   -48.481636778911913
C             Planetocentric Longitude (deg) =   -123.39881874871132
C             Range                    (km)  =    388.97510005267645
C
C            Surface representation: MGS/MOLA topography, 4 pixel/deg
C
C             Radius                   (km)  =    3387.6403704507966
C             Planetocentric Latitude  (deg) =   -48.496386688872484
C             Planetocentric Longitude (deg) =   -123.40074354811055
C             Range                    (km)  =    386.13616443321536
C
C
C        Corner vector 3
C
C          Vector in MGS_MOC_NA frame =
C           -0.00000185713838    0.00380156226592    0.99999277403434
C
C          Intercept:
C
C            Surface representation: Ellipsoid
C
C             Radius                   (km)  =    3384.9396897286833
C             Planetocentric Latitude  (deg) =   -48.481662348858336
C             Planetocentric Longitude (deg) =   -123.39882195503854
C             Range                    (km)  =    388.97464113550637
C
C            Surface representation: MGS/MOLA topography, 4 pixel/deg
C
C             Radius                   (km)  =    3387.6403603146173
C             Planetocentric Latitude  (deg) =   -48.496412042429789
C             Planetocentric Longitude (deg) =   -123.40074672915324
C             Range                    (km)  =    386.13571069851952
C
C
C        Corner vector 4
C
C          Vector in MGS_MOC_NA frame =
C           -0.00000185713838   -0.00380156226592    0.99999277403434
C
C          Intercept:
C
C            Surface representation: Ellipsoid
C
C             Radius                   (km)  =    3384.9411269137699
C             Planetocentric Latitude  (deg) =   -48.477507940479093
C             Planetocentric Longitude (deg) =   -123.47407797517752
C             Range                    (km)  =    388.98262331952731
C
C            Surface representation: MGS/MOLA topography, 4 pixel/deg
C
C             Radius                   (km)  =    3387.6408166344654
C             Planetocentric Latitude  (deg) =   -48.492284916898356
C             Planetocentric Longitude (deg) =   -123.47541506563026
C             Range                    (km)  =    386.14464664863749
C
C
C        Boresight vector
C
C          Vector in MGS_MOC_NA frame =
C            0.00000000000000    0.00000000000000    1.00000000000000
C
C          Intercept:
C
C            Surface representation: Ellipsoid
C
C             Radius                   (km)  =    3384.9404100068609
C             Planetocentric Latitude  (deg) =   -48.479580262226833
C             Planetocentric Longitude (deg) =   -123.43644973546644
C             Range                    (km)  =    388.97571440620783
C
C            Surface representation: MGS/MOLA topography, 4 pixel/deg
C
C             Radius                   (km)  =    3387.6402755067679
C             Planetocentric Latitude  (deg) =   -48.494341863340743
C             Planetocentric Longitude (deg) =   -123.43808042359795
C             Range                    (km)  =    386.13761526562394
C
C
C 
C     2) Use SUBPNT to find the sub-spacecraft point on Mars for the
C        Mars Reconnaissance Orbiter spacecraft (MRO) at a specified
C        time, using the "near point: ellipsoid" computation method.
C        Use both LT+S and CN+S aberration corrections to illustrate
C        the differences.
C
C        Convert the spacecraft to sub-observer point vector obtained
C        from SUBPNT into the MRO_HIRISE_LOOK_DIRECTION reference frame
C        at the observation time. Perform a consistency check with this
C        vector: compare the Mars surface intercept of the ray
C        emanating from the spacecraft and pointed along this vector
C        with the sub-observer point.
C
C        Perform the sub-observer point and surface intercept
C        computations using both triaxial ellipsoid and topographic
C        surface models. 
C
C        For this example, the topographic model is based on the MGS
C        MOLA DEM megr90n000eb, which has a resolution of 16
C        pixels/degree. Eight DSKs, each covering longitude and
C        latitude ranges of 90 degrees, were made from this data set.
C        For the region covered by a given DSK, a grid of approximately
C        1500 x 1500 interpolated heights was produced, and this grid
C        was tessellated using approximately 4.5 million triangular
C        plates, giving a total plate count of about 36 million for the
C        entire DSK set.
C
C        All DSKs in the set use the surface ID code 499001, so there
C        is no need to specify the surface ID in the METHOD strings
C        passed to SINCPT and SUBPNT.
C
C        Use the meta-kernel shown below to load the required SPICE
C        kernels.
C
C
C           KPL/MK
C
C           File: sincpt_ex2.tm
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
C              de430.bsp                        Planetary ephemeris
C              mar097.bsp                       Mars satellite ephemeris
C              pck00010.tpc                     Planet orientation and
C                                               radii
C              naif0011.tls                     Leapseconds
C              mro_psp4_ssd_mro95a.bsp          MRO ephemeris
C              mro_v11.tf                       MRO frame specifications
C              mro_sclkscet_00022_65536.tsc     MRO SCLK coefficients
C                                               parameters
C              mro_sc_psp_070925_071001.bc      MRO attitude
C              megr90n000eb_*_plate.bds         Plate model DSKs based 
C                                               on MEGDR DEM, resolution
C                                               16 pixels/degree.
C
C           \begindata
C
C              KERNELS_TO_LOAD = ( 
C
C                 'de430.bsp',
C                 'mar097.bsp',
C                 'pck00010.tpc',
C                 'naif0011.tls',
C                 'mro_psp4_ssd_mro95a.bsp',
C                 'mro_v11.tf',
C                 'mro_sclkscet_00022_65536.tsc',
C                 'mro_sc_psp_070925_071001.bc',
C                 'megr90n000eb_LL000E00N_UR090E90N_plate.bds'
C                 'megr90n000eb_LL000E90S_UR090E00S_plate.bds'
C                 'megr90n000eb_LL090E00N_UR180E90N_plate.bds'
C                 'megr90n000eb_LL090E90S_UR180E00S_plate.bds'
C                 'megr90n000eb_LL180E00N_UR270E90N_plate.bds'
C                 'megr90n000eb_LL180E90S_UR270E00S_plate.bds'
C                 'megr90n000eb_LL270E00N_UR360E90N_plate.bds'
C                 'megr90n000eb_LL270E90S_UR360E00S_plate.bds'  )
C                                  
C           \begintext
C
C
C
C       Example code begins here.
C          
C
C           PROGRAM EX2
C           IMPLICIT NONE
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      DPR
C           DOUBLE PRECISION      VDIST
C           DOUBLE PRECISION      VNORM
C
C     C
C     C     Local parameters
C     C
C           CHARACTER*(*)         META
C           PARAMETER           ( META   = 'sincpt_ex2.tm' )
C
C           CHARACTER*(*)         F1
C           PARAMETER           ( F1     = '(A,F21.9)' )
C
C           CHARACTER*(*)         F2
C           PARAMETER           ( F2     = '(A)' )
C
C           INTEGER               FRNMLN
C           PARAMETER           ( FRNMLN = 32 )
C
C           INTEGER               MTHLEN
C           PARAMETER           ( MTHLEN = 50 )
C
C           INTEGER               CORLEN
C           PARAMETER           ( CORLEN = 5 )
C
C           INTEGER               NCORR
C           PARAMETER           ( NCORR  = 2 )
C
C           INTEGER               NMETH
C           PARAMETER           ( NMETH  = 2 )
C
C     C
C     C     Local variables
C     C
C           CHARACTER*(CORLEN)    ABCORR ( NCORR )
C           CHARACTER*(FRNMLN)    FIXREF
C           CHARACTER*(FRNMLN)    HIREF
C           CHARACTER*(MTHLEN)    SINMTH ( NMETH )
C           CHARACTER*(MTHLEN)    SUBMTH ( NMETH )
C
C           DOUBLE PRECISION      ALT
C           DOUBLE PRECISION      ET
C           DOUBLE PRECISION      LAT
C           DOUBLE PRECISION      LON
C           DOUBLE PRECISION      MROVEC ( 3 )
C           DOUBLE PRECISION      RADIUS
C           DOUBLE PRECISION      SPOINT ( 3 )
C           DOUBLE PRECISION      SRFVEC ( 3 )
C           DOUBLE PRECISION      TRGEPC
C           DOUBLE PRECISION      XFORM  ( 3, 3 )
C           DOUBLE PRECISION      XEPOCH
C           DOUBLE PRECISION      XPOINT ( 3 )
C           DOUBLE PRECISION      XVEC   ( 3 )
C
C           INTEGER               I
C           INTEGER               J
C
C           LOGICAL               FOUND
C
C     C
C     C     Initial values
C     C
C           DATA                  ABCORR / 'LT+S', 'CN+S'            /
C           DATA                  FIXREF / 'IAU_MARS'                /
C           DATA                  SINMTH / 'Ellipsoid',
C          .                               'DSK/Unprioritized'       /
C           DATA                  SUBMTH / 'Ellipsoid/Near point',
C          .                               'DSK/Unprioritized/Nadir' /
C
C     C
C     C     Load kernel files via the meta-kernel.
C     C
C           CALL FURNSH ( META )
C
C     C
C     C     Convert the TDB request time string to seconds past
C     C     J2000, TDB.
C     C
C           CALL STR2ET ( '2007 SEP 30 00:00:00 TDB', ET )
C
C     C
C     C     Compute the sub-spacecraft point using the
C     C     "NEAR POINT: ELLIPSOID" definition.
C     C     Compute the results using both LT+S and CN+S
C     C     aberration corrections.
C     C
C     C     Repeat the computation for each method.
C     C
C     C
C           DO I = 1, NMETH
C
C              WRITE(*,F2) ' '
C              WRITE(*,F2) 'Sub-observer point computation method = '
C          .               // SUBMTH(I)
C
C              DO J = 1, NCORR
C
C                 CALL SUBPNT ( SUBMTH(I),
C          .                    'Mars', ET,     FIXREF, ABCORR(J),
C          .                    'MRO',  SPOINT, TRGEPC, SRFVEC    )
C     C
C     C           Compute the observer's altitude above SPOINT.
C     C
C                 ALT = VNORM ( SRFVEC )
C     C
C     C           Express SRFVEC in the MRO_HIRISE_LOOK_DIRECTION
C     C           reference frame at epoch ET. Since SRFVEC is
C     c           expressed relative to the IAU_MARS frame at
C     C           TRGEPC, we must call PXFRM2 to compute the position
C     C           transformation matrix from IAU_MARS at TRGEPC to the
C     C           MRO_HIRISE_LOOK_DIRECTION frame at time ET.
C     C
C     C           To make code formatting a little easier, we'll
C     C           store the long MRO reference frame name in a
C     C           variable:
C     C
C                 HIREF = 'MRO_HIRISE_LOOK_DIRECTION'
C
C                 CALL PXFRM2 ( FIXREF, HIREF,  TRGEPC, ET, XFORM )
C                 CALL MXV    ( XFORM,  SRFVEC, MROVEC )
C
C     C
C     C           Convert rectangular coordinates to planetocentric
C     C           latitude and longitude. Convert radians to degrees.
C     C
C                 CALL RECLAT ( SPOINT, RADIUS, LON, LAT  )
C
C                 LON = LON * DPR ()
C                 LAT = LAT * DPR ()
C     C
C     C           Write the results.
C     C
C                 WRITE(*,F2) ' '
C                 WRITE(*,F2) '   Aberration correction = '//ABCORR(J)
C                 WRITE(*,F1) ' '
C                 WRITE(*,F2) '      MRO-to-sub-observer vector in'
C                 WRITE(*,F2) '      MRO HIRISE look direction frame'
C                 WRITE(*,F1) '        X-component             (km) = ',
C          .                  MROVEC(1)
C                 WRITE(*,F1) '        Y-component             (km) = ',
C          .                  MROVEC(2)
C                 WRITE(*,F1) '        Z-component             (km) = ',
C          .               MROVEC(3)
C                 WRITE(*,F1) '      Sub-observer point radius (km) = ',
C          .                  RADIUS
C                 WRITE(*,F1) '      Planetocentric latitude  (deg) = ',
C          .                  LAT
C                 WRITE(*,F1) '      Planetocentric longitude (deg) = ',
C          .                  LON
C                 WRITE(*,F1) '      Observer altitude         (km) = ',
C          .                  ALT
C
C     C
C     C           Consistency check: find the surface intercept on
C     C           Mars of the ray emanating from the spacecraft and
C     C           having direction vector MROVEC in the MRO HIRISE
C     C           reference frame at ET. Call the intercept point
C     C           XPOINT. XPOINT should coincide with SPOINT, up to
C     C           a small round-off error.
C     C
C                 CALL SINCPT ( SINMTH(I), 'Mars', ET,    FIXREF,
C          .                    ABCORR(J), 'MRO',  HIREF, MROVEC,
C          .                    XPOINT,    XEPOCH, XVEC,  FOUND  )
C
C                 IF ( .NOT. FOUND ) THEN
C                    WRITE (*,F1) 'Bug: no intercept'
C                 ELSE
C     C
C     C              Report the distance between XPOINT and SPOINT.
C     C
C                    WRITE (*,* ) ' '
C                    WRITE (*,F1) '   Intercept comparison '
C          .         //           'error (km) = ',
C          .                      VDIST( XPOINT, SPOINT )
C                 END IF
C
C              END DO
C
C           END DO
C
C           END
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit 
C     platform, the output was:
C           
C
C        Sub-observer point computation method = Ellipsoid/Near point
C
C           Aberration correction = LT+S
C
C              MRO-to-sub-observer vector in
C              MRO HIRISE look direction frame
C                X-component             (km) =           0.286933229
C                Y-component             (km) =          -0.260425939
C                Z-component             (km) =         253.816326385
C              Sub-observer point radius (km) =        3388.299078378
C              Planetocentric latitude  (deg) =         -38.799836378
C              Planetocentric longitude (deg) =        -114.995297227
C              Observer altitude         (km) =         253.816622175
C
C           Intercept comparison error (km) =           0.000002144
C
C           Aberration correction = CN+S
C
C              MRO-to-sub-observer vector in
C              MRO HIRISE look direction frame
C                X-component             (km) =           0.286933107
C                Y-component             (km) =          -0.260426683
C                Z-component             (km) =         253.816315915
C              Sub-observer point radius (km) =        3388.299078376
C              Planetocentric latitude  (deg) =         -38.799836382
C              Planetocentric longitude (deg) =        -114.995297449
C              Observer altitude         (km) =         253.816611705
C
C           Intercept comparison error (km) =           0.000000001
C
C        Sub-observer point computation method = DSK/Unprioritized/Nadir
C
C           Aberration correction = LT+S
C
C              MRO-to-sub-observer vector in
C              MRO HIRISE look direction frame
C                X-component             (km) =           0.282372596
C                Y-component             (km) =          -0.256289313
C                Z-component             (km) =         249.784871247
C              Sub-observer point radius (km) =        3392.330239436
C              Planetocentric latitude  (deg) =         -38.800230156
C              Planetocentric longitude (deg) =        -114.995297338
C              Observer altitude         (km) =         249.785162334
C
C           Intercept comparison error (km) =           0.000002412
C
C           Aberration correction = CN+S
C
C              MRO-to-sub-observer vector in
C              MRO HIRISE look direction frame
C                X-component             (km) =           0.282372464
C                Y-component             (km) =          -0.256290075
C                Z-component             (km) =         249.784860121
C              Sub-observer point radius (km) =        3392.330239564
C              Planetocentric latitude  (deg) =         -38.800230162
C              Planetocentric longitude (deg) =        -114.995297569
C              Observer altitude         (km) =         249.785151209
C
C           Intercept comparison error (km) =           0.000000001
C
C
C$ Restrictions
C
C     A cautionary note: if aberration corrections are used, and 
C     if DREF is the target body-fixed frame, the epoch at which that
C     frame is evaluated is offset from ET by the light time between
C     the observer and the *center* of the target body. This light time
C     normally will differ from the light time between the observer and
C     intercept point. Consequently the orientation of the target
C     body-fixed frame at TRGEPC will not match that of the target
C     body-fixed frame at the epoch associated with DREF. As a result,
C     various derived quantities may not be as expected: for example,
C     SRFVEC would not be parallel to DVEC.
C
C     In many applications the errors arising from this frame
C     discrepancy may be insignificant; however a safe approach is to
C     always use as DREF a frame other than the target body-fixed
C     frame.
C     
C$ Literature_References
C
C     None.
C
C$ Author_and_Institution
C
C     N.J. Bachman   (JPL)
C     S.C. Krening   (JPL)
C     B.V. Semenov   (JPL)
C
C$ Version
C
C-    SPICELIB Version 3.0.0, 04-APR-2017 (NJB)
C
C       01-FEB-2016 (NJB)
C
C        Upgraded to support surfaces represented by DSKs. 
C
C        Updated kernels are used in header example programs.
C     
C-    SPICELIB Version 2.0.0, 31-MAR-2014 (NJB) (SCK) (BVS)
C
C        Bug fix: FIRST is now set to .FALSE. at the completion
C        of a successful initialization pass. This does not affect
C        the routine's outputs but improves efficiency.
C
C        Bug fix: redundant call to SPKSSB was removed. This does not
C        affect the routine's outputs but improves efficiency.
C
C        References to the new PXFRM2 routine were added, which changed
C        the Detailed Output section and the second example. Some header
C        comment corrections were made.
C
C        Upgrade: this routine now uses ZZVALCOR rather than
C        ZZPRSCOR, simplifying the implementation.
C
C        Upgrade: this routine now saves the input body names and
C        ZZBODTRN state counters and does name-ID conversions only if
C        the counters have changed.
C
C        Upgrade: this routine now saves the input frame names and POOL
C        state counters and does frame name-ID conversions only if the
C        counters have changed.
C
C-    SPICELIB Version 1.2.0, 07-APR-2010 (NJB)
C
C        Code style improvement: re-use of variables in 
C        FRINFO calls has been eliminated. There is no impact 
C        of the behavior of the routine.
C
C-    SPICELIB Version 1.1.0, 17-MAR-2009 (NJB)(EDW) 
C
C        Bug fix: quick test for non-intersection is
C        no longer performed when observer-target distance
C        is less than target's maximum radius.
C
C        Typos in the Detailed Input section's description of DREF
C        were corrected.
C
C        In the header examples, meta-kernel names were updated to use
C        the suffix 
C
C           ".tm"
C
C        Incorrect frame name FIXFRM was changed to FIXREF in
C        documentation.
C
C        Typo correction in Required_Reading, changed FRAME 
C        to FRAMES.
C
C-    SPICELIB Version 1.0.0, 02-MAR-2008 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     find surface intercept point
C     find intersection of ray and target body surface
C     find intercept of ray on target body surface
C
C-&
 

C$ Revisions
C
C-    SPICELIB Version 3.0.0, 01-FEB-2016 (NJB)
C
C        Upgraded to support surfaces represented by DSKs.
C     
C        The routine was re-written so as to use a private
C        routine to implement the intersection algorithm.
C        That routine has been generalized so that it does
C        not depend on the target surface representation: it
C        uses callback routines to compute ray-surface intercepts
C        for a specified ray and time, the surface tangency point
C        for a given ray, and the radius of an outer bounding
C        sphere for the target.
C
C-&

 
C
C     SPICELIB functions
C
      LOGICAL               EQSTR
      LOGICAL               FAILED
      LOGICAL               RETURN
      LOGICAL               VZERO

C
C     EXTERNAL routines
C 
      EXTERNAL              ZZMAXRAD
      EXTERNAL              ZZRAYSFX
      EXTERNAL              ZZRAYNP

C
C     Local parameters
C     
      CHARACTER*(*)         RNAME
      PARAMETER           ( RNAME  = 'SINCPT' )

C
C     Saved body name length.
C
      INTEGER               MAXL
      PARAMETER           ( MAXL   = 36 )

C
C     Saved frame name length.
C
      INTEGER               FRNMLN
      PARAMETER           ( FRNMLN = 32 )

C
C     Local variables
C
      CHARACTER*(CVTLEN)    PNTDEF
      CHARACTER*(CORLEN)    PRVCOR
      CHARACTER*(MTHLEN)    PRVMTH
      CHARACTER*(SHPLEN)    SHPSTR
      CHARACTER*(SUBLEN)    SUBTYP
      CHARACTER*(TMTLEN)    TRMSTR

      INTEGER               DCENTR
      INTEGER               DCLASS
      INTEGER               DFRCDE
      INTEGER               DTYPID
      INTEGER               FXCENT
      INTEGER               FXCLSS
      INTEGER               FXFCDE
      INTEGER               FXTYID
      INTEGER               NSURF
      INTEGER               OBSCDE
      INTEGER               SHAPE
      INTEGER               SRFLST ( MAXSRF )
      INTEGER               TRGCDE
      
      LOGICAL               ATTBLK ( NABCOR )
      LOGICAL               FIRST
      LOGICAL               FND
      LOGICAL               PRI
      LOGICAL               SURFUP
      LOGICAL               USECN
      LOGICAL               USELT
      LOGICAL               USESTL
      LOGICAL               XMIT

C
C     Saved name/ID item declarations.
C
      INTEGER               SVCTR1 ( CTRSIZ )
      CHARACTER*(MAXL)      SVTARG
      INTEGER               SVTCDE
      LOGICAL               SVFND1

      INTEGER               SVCTR2 ( CTRSIZ )
      CHARACTER*(MAXL)      SVOBSR
      INTEGER               SVOBSC
      LOGICAL               SVFND2

C
C     Saved frame name/ID item declarations.
C
      INTEGER               SVCTR3 ( CTRSIZ )
      CHARACTER*(FRNMLN)    SVFREF
      INTEGER               SVFXFC

      INTEGER               SVCTR4 ( CTRSIZ )
      CHARACTER*(FRNMLN)    SVDREF
      INTEGER               SVDFRC

C
C     Saved surface name/ID item declarations.
C
      INTEGER               SVCTR5 ( CTRSIZ )

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  NSURF
      SAVE                  PRI
      SAVE                  PRVCOR
      SAVE                  PRVMTH
      SAVE                  SHAPE
      SAVE                  SRFLST
      SAVE                  USECN
      SAVE                  USELT
      SAVE                  USESTL
      SAVE                  XMIT

C
C     Saved name/ID items.
C
      SAVE                  SVCTR1
      SAVE                  SVTARG
      SAVE                  SVTCDE
      SAVE                  SVFND1

      SAVE                  SVCTR2
      SAVE                  SVOBSR
      SAVE                  SVOBSC
      SAVE                  SVFND2

C
C     Saved frame name/ID items.
C
      SAVE                  SVCTR3
      SAVE                  SVFREF
      SAVE                  SVFXFC

      SAVE                  SVCTR4
      SAVE                  SVDREF
      SAVE                  SVDFRC

C
C     Saved surface name/ID items.
C
      SAVE                  SVCTR5
 
C
C     Initial values
C
      DATA                  FIRST  / .TRUE.  /
      DATA                  PRVCOR / ' '     /
      DATA                  PRVMTH / ' '     /
      DATA                  USECN  / .FALSE. / 
      DATA                  USELT  / .FALSE. / 
      DATA                  USESTL / .FALSE. / 
      DATA                  XMIT   / .FALSE. / 

C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( RNAME )
  
C
C     Nothing has been found yet.
C
      FOUND = .FALSE.

C
C     Counter initialization is done separately.
C
      IF ( FIRST ) THEN
C
C        Initialize counters.
C
         CALL ZZCTRUIN( SVCTR1 )
         CALL ZZCTRUIN( SVCTR2 )
         CALL ZZCTRUIN( SVCTR3 )
         CALL ZZCTRUIN( SVCTR4 )
         CALL ZZCTRUIN( SVCTR5 )

      END IF

      IF (  FIRST  .OR.  ( ABCORR .NE. PRVCOR )  ) THEN
C
C        Make sure the results of this block won't be reused
C        if we bail out due to an error.
C
         PRVCOR = ' '

C
C        The aberration correction flag differs from the value it
C        had on the previous call, if any. Analyze the new flag.
C
         CALL ZZVALCOR ( ABCORR, ATTBLK )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF
C
C        Set logical flags indicating the attributes of the requested
C        correction:
C
C           XMIT is .TRUE. when the correction is for transmitted
C           radiation.
C
C           USELT is .TRUE. when any type of light time correction
C           (normal or converged Newtonian) is specified.
C
C           USECN indicates converged Newtonian light time correction.
C
C           USESTL indicates stellar aberration corrections.
C
C
C        The above definitions are consistent with those used by
C        ZZPRSCOR.
C 
         XMIT    =  ATTBLK ( XMTIDX )
         USELT   =  ATTBLK ( LTIDX  )
         USECN   =  ATTBLK ( CNVIDX )
         USESTL  =  ATTBLK ( STLIDX )
C
C        The aberration correction flag is valid; save it.
C
         PRVCOR = ABCORR

      END IF

C
C     Obtain integer codes for the target and observer.
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
      
      
      CALL ZZBODS2C ( SVCTR2, SVOBSR, SVOBSC, SVFND2,
     .                OBSRVR, OBSCDE, FND    )
      
      IF ( .NOT. FND ) THEN

         CALL SETMSG ( 'The observer, '
     .   //            '''#'', is not a recognized name for an '
     .   //            'ephemeris object. The cause of this '
     .   //            'problem may be that you need an updated '
     .   //            'version of the SPICE Toolkit, or that you '
     .   //            'failed to load a kernel containing a '
     .   //            'name-ID mapping for this body.'           )
         CALL ERRCH  ( '#', OBSRVR                                )
         CALL SIGERR ( 'SPICE(IDCODENOTFOUND)'                    )
         CALL CHKOUT ( RNAME                                      )
         RETURN
      
      END IF
      
 
C
C     Check the input body codes. If they are equal, signal
C     an error.
C
      IF ( OBSCDE .EQ. TRGCDE ) THEN
 
         CALL SETMSG ( 'In computing the surface intercept point, ' //
     .                 'the observing body and target body are the '//
     .                 'same. Both are #.'                          )
         CALL ERRCH  ( '#',  OBSRVR                                 )
         CALL SIGERR ( 'SPICE(BODIESNOTDISTINCT)'                   )
         CALL CHKOUT ( RNAME                                        )
         RETURN
 
      END IF

C
C     Determine the attributes of the frame designated by FIXREF.
C
      CALL ZZNAMFRM ( SVCTR3, SVFREF, SVFXFC, FIXREF, FXFCDE )

      CALL FRINFO ( FXFCDE, FXCENT, FXCLSS, FXTYID, FND )

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
C     Check for a zero ray direction vector.
C     
      IF ( VZERO(DVEC) ) THEN

         CALL SETMSG ( 'Input ray direction was the zero vector; this '
     .   //            'vector must be non-zero.'                      )
         CALL SIGERR ( 'SPICE(ZEROVECTOR)'                             )
         CALL CHKOUT ( RNAME                                           )
         RETURN

      END IF
   
C
C     Determine the attributes of the frame designated by DREF.
C
      CALL ZZNAMFRM ( SVCTR4, SVDREF, SVDFRC, DREF, DFRCDE )

      CALL FRINFO ( DFRCDE, DCENTR, DCLASS, DTYPID, FND )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

      IF ( .NOT. FND ) THEN

         CALL SETMSG ( 'Reference frame # is not recognized by ' //
     .                 'the SPICE frame subsystem. Possibly '    //
     .                 'a required frame definition kernel has ' //
     .                 'not been loaded.'                        )
         CALL ERRCH  ( '#',  DREF                                )
         CALL SIGERR ( 'SPICE(NOFRAME)'                          )
         CALL CHKOUT ( RNAME                                     )
         RETURN

      END IF

C
C     Check whether the surface name/ID mapping has been updated.
C
      CALL ZZSRFTRK ( SVCTR5, SURFUP )

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
     .                   PRI,    NSURF,  SRFLST, PNTDEF, TRMSTR )
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
     .      //            'SUBSLR, but is not applicable for SINCPT.' )
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


      IF ( SHAPE .EQ. ELLSHP ) THEN
C
C        Initialize the intercept algorithm to use the reference
C        ellipsoid of the target body. 
C        
         CALL ZZSUELIN ( TRGCDE )


      ELSE IF ( SHAPE .EQ. DSKSHP ) THEN
C
C        This is the DSK case.
C
C        If the method string listed a set of surface IDs, NSURF is
C        positive and SRFLST contains those IDs.
C
C        Initialize the intercept algorithm to use a DSK
C        model for the surface of the target body. 
C        
         CALL ZZSUDSKI ( TRGCDE, NSURF, SRFLST, FXFCDE )

      ELSE
C
C        This is a backstop check.
C
         CALL SETMSG ( '[2] Returned shape value from method '
     .   //            'string was <#>.'                      )
         CALL ERRCH  ( '#', SHPSTR                            )
         CALL SIGERR ( 'SPICE(BUG)'                           )
         CALL CHKOUT ( RNAME                                  )
         RETURN

      END IF

C
C     Perform the intercept computation.
C
      CALL ZZSFXCOR ( ZZRAYNP, ZZMAXRAD, ZZRAYSFX, TRGCDE,  
     .                ET,      ABCORR,   USELT,    USECN,        
     .                USESTL,  XMIT,     FIXREF,   OBSCDE,   
     .                DFRCDE,  DCLASS,   DCENTR,   DVEC,    
     .                SPOINT,  TRGEPC,   SRFVEC,   FOUND  )
     
      CALL CHKOUT ( RNAME )
      RETURN
      END
