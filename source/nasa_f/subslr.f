C$Procedure SUBSLR ( Sub-solar point )
 
      SUBROUTINE SUBSLR ( METHOD, TARGET, ET,     FIXREF,  
     .                    ABCORR, OBSRVR, SPOINT, TRGEPC, SRFVEC )
      
C$ Abstract
C
C     Compute the rectangular coordinates of the sub-solar point on
C     a target body at a specified epoch, optionally corrected for
C     light time and stellar aberration.
C
C     The surface of the target body may be represented by a triaxial
C     ellipsoid or by topographic data provided by DSK files.
C
C     This routine supersedes SUBSOL.
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
C     NAIF_IDS
C     PCK
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
      DOUBLE PRECISION      SPOINT ( 3 )
      DOUBLE PRECISION      TRGEPC
      DOUBLE PRECISION      SRFVEC ( 3 )
 
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
C     SPOINT     O   Sub-solar point on the target body.
C     TRGEPC     O   Sub-solar point epoch.
C     SRFVEC     O   Vector from observer to sub-solar point.
C
C$ Detailed_Input
C
C     METHOD   is a short string providing parameters defining
C              the computation method to be used. In the syntax
C              descriptions below, items delimited by brackets
C              are optional.
C
C              METHOD may be assigned the following values:
C
C                 'NEAR POINT/ELLIPSOID'
C
C                    The sub-solar point computation uses a triaxial
C                    ellipsoid to model the surface of the target body.
C                    The sub-solar point is defined as the nearest
C                    point on the target relative to the sun.
C
C                    The word "NADIR" may be substituted for the phrase
C                    "NEAR POINT" in the string above. 
C
C                    For backwards compatibility, the older syntax
C
C                       'Near point: ellipsoid'
C
C                    is accepted as well.
C
C
C                 'INTERCEPT/ELLIPSOID'
C
C                    The sub-solar point computation uses a triaxial
C                    ellipsoid to model the surface of the target body.
C                    The sub-solar point is defined as the target
C                    surface intercept of the line containing the sun
C                    and the target's center.
C
C                    For backwards compatibility, the older syntax
C
C                       'Intercept: ellipsoid'
C
C                    is accepted as well.
C
C
C                 'NADIR/DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                    The sub-solar point computation uses DSK data to
C                    model the surface of the target body. The
C                    sub-solar point is defined as the intercept, on
C                    the surface represented by the DSK data, of the
C                    line containing the sun and the nearest point on
C                    the target's reference ellipsoid. If multiple such
C                    intercepts exist, the one closest to the sun is
C                    selected.
C
C                    Note that this definition of the sub-solar point
C                    is not equivalent to the "nearest point on the
C                    surface to the sun."  The phrase "NEAR POINT" may
C                    NOT be substituted for "NADIR" in the string
C                    above.
C
C                    The surface list specification is optional. The
C                    syntax of the list is
C
C                       <surface 1> [, <surface 2>...]
C
C                    If present, it indicates that data only for the
C                    listed surfaces are to be used; however, data
C                    need not be available for all surfaces in the
C                    list. If absent, loaded DSK data for any surface
C                    associated with the target body are used.
C
C                    The surface list may contain surface names or
C                    surface ID codes. Names containing blanks must
C                    be delimited by double quotes, for example
C
C                       SURFACES = "Mars MEGDR 128 PIXEL/DEG"
C
C                    If multiple surfaces are specified, their names
C                    or IDs must be separated by commas.
C
C                    See the Particulars section below for details
C                    concerning use of DSK data.
C
C
C                 'INTERCEPT/DSK/UNPRIORITIZED[/SURFACES = 
C                                              <surface list>]'
C
C                    The sub-solar point computation uses DSK data to
C                    model the surface of the target body. The
C                    sub-solar point is defined as the target surface
C                    intercept of the line containing the sun and the
C                    target's center.
C
C                    If multiple such intercepts exist, the one closest
C                    to the sun is selected.
C
C                    The surface list specification is optional. The
C                    syntax of the list is identical to that for the
C                    NADIR option described above.
C
C
C                 Neither case nor white space are significant in
C                 METHOD, except within double-quoted strings. For
C                 example, the string ' eLLipsoid/nearpoint ' is valid.
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
C     TARGET      is the name of the target body. The target body is 
C                 an ephemeris object (its trajectory is given by
C                 SPK data), and is an extended object.
C
C                 The string TARGET is case-insensitive, and leading
C                 and trailing blanks in TARGET are not significant.
C                 Optionally, you may supply a string containing the
C                 integer ID code for the object. For example both
C                 'MOON' and '301' are legitimate strings that indicate
C                 the Moon is the target body.
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
C                 the epoch at which the position and orientation of
C                 the target body and the position of the Sun are
C                 computed.
C
C                 When aberration corrections are used, ET is the epoch
C                 at which the observer's state relative to the solar
C                 system barycenter is computed; in this case the
C                 position and orientation of the target body are
C                 computed at ET-LT, where LT is the one-way light time
C                 between the sub-solar point and the observer. See the
C                 description of ABCORR below for details.
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
C                 The output sub-solar point SPOINT and the
C                 observer-to-sub-solar point vector SRFVEC will be
C                 expressed relative to this reference frame.
C
C                 
C     ABCORR      indicates the aberration correction to be applied
C                 when computing the target position and orientation
C                 and the position of the Sun.
C
C                 For remote sensing applications, where the apparent
C                 sub-solar point seen by the observer is desired,
C                 normally either of the corrections
C              
C                    'LT+S' 
C                    'CN+S'
C     
C                 should be used. These and the other supported options
C                 are described below. ABCORR may be any of the 
C                 following:
C
C                    'NONE'     Apply no correction. Return the
C                               geometric sub-solar point on the target
C                               body.
C
C                 Let LT represent the one-way light time between the
C                 observer and the sub-solar point (note: NOT between
C                 the observer and the target body's center). The
C                 following values of ABCORR apply to the "reception"
C                 case in which photons depart from the sub-solar
C                 point's location at the light-time corrected epoch
C                 ET-LT and *arrive* at the observer's location at ET:
C
C                    'LT'       Correct for one-way light time (also
C                               called "planetary aberration") using a
C                               Newtonian formulation. This correction
C                               yields the location of sub-solar
C                               point at the moment it emitted photons
C                               arriving at the observer at ET.
C 
C                               The light time correction uses an
C                               iterative solution of the light time
C                               equation. The solution invoked by the
C                               'LT' option uses one iteration.
C
C                               The target position and orientation as
C                               seen by the observer are corrected for
C                               light time. The position of the Sun
C                               relative to the target is corrected for
C                               one-way light time between the Sun and
C                               target.
C
C                    'LT+S'     Correct for one-way light time and
C                               stellar aberration using a Newtonian
C                               formulation. This option modifies the
C                               sub-solar point obtained with the 'LT'
C                               option to account for the observer's
C                               velocity relative to the solar system
C                               barycenter. These corrections yield
C                               the apparent sub-solar point.
C
C                    'CN'       Converged Newtonian light time
C                               correction. In solving the light time
C                               equation, the 'CN' correction iterates
C                               until the solution converges. Both the
C                               position and rotation of the target
C                               body, and the position of the Sun, are
C                               corrected for light time.
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
C                 Neither case nor white space are significant in
C                 ABCORR. For example, the string 
C
C                   'Lt + s'
C
C                 is valid.
C
C
C     OBSRVR      is the name of the observing body. The observing body
C                 is an ephemeris object: it typically is a spacecraft,
C                 the earth, or a surface point on the earth. OBSRVR is
C                 case-insensitive, and leading and trailing blanks in
C                 OBSRVR are not significant. Optionally, you may
C                 supply a string containing the integer ID code for
C                 the object. For example both 'MOON' and '301' are
C                 legitimate strings that indicate the Moon is the
C                 observer.
C
C                 The observer may coincide with the target.
C
C$ Detailed_Output
C
C
C     SPOINT      is the sub-solar point on the target body. 
C
C                 For target shapes modeled by ellipsoids, the
C                 sub-solar point is defined either as the point on the
C                 target body that is closest to the sun, or the target
C                 surface intercept of the line from the sun to the
C                 target's center.
C
C                 For target shapes modeled by topographic data
C                 provided by DSK files, the sub-solar point is defined
C                 as the target surface intercept of the line from the
C                 sun to either the nearest point on the reference
C                 ellipsoid, or to the target's center. If multiple
C                 such intercepts exist, the one closest to the sun is
C                 selected.
C
C                 The input argument METHOD selects the target shape
C                 model and sub-solar point definition to be used.
C 
C                 SPOINT is expressed in Cartesian coordinates,
C                 relative to the body-fixed target frame designated by
C                 FIXREF. The body-fixed target frame is evaluated at
C                 the sub-solar point epoch TRGEPC (see description
C                 below).
C
C                 When aberration corrections are used, SPOINT is
C                 computed using target body position and orientation
C                 that have been adjusted for the corrections
C                 applicable to SPOINT itself rather than to the target
C                 body's center. In particular, if the stellar
C                 aberration correction applicable to SPOINT is
C                 represented by a shift vector S, then the light-time
C                 corrected position of the target is shifted by S
C                 before the sub-solar point is computed.
C                 
C                 The components of SPOINT have units of km.
C
C
C     TRGEPC      is the "sub-solar point epoch." TRGEPC is defined as
C                 follows: letting LT be the one-way light time between
C                 the observer and the sub-solar point, TRGEPC is
C                 either the epoch ET-LT or ET depending on whether the
C                 requested aberration correction is, respectively, for
C                 received radiation or omitted. LT is computed using
C                 the method indicated by ABCORR.
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
C         If transmission aberration corrections are specified, the
C         error SPICE(NOTSUPPORTED) is signaled.
C
C     2)  If either the target or observer input strings cannot be
C         converted to an integer ID code, the error
C         SPICE(IDCODENOTFOUND) is signaled.
C
C     3)  If the input target body-fixed frame FIXREF is not
C         recognized, the error SPICE(NOFRAME) is signaled. A frame
C         name may fail to be recognized because a required frame
C         specification kernel has not been loaded; another cause is a
C         misspelling of the frame name.
C
C     4)  If the input frame FIXREF is not centered at the target body,
C         the error SPICE(INVALIDFRAME) is signaled.
C
C     5)  If the input argument METHOD is not recognized, the error
C         SPICE(INVALIDMETHOD) is signaled by this routine, or the
C         error is diagnosed by a routine in the call tree of this
C         routine.
C
C         If the sub-solar point type is not specified or is not
C         recognized, the error SPICE(INVALIDSUBTYPE) is signaled.
C
C     6)  If insufficient ephemeris data have been loaded prior to
C         calling SUBSLR, the error will be diagnosed and signaled by a
C         routine in the call tree of this routine. Note that when
C         light time correction is used, sufficient ephemeris data must
C         be available to propagate the states of observer, target, and
C         the Sun to the solar system barycenter.
C
C     7)  If the computation method specifies an ellipsoidal target
C         shape and triaxial radii of the target body have not been
C         loaded into the kernel pool prior to calling SUBSLR, the
C         error will be diagnosed and signaled by a routine in the call
C         tree of this routine.
C
C     8)  The target must be an extended body, and must have a shape
C         for which a sub-solar point can be defined.
C
C         If the target body's shape is modeled as an ellipsoid, and if
C         any of the radii of the target body are non-positive, the
C         error will be diagnosed and signaled by routines in the call
C         tree of this routine.
C
C         If the target body's shape is modeled by DSK data, the shape
C         must be such that the specified sub-solar point definition is
C         applicable. For example, if the target shape is a torus, both
C         the NADIR and INTERCEPT definitions might be inapplicable,
C         depending on the relative locations of the sun and target.
C
C     9)  If PCK data specifying the target body-fixed frame
C         orientation have not been loaded prior to calling SUBSLR,
C         the error will be diagnosed and signaled by a routine in the
C         call tree of this routine.
C
C     10) If METHOD specifies that the target surface is represented by
C         DSK data, and no DSK files are loaded for the specified
C         target, the error is signaled by a routine in the call tree
C         of this routine.
C
C     12) If METHOD specifies that the target surface is represented
C         by DSK data, and the ray from the observer to the
C         sub-observer point doesn't intersect the target body's
C         surface, the error SPICE(SUBPOINTNOTFOUND) will be signaled.
C
C     13) In some very rare cases, the surface intercept on the 
C         target body's reference ellipsoid of the observer to target
C         center vector may not be computable. In these cases the
C         error SPICE(DEGENERATECASE) is signaled.
C
C     14) If the target body is the sun, the error SPICE(INVALIDTARGET)
C         is signaled.
C
C$ Files
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - SPK data: ephemeris data for target, observer, and Sun must
C          be loaded. If aberration corrections are used, the states of
C          target, observer, and the Sun relative to the solar system
C          barycenter must be calculable from the available ephemeris
C          data. Typically ephemeris data are made available by loading
C          one or more SPK files via FURNSH.
C
C        - PCK data: rotation data for the target body must be
C          loaded. These may be provided in a text or binary PCK file.
C
C        - Shape data for the target body:
C                
C            PCK data: 
C
C               If the target body shape is modeled as an ellipsoid,
C               triaxial radii for the target body must be loaded into
C               the kernel pool. Typically this is done by loading a
C               text PCK file via FURNSH.
C
C               Triaxial radii are also needed if the target shape is
C               modeled by DSK data, but the DSK NADIR method is
C               selected.
C
C            DSK data: 
C
C               If the target shape is modeled by DSK data, DSK files
C               containing topographic data for the target body must be
C               loaded. If a surface list is specified, data for at
C               least one of the listed surfaces must be loaded.
C
C     The following data may be required:
C
C        - Frame data: if a frame definition is required to convert the
C          observer and target states to the body-fixed frame of the
C          target, that definition must be available in the kernel
C          pool. Typically the definition is supplied by loading a
C          frame kernel via FURNSH.
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
C     In all cases, kernel data are normally loaded once per program
C     run, NOT every time this routine is called.
C
C$ Particulars
C
C     There are two different popular ways to define the sub-solar
C     point: "nearest point on target to the Sun" or "target surface
C     intercept of the line containing the Sun and target." These
C     coincide when the target is spherical and generally are distinct
C     otherwise.
C
C     This routine computes light time corrections using light time
C     between the observer and the sub-solar point, as opposed to the
C     center of the target. Similarly, stellar aberration corrections
C     done by this routine are based on the direction of the vector
C     from the observer to the light-time corrected sub-solar point,
C     not to the target center. This technique avoids errors due to the
C     differential between aberration corrections across the target
C     body. Therefore it's valid to use aberration corrections with
C     this routine even when the observer is very close to the
C     sub-solar point, in particular when the observer to sub-solar
C     point distance is much less than the observer to target center
C     distance.
C     
C     When comparing sub-solar point computations with results from
C     sources other than SPICE, it's essential to make sure the same
C     geometric definitions are used.
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
C           NADIR/DSK/UNPRIORITIZED/<surface list>
C           DSK/NADIR/<surface list>/UNPRIORITIZED
C           UNPRIORITIZED/<surface list>/DSK/NADIR
C
C        The simplest form of the METHOD argument specifying use of
C        DSK data is one that lacks a surface list, for example:
C
C           'NADIR/DSK/UNPRIORITIZED'
C           'INTERCEPT/DSK/UNPRIORITIZED'
C
C        For applications in which all loaded DSK data for the target
C        body are for a single surface, and there are no competing
C        segments, the above strings suffice. This is expected to be
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
C        'NADIR/DSK/UNPRIORITIZED/SURFACES= "Mars MEGDR 64 PIXEL/DEG",3'
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
C
C     The numerical results shown for this example may differ across
C     platforms. The results depend on the SPICE kernels used as input,
C     the compiler and supporting libraries, and the machine specific
C     arithmetic implementation.
C
C
C     1) Find the sub-solar point on Mars as seen from the Earth for a
C        specified time. 
C
C        Compute the sub-solar point using both triaxial ellipsoid
C        and topographic surface models. Topography data are provided by
C        a DSK file. For the ellipsoid model, use both the "intercept"
C        and "near point" sub-observer point definitions; for the DSK
C        case, use both the "intercept" and "nadir" definitions.
C
C        Display the locations of both the sun and the sub-solar
C        point relative to the center of Mars, in the IAU_MARS
C        body-fixed reference frame, using both planetocentric and
C        planetographic coordinates.
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
C           File: subslr_ex1.tm
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
C                                  'megr90n000cb_plate.bds' )
C           \begintext
C
C
C
C       Example code begins here.
C
C
C           PROGRAM EX1
C           IMPLICIT NONE
C     C
C     C     SPICELIB functions
C     C
C           DOUBLE PRECISION      DPR
C     C
C     C     Local parameters
C     C
C           CHARACTER*(*)         META
C           PARAMETER           ( META   = 'subslr_ex1.tm' )
C
C           CHARACTER*(*)         FM
C           PARAMETER           ( FM     =  '(A,F21.9)' )
C
C           INTEGER               MTHLEN
C           PARAMETER           ( MTHLEN = 50 )
C
C           INTEGER               NMETH
C           PARAMETER           ( NMETH  = 4 )
C
C     C
C     C     Local variables
C     C
C           CHARACTER*(MTHLEN)    METHOD ( NMETH )
C
C           DOUBLE PRECISION      ET
C           DOUBLE PRECISION      F
C           DOUBLE PRECISION      RADII  ( 3 )
C           DOUBLE PRECISION      RE
C           DOUBLE PRECISION      RP
C           DOUBLE PRECISION      SPCLAT
C           DOUBLE PRECISION      SPCLON
C           DOUBLE PRECISION      SPCRAD
C           DOUBLE PRECISION      SPGALT
C           DOUBLE PRECISION      SPGLAT
C           DOUBLE PRECISION      SPGLON
C           DOUBLE PRECISION      SPOINT ( 3 )
C           DOUBLE PRECISION      SRFVEC ( 3 )
C           DOUBLE PRECISION      SUNLT
C           DOUBLE PRECISION      SUNPOS ( 3 )
C           DOUBLE PRECISION      SUNST  ( 6 )
C           DOUBLE PRECISION      SUPCLN
C           DOUBLE PRECISION      SUPCLT
C           DOUBLE PRECISION      SUPCRD
C           DOUBLE PRECISION      SUPGAL
C           DOUBLE PRECISION      SUPGLN
C           DOUBLE PRECISION      SUPGLT
C           DOUBLE PRECISION      TRGEPC
C
C           INTEGER               I
C           INTEGER               N
C     C
C     C     Saved variables
C     C
C           SAVE                  METHOD
C     C
C     C     Initial values
C     C
C           DATA                  METHOD / 'Intercept/ellipsoid',
C          .                               'Near point/ellipsoid',
C          .                      'Intercept/DSK/Unprioritized',
C          .                      'Nadir/DSK/Unprioritized'      /
C     C
C     C     Load kernel files via the meta-kernel.
C     C
C           CALL FURNSH ( META )
C
C     C
C     C     Convert the UTC request time to ET (seconds past
C     C     J2000, TDB).
C     C
C           CALL STR2ET ( '2008 AUG 11 00:00:00', ET )
C
C     C
C     C     Look up the target body's radii. We'll use these to
C     C     convert Cartesian to planetographic coordinates. Use
C     C     the radii to compute the flattening coefficient of
C     C     the reference ellipsoid.
C     C
C           CALL BODVRD ( 'MARS', 'RADII', 3, N, RADII )
C
C     C
C     C     Let RE and RP be, respectively, the equatorial and
C     C     polar radii of the target.
C     C
C           RE = RADII( 1 )
C           RP = RADII( 3 )
C
C           F  = ( RE - RP ) / RE
C
C     C
C     C     Compute the sub-solar point using light time and stellar
C     C     aberration corrections. Use the "target surface intercept"
C     C     definition of sub-solar point on the first loop
C     C     iteration, and use the "near point" definition on the
C     C     second.
C     C
C           DO I = 1, NMETH
C
C              CALL SUBSLR ( METHOD(I),
C          .                'MARS',  ET,     'IAU_MARS', 'CN+S',
C          .                'EARTH', SPOINT, TRGEPC,     SRFVEC )
C     C
C     C        Convert the sub-solar point's rectangular coordinates
C     C        to planetographic longitude, latitude and altitude.
C     C        Convert radians to degrees.
C     C
C              CALL RECPGR ( 'MARS', SPOINT, RE,    F,
C          .                 SPGLON, SPGLAT, SPGALT   )
C
C              SPGLON = SPGLON * DPR ()
C              SPGLAT = SPGLAT * DPR ()
C
C     C
C     C        Convert sub-solar point's rectangular coordinates to
C     C        planetocentric radius, longitude, and latitude. Convert
C     C        radians to degrees.
C     C
C              CALL RECLAT ( SPOINT, SPCRAD, SPCLON, SPCLAT )
C
C              SPCLON = SPCLON * DPR ()
C              SPCLAT = SPCLAT * DPR ()
C
C     C
C     C        Compute the Sun's apparent position relative to the
C     C        sub-solar point at TRGEPC. Add the position of
C     C        the sub-solar point relative to the target's center
C     C        to obtain the position of the sun relative to the
C     C        target's center. Express the latter position in
C     C        planetographic coordinates.
C     C
C              CALL SPKCPO ( 'SUN',  TRGEPC, 'IAU_MARS', 'OBSERVER',
C          .                 'CN+S', SPOINT, 'MARS',     'IAU_MARS',
C          .                 SUNST,  SUNLT                            )
C
C              CALL VADD ( SUNST, SPOINT, SUNPOS )
C
C              CALL RECPGR ( 'MARS', SUNPOS, RE,    F,
C          .                 SUPGLN, SUPGLT, SUPGAL   )
C
C              SUPGLN = SUPGLN * DPR ()
C              SUPGLT = SUPGLT * DPR ()
C
C     C
C     C        Convert the Sun's rectangular coordinates to
C     C        planetocentric radius, longitude, and latitude.
C     C        Convert radians to degrees.
C     C
C              CALL RECLAT ( SUNPOS, SUPCRD, SUPCLN, SUPCLT )
C
C              SUPCLN = SUPCLN * DPR ()
C              SUPCLT = SUPCLT * DPR ()
C
C     C
C     C        Write the results.
C     C
C              WRITE(*,FM) ' '
C              WRITE(*,* ) 'Computation method = ', METHOD(I)
C              WRITE(*,FM) ' '
C              WRITE(*,FM)
C          .   '  Sub-solar point altitude            (km) = ', SPGALT
C              WRITE(*,FM)
C          .   '  Sub-solar planetographic longitude (deg) = ', SPGLON
C              WRITE(*,FM)
C          .   '  Sun''s planetographic longitude     (deg) = ', SUPGLN
C              WRITE(*,FM)
C          .   '  Sub-solar planetographic latitude  (deg) = ', SPGLAT
C              WRITE(*,FM)
C          .   '  Sun''s planetographic latitude      (deg) = ', SUPGLT
C              WRITE(*,FM)
C          .   '  Sub-solar planetocentric longitude (deg) = ', SPCLON
C              WRITE(*,FM)
C          .   '  Sun''s planetocentric longitude     (deg) = ', SUPCLN
C              WRITE(*,FM)
C          .   '  Sub-solar planetocentric latitude  (deg) = ', SPCLAT
C              WRITE(*,FM)
C          .   '  Sun''s planetocentric latitude      (deg) = ', SUPCLT
C              WRITE(*,FM) ' '
C
C           END DO
C
C           END
C
C
C
C     When this program was executed on a PC/Linux/gfortran 64-bit
C     platform, the output was:
C
C
C      Computation method = Intercept/ellipsoid
C
C       Sub-solar point altitude            (km) =           0.000000000
C       Sub-solar planetographic longitude (deg) =         175.810675508
C       Sun's planetographic longitude     (deg) =         175.810675508
C       Sub-solar planetographic latitude  (deg) =          23.668550281
C       Sun's planetographic latitude      (deg) =          23.420823362
C       Sub-solar planetocentric longitude (deg) =        -175.810675508
C       Sun's planetocentric longitude     (deg) =        -175.810675508
C       Sub-solar planetocentric latitude  (deg) =          23.420819936
C       Sun's planetocentric latitude      (deg) =          23.420819936
C
C
C      Computation method = Near point/ellipsoid
C
C       Sub-solar point altitude            (km) =          -0.000000000
C       Sub-solar planetographic longitude (deg) =         175.810675408
C       Sun's planetographic longitude     (deg) =         175.810675408
C       Sub-solar planetographic latitude  (deg) =          23.420823362
C       Sun's planetographic latitude      (deg) =          23.420823362
C       Sub-solar planetocentric longitude (deg) =        -175.810675408
C       Sun's planetocentric longitude     (deg) =        -175.810675408
C       Sub-solar planetocentric latitude  (deg) =          23.175085578
C       Sun's planetocentric latitude      (deg) =          23.420819936
C
C
C      Computation method = Intercept/DSK/Unprioritized
C
C       Sub-solar point altitude            (km) =          -4.052254284
C       Sub-solar planetographic longitude (deg) =         175.810675512
C       Sun's planetographic longitude     (deg) =         175.810675512
C       Sub-solar planetographic latitude  (deg) =          23.668848891
C       Sun's planetographic latitude      (deg) =          23.420823362
C       Sub-solar planetocentric longitude (deg) =        -175.810675512
C       Sun's planetocentric longitude     (deg) =        -175.810675512
C       Sub-solar planetocentric latitude  (deg) =          23.420819936
C       Sun's planetocentric latitude      (deg) =          23.420819936
C
C
C      Computation method = Nadir/DSK/Unprioritized
C
C       Sub-solar point altitude            (km) =          -4.022302438
C       Sub-solar planetographic longitude (deg) =         175.810675412
C       Sun's planetographic longitude     (deg) =         175.810675412
C       Sub-solar planetographic latitude  (deg) =          23.420823362
C       Sun's planetographic latitude      (deg) =          23.420823362
C       Sub-solar planetocentric longitude (deg) =        -175.810675412
C       Sun's planetocentric longitude     (deg) =        -175.810675412
C       Sub-solar planetocentric latitude  (deg) =          23.174793924
C       Sun's planetocentric latitude      (deg) =          23.420819936
C
C
C$ Restrictions
C
C    None.
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
C        Added FAILED tests.
C
C       14-JUL-2016 (NJB)
C
C        Now uses surface mapping tracking capability.
C        Updated header.
C
C       09-FEB-2015 (NJB) 
C 
C        Support for surface specification was added. 
C        Header was updated to document DSK features.
C      
C       24-DEC-2014 (NJB) 
C
C        Updated to support DSK data.
C
C-    SPICELIB Version 2.0.0, 31-MAR-2014 (NJB) (SCK) (BVS)
C
C        Bug fix: stellar aberration is no longer applied to the
C        observer-to-estimated sub-solar point vector while solving for
C        the sub-solar point. This correction involved unnecessary code
C        but did not affect this routine's outputs.
C
C        Bug fix: FIRST is now set to .FALSE. at the completion of a
C        successful initialization pass. This does not affect the
C        routine's outputs but improves efficiency.
C
C        Exceptions removed: the observer and target are now 
C        permitted to coincide.
C
C        Upgrade: the algorithm for finding the apparent state of the
C        sun as seen from the estimated sub-solar point has been
C        improved.
C
C        Upgrade: this routine now uses ZZVALCOR rather than ZZPRSCOR,
C        simplifying the implementation.
C
C        The header example program was updated to reflect the new
C        method of computing the apparent sun location, and the set
C        of kernels referenced by the example meta-kernel were updated.
C        The display of the program's output was updated accordingly.
C
C        References to the new PXFRM2 routine were added, which changed
C        the Detailed Output section.
C
C        Updated to save the input body names and ZZBODTRN state
C        counters and to do name-ID conversions only if the counters
C        have changed.
C
C        Updated to save the input frame name and POOL state counter
C        and to do frame name-ID conversion only if the counter has
C        changed.
C
C        Updated to call LJUCRS instead of CMPRSS/UCASE. 
C
C-    SPICELIB Version 1.1.0, 18-MAY-2010 (NJB) 
C
C        Bug fix: calls to FAILED() have been added after
C        SPK calls, target radius lookup, near point
C        and surface intercept computations.
C
C-    SPICELIB Version 1.0.1, 17-MAR-2009 (NJB) 
C
C        Typo correction: changed FIXFRM to FIXREF in header
C        documentation. Meta-kernel name suffix was changed to
C        ".tm" in header code example.
C
C        Typo correction in Required_Reading, changed 
C        FRAME to FRAMES.
C
C-    SPICELIB Version 1.0.0, 02-MAR-2008 (NJB) 
C
C-&
 
C$ Index_Entries
C
C     find sub-solar point on target body
C     find nearest point to sun on target body
C
C-&
 

C$ Revisions
C
C     None.
C
C-&

 
C
C     SPICELIB functions
C
      DOUBLE PRECISION      CLIGHT
      DOUBLE PRECISION      TOUCHD
      DOUBLE PRECISION      VDIST
      
      LOGICAL               EQSTR
      LOGICAL               FAILED
      LOGICAL               RETURN
 
C
C     Local parameters
C
      CHARACTER*(*)         RNAME
      PARAMETER           ( RNAME  =  'SUBSLR' )

C
C     This value will become system-dependent when systems
C     using 128-bit d.p. numbers are supported by SPICELIB.
C     CNVLIM, when added to 1.0D0, should yield 1.0D0. 
C
      DOUBLE PRECISION      CNVLIM
      PARAMETER           ( CNVLIM = 1.D-17 )
     

      INTEGER               MAXITR
      PARAMETER           ( MAXITR = 10 )

      INTEGER               SUN
      PARAMETER           ( SUN    = 10 )
 
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
 
      DOUBLE PRECISION      ALTSUN
      DOUBLE PRECISION      DVEC   ( 3 )
      DOUBLE PRECISION      ETDIFF
      DOUBLE PRECISION      J2POS  ( 3 )
      DOUBLE PRECISION      LT
      DOUBLE PRECISION      LTDIFF
      DOUBLE PRECISION      OBSPOS ( 3 )
      DOUBLE PRECISION      OBSRNG
      DOUBLE PRECISION      PREVET
      DOUBLE PRECISION      PREVLT
      DOUBLE PRECISION      RADII  ( 3 )
      DOUBLE PRECISION      S
      DOUBLE PRECISION      SLT
      DOUBLE PRECISION      SPOS   ( 3 )
      DOUBLE PRECISION      SSBOST ( 6 )
      DOUBLE PRECISION      SSBTST ( 6 )
      DOUBLE PRECISION      SSLRLT
      DOUBLE PRECISION      SSLRST ( 6 )
      DOUBLE PRECISION      SUNST  ( 6 )
      DOUBLE PRECISION      TPOS   ( 3 )
      DOUBLE PRECISION      XFORM  ( 3, 3 )

      INTEGER               FIXCID
      INTEGER               FIXCLS
      INTEGER               FIXCTR
      INTEGER               FIXFID
      INTEGER               I
      INTEGER               NITR
      INTEGER               NRADII
      INTEGER               NSURF
      INTEGER               OBSCDE
      INTEGER               SHAPE
      INTEGER               SRFLST ( MAXSRF )
      INTEGER               TRGCDE
     
      LOGICAL               ATTBLK ( NABCOR )
      LOGICAL               FIRST
      LOGICAL               FND
      LOGICAL               NEAR
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
      INTEGER               SVREFC

C
C     Saved surface name/ID item declarations.
C
      INTEGER               SVCTR4 ( CTRSIZ )

C
C     Saved variables
C
      SAVE                  FIRST
      SAVE                  NEAR
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
      SAVE                  SVREFC

C
C     Saved surface name/ID items.
C
      SAVE                  SVCTR4

C
C     Initial values
C
      DATA                  FIRST   / .TRUE. /
      DATA                  NEAR    / .TRUE. /
      DATA                  PRVCOR  / ' '    /
      DATA                  PRVMTH  / ' '    /
      DATA                  SHAPE   / ELLSHP /


C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( RNAME )
  
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
C        Reject an aberration correction flag calling for transmission
C        corrections.
C
         IF ( ATTBLK(XMTIDX) ) THEN

            CALL SETMSG ( 'Aberration correction flag # calls for '
     .      //            'transmission-style corrections.'        )
            CALL ERRCH  ( '#', ABCORR                              )
            CALL SIGERR ( 'SPICE(NOTSUPPORTED)'                    )
            CALL CHKOUT ( RNAME                                    )
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
C        The above definitions are consistent with those used by
C        ZZPRSCOR.
C 
         XMIT    =  ATTBLK ( XMTIDX )
         USELT   =  ATTBLK ( LTIDX  )
         USECN   =  ATTBLK ( CNVIDX )
         USESTL  =  ATTBLK ( STLIDX )

C
C        The aberration correction flag is recognized; save it.
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
C     The target body may not be the sun.
C
      IF ( TRGCDE .EQ. SUN ) THEN
 
         CALL SETMSG ( 'The target body is the sun; the sub-solar '   
     .   //            'point is undefined for this case.'          )
         CALL SIGERR ( 'SPICE(INVALIDTARGET)'                       )
         CALL CHKOUT ( RNAME                                        )
         RETURN
 
      END IF
      
C
C     Determine the attributes of the frame designated by FIXREF.
C
      CALL ZZNAMFRM ( SVCTR3, SVFREF, SVREFC, FIXREF, FIXFID )

      CALL FRINFO ( FIXFID, FIXCTR, FIXCLS, FIXCID, FND )

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
      IF ( FIXCTR .NE. TRGCDE ) THEN

         CALL SETMSG ( 'Reference frame # is not centered at the ' 
     .   //            'target body #. The ID code of the frame '
     .   //            'center is #.'                             )
         CALL ERRCH  ( '#',  FIXREF                               )
         CALL ERRCH  ( '#',  TARGET                               )
         CALL ERRINT ( '#',  FIXCTR                               )
         CALL SIGERR ( 'SPICE(INVALIDFRAME)'                      )
         CALL CHKOUT ( RNAME                                      )
         RETURN

      END IF

C
C     Check whether the surface name/ID mapping has been updated.
C
      CALL ZZSRFTRK ( SVCTR4, SURFUP )

C
C     If necessary, parse the method specification. PRVMTH
C     and the derived variables NEAR and SHAPE start out with
C     valid values. PRVMTH records the last valid value of
C     METHOD; NEAR and SHAPE are the corresponding variables.
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


         IF ( SUBTYP .EQ. ' ' ) THEN

            CALL SETMSG ( 'Sub-solar point type is required '
     .      //            'but was not found in the method string #.' )
            CALL ERRCH  ( '#', METHOD                                 )
            CALL SIGERR ( 'SPICE(INVALIDSUBTYPE)'                     )
            CALL CHKOUT ( RNAME                                       )
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
            CALL SETMSG ( 'Returned shape value from method '
     .      //            'string was <#>.'                  )
            CALL ERRCH  ( '#', SHPSTR                        )
            CALL SIGERR ( 'SPICE(BUG)'                       )
            CALL CHKOUT ( RNAME                              )
            RETURN

         END IF


         IF ( SHAPE .EQ. ELLSHP ) THEN
C
C           Allow both "near point" and "nadir" expressions
C           the ellipsoid case, since these are equivalent.
C
            NEAR =      EQSTR( SUBTYP, 'NEAR POINT' )
     .             .OR. EQSTR( SUBTYP, 'NADIR'      )

         ELSE
C
C           "near point" is not supported for DSKs.
C
            NEAR = EQSTR( SUBTYP, 'NADIR' )

         END IF


         IF ( .NOT. NEAR ) THEN

            IF ( .NOT. EQSTR( SUBTYP, 'INTERCEPT' ) ) THEN

               CALL SETMSG ( 'Invalid sub-solar point type <#> '
     .         //            'was found in the method string #.' )
               CALL ERRCH  ( '#', SUBTYP                         )
               CALL ERRCH  ( '#', METHOD                         )
               CALL SIGERR ( 'SPICE(INVALIDSUBTYPE)'             )
               CALL CHKOUT ( RNAME                               )
               RETURN

            END IF

         END IF

C
C        Save the current value of METHOD. 
C
         PRVMTH = METHOD

      END IF

C
C     At this point, the first pass actions were successful.
C
      FIRST = .FALSE.


      IF ( SHAPE .EQ. DSKSHP ) THEN
C
C        This is the DSK case.
C
C        Initialize the intercept algorithm to use a DSK
C        model for the surface of the target body. 
C        
         CALL ZZSUDSKI ( TRGCDE, NSURF, SRFLST, FIXFID )


      ELSE IF ( SHAPE .NE. ELLSHP ) THEN

         CALL SETMSG ( 'Computation method argument was <#>; this ' 
     .   //            'string must specify a supported shape '
     .   //            'model and computation type. See the '
     .   //            'header of SUBSLR for details.'            )
         CALL ERRCH  ( '#',  METHOD                               )
         CALL SIGERR ( 'SPICE(INVALIDMETHOD)'                     )
         CALL CHKOUT ( RNAME                                      )
         RETURN

      END IF

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF


C
C     Get the sign S prefixing LT in the expression for TRGEPC.
C     When light time correction is not used, setting S = 0
C     allows us to seamlessly set TRGEPC equal to ET.
C
      IF ( USELT ) THEN
         S = -1.D0         
      ELSE
         S =  0.D0
      END IF
 
C
C     Determine the position of the observer in target body-fixed
C     coordinates. This is a first estimate.
C
C         -  Call SPKEZP to compute the position of the target body as
C            seen from the observing body and the light time (LT)
C            between them. We request that the coordinates of POS be
C            returned relative to the body fixed reference frame
C            associated with the target body, using aberration
C            corrections specified by the input argument ABCORR.
C
C         -  Call VMINUS to negate the direction of the vector (OBSPOS)
C            so it will be the position of the observer as seen from
C            the target body in target body fixed coordinates.
C
C            Note that this result is not the same as the result of
C            calling SPKEZP with the target and observer switched. We
C            computed the vector FROM the observer TO the target in
C            order to get the proper light time and stellar aberration
C            corrections (if requested). Now we need the inverse of
C            that corrected vector in order to compute the sub-solar
C            point.
C
      CALL SPKEZP ( TRGCDE, ET, FIXREF, ABCORR, OBSCDE, TPOS, LT )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

C
C     Negate the target's position to obtain the position of the
C     observer relative to the target.
C
      CALL VMINUS ( TPOS, OBSPOS )

C
C     Make a first estimate of the target epoch.
C
      TRGEPC = ET + S*LT

C
C     Find the sub-solar point given the target epoch, observer-target
C     position, and target body orientation we've already computed. If
C     we're not using light time correction, this is all we need do.
C     Otherwise, our result will give us an initial estimate of the
C     target epoch, which we'll then improve.
C
C     Get the radii of the target body from the kernel pool.
C
      CALL BODVCD ( TRGCDE, 'RADII', 3, NRADII, RADII )

C
C     Get the position of the Sun SPOS as seen from the target
C     in the target body-fixed frame at TRGEPC. 
C
      CALL SPKEZP ( SUN, TRGEPC, FIXREF, ABCORR, TRGCDE, SPOS, SLT )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

C
C     Make a first estimate of the sub-solar point. The algorithm
C     we use depends on the sub-solar point definition.
C
      IF ( NEAR ) THEN
C
C        Locate the nearest point to the Sun on the target's
C        reference ellipsoid.
C
         CALL NEARPT ( SPOS,   RADII(1), RADII(2), RADII(3),
     .                 SPOINT, ALTSUN                       )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        If the target is an ellipsoid, the NEARPT call above does the
C        trick. For DSKs, we define a ray emanating from the sun and
C        passing through the near point on the reference ellipsoid. The
C        closest ray-DSK surface intercept to the sun is the initial
C        estimate of the sub-solar point.
C
         IF ( SHAPE .EQ. DSKSHP ) THEN           
C
C           Generate the ray direction; find the DSK intercept.
C
            CALL VSUB ( SPOINT, SPOS, DVEC )

            CALL ZZSBFXR ( TRGCDE, NSURF, SRFLST, TRGEPC,    
     .                     FIXFID, SPOS,  DVEC,   SPOINT, FND )
 
            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF

            IF ( .NOT. FND ) THEN

               CALL SETMSG ( 'No sub-solar point was found on '
     .         //            'the surface defined by DSK data. '
     .         //            'Observer is #; target is #. This '
     .         //            'problem can occur for bodies having '
     .         //            'shapes not well modeled by '
     .         //            'ellipsoids. Consider using the '
     .         //            '"Intercept: DSK" computation method.' )
               CALL ERRCH  ( '#', OBSRVR                            )
               CALL ERRCH  ( '#', TARGET                            )
               CALL SIGERR ( 'SPICE(SUBPOINTNOTFOUND)'              )
               CALL CHKOUT ( RNAME                                  )
               RETURN

            END IF

         END IF


      ELSE
C
C        This is the case for the "intercept" sub-solar point
C        definition.
C
C        Generate the ray direction.
C
         CALL VMINUS ( SPOS, DVEC )

         IF ( SHAPE .EQ. ELLSHP ) THEN
C
C           Locate the surface intercept of the ray from the 
C           sun to the target center.
C
            CALL SURFPT ( SPOS,   DVEC, RADII(1), RADII(2), RADII(3), 
     .                    SPOINT, FND                                )
               
            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF

            IF ( .NOT. FND ) THEN
C
C              If there's no intercept, we have a numerical problem.
C              
               CALL SETMSG ( 'No intercept of sun-target '
     .         //            'ray was found.'              )
               CALL SIGERR ( 'SPICE(DEGENERATECASE)'       )
               CALL CHKOUT ( RNAME                         )
               RETURN

            END IF


         ELSE
C
C           Find the DSK intercept.
C
            CALL ZZSBFXR ( TRGCDE, NSURF, SRFLST, TRGEPC,    
     .                     FIXFID, SPOS,  DVEC,   SPOINT, FND )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF

            IF ( .NOT. FND ) THEN

               CALL SETMSG ( 'No sub-solar point was found on '
     .         //            'the surface defined by DSK data. '
     .         //            'Observer is #; target is #. This '
     .         //            'problem can occur for a body having '
     .         //            'an irregular shape such that the '
     .         //            'origin of the body-fixed reference '
     .         //            'frame is outside of the body. A tor'
     .         //            'us is an example of such a shape.'   )
               CALL ERRCH  ( '#', OBSRVR                           )
               CALL ERRCH  ( '#', TARGET                           )
               CALL SIGERR ( 'SPICE(SUBPOINTNOTFOUND)'             )
               CALL CHKOUT ( RNAME                                 )
               RETURN

            END IF

         END IF

      END IF

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

C
C     If we're not using light time and stellar aberration
C     corrections, we're almost done now. Note that we need only
C     check for use of light time corrections, because use of
C     stellar aberration corrections alone has been prevented by an
C     earlier check.
C        
      IF ( .NOT. USELT ) THEN
C
C        The TRGEPC value we'll return is simply the input time. The
C        previous call to SPKEZP call yielded the vector OBSPOS. SPOINT
C        was set immediately above. The only output left to compute is
C        SRFVEC.
C
         TRGEPC = ET

         CALL VSUB ( SPOINT, OBSPOS, SRFVEC )

         CALL CHKOUT ( RNAME )
         RETURN

      END IF

C
C     Compute the range from the observer to the sub-solar
C     point. We'll use this range for a light time estimate.
C
      OBSRNG = VDIST ( OBSPOS, SPOINT )

C
C     Compute the one-way light time and target epoch based on our
C     first computation of SPOINT. The coefficient S has been
C     set to give us the correct answer for each aberration
C     correction case.
C
      LT     = OBSRNG / CLIGHT()
      TRGEPC = ET  + S*LT

C
C     We'll now make an improved sub-solar point estimate using the
C     previous estimate of the sub-solar point. The number of
C     iterations depends on the light time correction type.
C
      IF ( USECN ) THEN
         NITR = MAXITR
      ELSE
         NITR = 1
      END IF

C
C     Get the J2000-relative state of the observer relative to
C     the solar system barycenter at ET.
C
      CALL SPKSSB ( OBSCDE, ET, 'J2000', SSBOST )

      IF ( FAILED() ) THEN
         CALL CHKOUT ( RNAME )
         RETURN
      END IF

C
C     Initialize the variables required to evaluate the 
C     loop termination condition.
C
      I      = 0
      LTDIFF = 1.D0
      ETDIFF = 1.D0
      PREVLT = LT
      PREVET = TRGEPC

      DO WHILE (       ( I      .LT.   NITR                ) 
     .           .AND. ( LTDIFF .GT. ( CNVLIM * ABS(LT) )  )
     .           .AND. ( ETDIFF .GT.   0.D0                )  )
C
C        Get the J2000-relative state of the target relative to
C        the solar system barycenter at the target epoch.
C
         CALL SPKSSB ( TRGCDE, TRGEPC, 'J2000', SSBTST )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        Find the position of the observer relative to the target.
C        Convert this vector from the J2000 frame to the target
C        frame at TRGEPC.
C
         CALL VSUB   ( SSBOST,  SSBTST, J2POS         )
         CALL PXFORM ( 'J2000', FIXREF, TRGEPC, XFORM )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         CALL MXV ( XFORM, J2POS, OBSPOS )

C
C        Note: if we're using stellar aberration correction, we do not
C        apply the stellar aberration correction of the estimated
C        sub-solar point as seen by the observer to the next estimate
C        of the sun-solar point. The location of this point depends on
C        the illumination of the target, which is a function of the
C        observer-surface point light time. This is the only way in
C        which the observer plays a role in determining the sub-solar
C        point.
C
C        Stellar aberration of the sun's position relative to the
C        sub-solar point *is* used.
C
C        First find the apparent position of the sun as seen from the
C        estimated sub-solar point.
C 

         CALL SPKCPO ( 'SUN',  TRGEPC, FIXREF, 'OBSERVER',
     .                 ABCORR, SPOINT, TARGET, FIXREF, SUNST, SLT )
 
C
C        Create the target-center to sun vector.
C
         CALL VADD ( SUNST, SPOINT, SPOS )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        Find the sub-solar point using the current estimated
C        geometry.
C
         IF ( NEAR ) THEN
C
C           Locate the nearest point to the sun on the target's
C           reference ellipsoid.
C
            CALL NEARPT ( SPOS,   RADII(1), RADII(2), RADII(3),
     .                    SPOINT, ALTSUN                        )

            IF ( FAILED() ) THEN
               CALL CHKOUT ( RNAME )
               RETURN
            END IF
C
C           If the target is an ellipsoid, the NEARPT call above does
C           the trick. For DSKs, we define a ray emanating from the sun
C           and passing through the near point on the reference
C           ellipsoid. The closest ray-DSK surface intercept to the sun
C           is the initial estimate of the sub-solar point.
C
            IF ( SHAPE .EQ. DSKSHP ) THEN
C
C              Locate the surface intercept of the ray from the 
C              Sun to the target center.
C
               CALL VSUB ( SPOINT, SPOS, DVEC )

               CALL ZZSBFXR ( TRGCDE, NSURF, SRFLST, TRGEPC,    
     .                        FIXFID, SPOS,  DVEC,   SPOINT, FND )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

               IF ( .NOT. FND ) THEN

                  CALL SETMSG ( 'No sub-solar point was found on '
     .            //            'the surface defined by DSK data. '
     .            //            'Observer is #; target is #. This '
     .            //            'problem can occur for bodies hav'
     .            //            'ing shapes not well modeled by '
     .            //            'ellipsoids. Consider using the '
     .            //            '"Intercept: DSK" computation '
     .            //            'method.'                        )
                  CALL ERRCH  ( '#', OBSRVR                      )
                  CALL ERRCH  ( '#', TARGET                      )
                  CALL SIGERR ( 'SPICE(SUBPOINTNOTFOUND)'        )
                  CALL CHKOUT ( RNAME                            )
                  RETURN

               END IF

            END IF


         ELSE

C           This is the "intercept" case.
C
C           Generate the ray direction.
C
            CALL VMINUS ( SPOS, DVEC )


            IF ( SHAPE .EQ. ELLSHP ) THEN
C
C              Locate the surface intercept of the ray from the 
C              sun to the target center.
C
               CALL SURFPT ( SPOS,     DVEC,   RADII(1), RADII(2), 
     .                       RADII(3), SPOINT, FND                )
               
               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

               IF ( .NOT. FND ) THEN
C
C                 If there's no intercept, we have a numerical problem.
C              
                  CALL SETMSG ( 'No intercept of sun-target '
     .            //            'ray was found.'              )
                  CALL SIGERR ( 'SPICE(DEGENERATECASE)'       )
                  CALL CHKOUT ( RNAME                         )
                  RETURN

               END IF

            ELSE
C
C              Find the DSK intercept.
C
               CALL ZZSBFXR ( TRGCDE, NSURF, SRFLST, TRGEPC,    
     .                        FIXFID, SPOS,  DVEC,   SPOINT, FND )

               IF ( FAILED() ) THEN
                  CALL CHKOUT ( RNAME )
                  RETURN
               END IF

               IF ( .NOT. FND ) THEN

                  CALL SETMSG ( 'No sub-solar point was found on '
     .            //            'the surface defined by DSK data. '
     .            //            'Observer is #; target is #. This '
     .            //            'problem can occur for a body having '
     .            //            'an irregular shape such that the '
     .            //            'origin of the body-fixed reference '
     .            //            'frame is outside of the body. A tor'
     .            //            'us is an example of such a shape.'   )
                  CALL ERRCH  ( '#', OBSRVR                           )
                  CALL ERRCH  ( '#', TARGET                           )
                  CALL SIGERR ( 'SPICE(SUBPOINTNOTFOUND)'             )
                  CALL CHKOUT ( RNAME                                 )
                  RETURN

               END IF

            END IF

         END IF

C
C        Update the observer to sub-solar point range.
C
         OBSRNG = VDIST ( OBSPOS, SPOINT )

C
C        Compute a new light time estimate and new target epoch.
C
         LT      =  OBSRNG / CLIGHT()
         TRGEPC  =  ET  + S*LT

C
C        At this point, we have new estimates of the sub-solar point
C        SPOINT, the observer altitude OBSRNG, the target epoch
C        TRGEPC, and the position of the observer relative to the
C        target OBSPOS.
C
C        We use the d.p. identity function TOUCHD to force the compiler
C        to create double precision arguments from the differences
C        LT-PREVLT and TRGEPC-PREVET. Some compilers will perform
C        extended-precision register arithmetic, which can prevent a
C        difference from rounding to zero. Simply storing the result of
C        the subtraction in a double precision variable doesn't solve
C        the problem, because that variable can be optimized out of
C        existence.
C
         LTDIFF  =   ABS( TOUCHD(LT     - PREVLT) )
         ETDIFF  =   ABS( TOUCHD(TRGEPC - PREVET) )
         PREVLT  =   LT
         PREVET  =   TRGEPC        
         I       =   I + 1

         CALL SPKSSB ( TRGCDE, TRGEPC, 'J2000', SSBTST )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

C
C        Find the position of the observer relative to the target.
C        Convert this vector from the J2000 frame to the target
C        frame at TRGEPC.
C
         CALL VSUB   ( SSBOST,  SSBTST, J2POS         )
         CALL PXFORM ( 'J2000', FIXREF, TRGEPC, XFORM )

         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         CALL MXV ( XFORM, J2POS, OBSPOS )

      END DO

C
C     SPOINT and TRGEPC have been set. Compute SRFVEC, using
C     both light time and stellar aberration corrections if
C     these have been requested by the caller. Note that
C     the position of the target body was computed without
C     stellar aberration corrections, so we may not be able
C     to use it for this computation.
C     
      IF ( USESTL ) THEN

         CALL SPKCPT ( SPOINT,   TARGET, FIXREF, ET,     FIXREF,
     .                 'TARGET', ABCORR, OBSRVR, SSLRST, SSLRLT )
         IF ( FAILED() ) THEN
            CALL CHKOUT ( RNAME )
            RETURN
         END IF

         CALL VEQU ( SSLRST, SRFVEC )

      ELSE

         CALL VSUB ( SPOINT, OBSPOS, SRFVEC )

      END IF

      CALL CHKOUT ( RNAME )
      RETURN
      END
