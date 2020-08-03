C$Procedure ILUMIN ( Illumination angles )
 
      SUBROUTINE ILUMIN ( METHOD, TARGET, ET,     FIXREF,  
     .                    ABCORR, OBSRVR, SPOINT, TRGEPC,  
     .                    SRFVEC, PHASE,  SOLAR,  EMISSN  )
 
C$ Abstract
C
C     Find the illumination angles (phase, solar incidence, and
C     emission) at a specified surface point of a target body.
C
C     This routine supersedes ILLUM.
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
C     ANGLES
C     GEOMETRY
C     ILLUMINATION
C
C$ Declarations
 
      IMPLICIT NONE

      CHARACTER*(*)         METHOD
      CHARACTER*(*)         TARGET
      DOUBLE PRECISION      ET
      CHARACTER*(*)         FIXREF
      CHARACTER*(*)         ABCORR
      CHARACTER*(*)         OBSRVR
      DOUBLE PRECISION      SPOINT ( 3 )
      DOUBLE PRECISION      TRGEPC
      DOUBLE PRECISION      SRFVEC ( 3 )
      DOUBLE PRECISION      PHASE
      DOUBLE PRECISION      SOLAR
      DOUBLE PRECISION      EMISSN
 
C$ Brief_I/O
C
C     Variable  I/O  Description
C     --------  ---  --------------------------------------------------
C     METHOD     I   Computation method.
C     TARGET     I   Name of target body.
C     ET         I   Epoch in ephemeris seconds past J2000 TDB.
C     FIXREF     I   Body-fixed, body-centered target body frame.
C     ABCORR     I   Desired aberration correction.
C     OBSRVR     I   Name of observing body.
C     SPOINT     I   Body-fixed coordinates of a target surface point.
C     TRGEPC     O   Target surface point epoch.
C     SRFVEC     O   Vector from observer to target surface point.
C     PHASE      O   Phase angle at the surface point.
C     SOLAR      O   Solar incidence angle at the surface point.
C     EMISSN     O   Emission angle at the surface point.
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
C                       The illumination angle computation uses a
C                       triaxial ellipsoid to model the surface of the
C                       target body. The ellipsoid's radii must be
C                       available in the kernel pool.
C
C
C                    'DSK/UNPRIORITIZED[/SURFACES = <surface list>]'
C
C                       The illumination angle computation uses
C                       topographic data to model the surface of the
C                       target body. These data must be provided by
C                       loaded DSK files.
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
C     ET          is the epoch, expressed as seconds past J2000 TDB,
C                 for which the apparent illumination angles at the
C                 specified surface point on the target body, as seen
C                 from the observing body, are to be computed.
C
C
C     FIXREF      is the name of the body-fixed, body-centered
C                 reference frame associated with the target body. The
C                 input surface point SPOINT and the output vector
C                 SRFVEC are expressed relative to this reference
C                 frame. The string FIXREF is case-insensitive, and
C                 leading and trailing blanks in FIXREF are not
C                 significant.
C
C
C     ABCORR      is the aberration correction to be used in computing
C                 the position and orientation of the target body and
C                 the location of the Sun.
C         
C                 For remote sensing applications, where the apparent
C                 illumination angles seen by the observer are desired,
C                 normally either of the corrections 
C              
C                    'LT+S' 
C                    'CN+S'
C     
C                 should be used. These and the other supported options
C                 are described below. ABCORR may be any of the 
C                 following:
C
C                    'NONE'     No aberration correction.
C
C                 Let LT represent the one-way light time between the
C                 observer and SPOINT (note: NOT between the observer
C                 and the target body's center). The following values
C                 of ABCORR apply to the "reception" case in which
C                 photons depart from SPOINT at the light-time
C                 corrected epoch ET-LT and *arrive* at the observer's
C                 location at ET:
C
C                    'LT'       Correct both the position of SPOINT as
C                               seen by the observer, and the position
C                               of the Sun as seen by the target, for
C                               light time. Correct the orientation of
C                               the target for light time.
C
C                    'LT+S'     Correct both the position of SPOINT as
C                               seen by the observer, and the position
C                               of the Sun as seen by the target, for
C                               light time and stellar aberration.
C                               Correct the orientation of the target
C                               for light time.
C
C                    'CN'       Converged Newtonian light time
C                               correction. In solving the light time
C                               equations for target and the Sun, the
C                               "CN" correction iterates until the
C                               solution converges.
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
C
C                 The following values of ABCORR apply to the
C                 "transmission" case in which photons *arrive* at
C                 SPOINT at the light-time corrected epoch ET+LT and
C                 *depart* from the observer's location at ET:
C
C                    'XLT'      "Transmission" case: correct for
C                               one-way light time using a Newtonian
C                               formulation. This correction yields the
C                               illumination angles at the moment that
C                               SPOINT receives photons emitted from the
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
C                               angles obtained with the 'XLT' option
C                               to account for the observer's and
C                               target's velocities relative to the
C                               solar system barycenter (the latter
C                               velocity is used in computing the
C                               direction to the apparent illumination
C                               source).
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
C                 OBSRVR may be not be identical to TARGET.
C
C
C     SPOINT      is a surface point on the target body, expressed in
C                 Cartesian coordinates, relative to the body-fixed
C                 target frame designated by FIXREF.
C
C                 SPOINT need not be visible from the observer's
C                 location at the epoch ET.
C
C                 The components of SPOINT have units of km.
C
C
C$ Detailed_Output
C
C
C     TRGEPC      is the "surface point epoch." TRGEPC is defined as
C                 follows: letting LT be the one-way light time between
C                 the observer and the input surface point SPOINT,
C                 TRGEPC is either the epoch ET-LT or ET depending on
C                 whether the requested aberration correction is,
C                 respectively, for received radiation or omitted. LT
C                 is computed using the method indicated by ABCORR.
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
C     PHASE       is the phase angle at SPOINT, as seen from OBSRVR at
C                 time ET. This is the angle between the negative of
C                 the vector SRFVEC and the SPOINT-Sun vector at
C                 TRGEPC. Units are radians. The range of PHASE is
C                 [0, pi].  
C
C     SOLAR       is the solar incidence angle at SPOINT, as seen from
C                 OBSRVR at time ET. This is the angle between the
C                 surface normal vector at SPOINT and the SPOINT-Sun
C                 vector at TRGEPC. Units are radians. The range of
C                 SOLAR is [0, pi].  
C
C     EMISSN      is the emission angle at SPOINT, as seen from OBSRVR
C                 at time ET. This is the angle between the surface
C                 normal vector at SPOINT and the negative of the
C                 vector SRFVEC. Units are radians. The range of EMISSN
C                 is [0, pi]. 
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
C     2)  If either the target or observer input strings cannot be
C         converted to an integer ID code, the error will be signaled
C         by a routine in the call tree of this routine.
C
C     3)  If OBSRVR and TARGET map to the same NAIF integer ID code,
C         the error will be signaled by a routine in the call tree of
C         this routine.
C
C     4)  If the input target body-fixed frame FIXREF is not
C         recognized, the error will be signaled by a routine in the
C         call tree of this routine. A frame name may fail to be
C         recognized because a required frame specification kernel has
C         not been loaded; another cause is a misspelling of the frame
C         name.
C
C     5)  If the input frame FIXREF is not centered at the target body,
C         the error will be signaled by a routine in the call tree of
C         this routine.
C
C     6)  If the input argument METHOD is not recognized, the error
C         will be signaled by a routine in the call tree of this
C         routine.
C
C     7)  If insufficient ephemeris data have been loaded prior to
C         calling ILUMIN, the error will be diagnosed and signaled by a
C         routine in the call tree of this routine. Note that when
C         light time correction is used, sufficient ephemeris data must
C         be available to propagate the states of observer, target, and
C         the Sun to the solar system barycenter.
C
C     8)  If the computation method specifies an ellipsoidal target
C         shape and triaxial radii of the target body have not been
C         loaded into the kernel pool prior to calling ILUMIN, the
C         error will be diagnosed and signaled by a routine in the call
C         tree of this routine.
C
C     9)  The target must be an extended body: if any of the radii of
C         the target body are non-positive, the error will be
C         diagnosed and signaled by routines in the call tree of this
C         routine.
C
C     10) If PCK data specifying the target body-fixed frame
C         orientation have not been loaded prior to calling ILUMIN,
C         the error will be diagnosed and signaled by a routine in the
C         call tree of this routine.
C
C     11) If METHOD specifies that the target surface is represented by
C         DSK data, and no DSK files are loaded for the specified
C         target, the error is signaled by a routine in the call tree
C         of this routine.
C         
C     12) If METHOD specifies that the target surface is represented
C         by DSK data, and data representing the portion of the surface
C         on which SPOINT is located are not available, an error will 
C         be signaled by a routine in the call tree of this routine.
C
C     13) If METHOD specifies that the target surface is represented
C         by DSK data, SPOINT must lie on the target surface, not above
C         or below it. A small tolerance is used to allow for round-off
C         error in the calculation determining whether SPOINT is on the
C         surface. If, in the DSK case, SPOINT is too far from the
C         surface, an error will be signaled by a routine in the call
C         tree of this routine.
C
C         If the surface is represented by a triaxial ellipsoid, SPOINT
C         is not required to be close to the ellipsoid; however, the
C         results computed by this routine will be unreliable if SPOINT
C         is too far from the ellipsoid.
C
C$ Files
C
C
C     Appropriate kernels must be loaded by the calling program before
C     this routine is called.
C
C     The following data are required:
C
C        - SPK data: ephemeris data for target, observer, and the
C          illumination source must be loaded. If aberration
C          corrections are used, the states of target, observer, and
C          the illumination source relative to the solar system
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
C
C$ Particulars
C
C
C     SPICELIB contains four routines that compute illumination angles:
C     
C        ILLUMF (same as ILLUMG, except that illumination
C                and visibility flags are returned.)
C
C        ILLUMG (same as ILUMIN, except that the caller
C                specifies the illumination source.)
C
C        ILUMIN (this routine)
C                
C        ILLUM  (deprecated)
C     
C     ILLUMF is the most capable of the set.
C
C
C     Illumination angles
C     ===================
C
C     The term "illumination angles" refers to the following set of
C     angles:
C
C
C        phase angle              Angle between the vectors from the
C                                 surface point to the observer and
C                                 from the surface point to the Sun.
C
C        solar incidence angle    Angle between the surface normal at
C                                 the specified surface point and the
C                                 vector from the surface point to the
C                                 Sun.
C
C        emission angle           Angle between the surface normal at
C                                 the specified surface point and the
C                                 vector from the surface point to the
C                                 observer.
C 
C     The diagram below illustrates the geometric relationships
C     defining these angles. The labels for the solar incidence,
C     emission, and phase angles are "s.i.", "e.", and "phase".
C
C
C                                                      *
C                                                     Sun
C
C                    surface normal vector
C                              ._                 _.
C                              |\                 /|  Sun vector
C                                \    phase      /
C                                 \   .    .    /
C                                 .            .
C                                   \   ___   /
C                              .     \/     \/
C                                    _\ s.i./
C                             .    /   \   /
C                             .   |  e. \ /
C         *             <--------------- *  surface point on
C      viewing            vector            target body
C      location           to viewing
C      (observer)         location
C
C
C     Note that if the target-observer vector, the target normal vector
C     at the surface point, and the target-sun vector are coplanar,
C     then phase is the sum of incidence and emission. This is rarely
C     true; usually
C
C        phase angle  <  solar incidence angle + emission angle
C
C     All of the above angles can be computed using light time
C     corrections, light time and stellar aberration corrections, or
C     no aberration corrections. In order to describe apparent
C     geometry as observed by a remote sensing instrument, both
C     light time and stellar aberration corrections should be used.
C     
C     The way aberration corrections are applied by this routine
C     is described below.
C
C        Light time corrections
C        ======================
C
C           Observer-target surface point vector
C           ------------------------------------
C
C           Let ET be the epoch at which an observation or remote
C           sensing measurement is made, and let ET - LT ("LT" stands
C           for "light time") be the epoch at which the photons
C           received at ET were emitted from the surface point SPOINT.
C           Note that the light time between the surface point and
C           observer will generally differ from the light time between
C           the target body's center and the observer.
C
C
C           Target body's orientation
C           -------------------------
C
C           Using the definitions of ET and LT above, the target body's
C           orientation at ET - LT is used. The surface normal is
C           dependent on the target body's orientation, so the body's
C           orientation model must be evaluated for the correct epoch.
C
C
C           Target body -- Sun vector
C           -------------------------
C
C           The surface features on the target body near SPOINT will
C           appear in a measurement made at ET as they were at ET-LT.
C           In particular, lighting on the target body is dependent on
C           the apparent location of the Sun as seen from the target
C           body at ET-LT. So, a second light time correction is used
C           to compute the position of the Sun relative to the surface
C           point.
C
C
C        Stellar aberration corrections
C        ==============================
C
C        Stellar aberration corrections are applied only if
C        light time corrections are applied as well.
C
C           Observer-target surface point body vector
C           -----------------------------------------
C
C           When stellar aberration correction is performed, the
C           direction vector SRFVEC is adjusted so as to point to the
C           apparent position of SPOINT: considering SPOINT to be an
C           ephemeris object, SRFVEC points from the observer's
C           position at ET to the light time and stellar aberration
C           corrected position of SPOINT.
C
C           Target body-Sun vector
C           ----------------------
C
C           The target body-Sun vector is the apparent position of the
C           Sun, corrected for light time and stellar aberration, as
C           seen from the target body at time ET-LT.  
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
C              'DSK/UNPRIORITIZED/SURFACES = '
C           // '"Mars MEGDR 64 PIXEL/DEG", 3'
C
C
C        Aberration corrections using DSK data
C        -------------------------------------
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
C     1) Find the phase, solar incidence, and emission angles at the
C        sub-solar and sub-spacecraft points on Mars as seen from the
C        Mars Global Surveyor spacecraft at a specified UTC time. Use
C        light time and stellar aberration corrections.
C
C        Use both an ellipsoidal Mars shape model and topographic data
C        provided by a DSK file.
C
C        Use the meta-kernel shown below to load the required SPICE
C        kernels.
C
C           KPL/MK
C
C           File: ilumin_ex1.tm
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
C                                  'mgs_ext12_ipng_mgs95j.bsp',
C                                  'megr90n000cb_plate.bds'      )
C           \begintext
C
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
C           DOUBLE PRECISION      DPR
C     C
C     C     Local parameters
C     C
C           CHARACTER*(*)         F1
C           PARAMETER           ( F1     = '(A,F15.9)' )
C
C           CHARACTER*(*)         F2
C           PARAMETER           ( F2     = '(A)' )
C
C           CHARACTER*(*)         F3
C           PARAMETER           ( F3     = '(A,2(2X,L))' )
C
C           CHARACTER*(*)         META
C           PARAMETER           ( META   = 'ilumin_ex1.tm' )
C
C           INTEGER               NAMLEN
C           PARAMETER           ( NAMLEN = 32 )
C
C           INTEGER               TIMLEN
C           PARAMETER           ( TIMLEN = 25 )
C
C           INTEGER               CORLEN
C           PARAMETER           ( CORLEN = 5 )
C
C           INTEGER               MTHLEN
C           PARAMETER           ( MTHLEN = 50 )
C
C           INTEGER               NMETH
C           PARAMETER           ( NMETH  = 2 )
C     C
C     C     Local variables
C     C
C           CHARACTER*(CORLEN)    ABCORR
C           CHARACTER*(NAMLEN)    FIXREF
C           CHARACTER*(MTHLEN)    ILUMTH ( NMETH )
C           CHARACTER*(NAMLEN)    OBSRVR
C           CHARACTER*(MTHLEN)    SUBMTH ( NMETH )
C           CHARACTER*(NAMLEN)    TARGET
C           CHARACTER*(TIMLEN)    UTC
C
C           DOUBLE PRECISION      ET
C           DOUBLE PRECISION      SRFVEC ( 3 )
C           DOUBLE PRECISION      SSCEMI
C           DOUBLE PRECISION      SSCPHS
C           DOUBLE PRECISION      SSCPT  ( 3 )
C           DOUBLE PRECISION      SSCSOL
C           DOUBLE PRECISION      SSLEMI
C           DOUBLE PRECISION      SSLPHS
C           DOUBLE PRECISION      SSLSOL
C           DOUBLE PRECISION      SSOLPT ( 3 )
C           DOUBLE PRECISION      TRGEPC
C
C           INTEGER               I
C
C
C     C
C     C     Initial values
C     C
C           DATA                  ILUMTH / 'Ellipsoid',
C          .                               'DSK/Unprioritized' /
C
C           DATA                  SUBMTH / 'Near Point/Ellipsoid',
C          .                            'DSK/Nadir/Unprioritized' /
C
C     C
C     C     Load kernel files.
C     C
C           CALL FURNSH ( META )
C     C
C     C     Convert the UTC request time string to seconds past
C     C     J2000 TDB.
C     C
C           UTC = '2003 OCT 13 06:00:00 UTC'
C
C           CALL UTC2ET ( UTC, ET )
C
C           WRITE (*,F2) ' '
C           WRITE (*,F2) 'UTC epoch is '//UTC
C     C
C     C     Assign observer and target names. The acronym MGS
C     C     indicates Mars Global Surveyor. See NAIF_IDS for a
C     C     list of names recognized by SPICE. Also set the
C     C     aberration correction flag.
C     C
C           TARGET = 'Mars'
C           OBSRVR = 'MGS'
C           FIXREF = 'IAU_MARS'
C           ABCORR = 'CN+S'
C
C           DO I = 1, NMETH
C     C
C     C        Find the sub-solar point on Mars as
C     C        seen from the MGS spacecraft at ET. Use the
C     C        "near point" style of sub-point definition
C     C        when the shape model is an ellipsoid, and use
C     C        the "nadir" style when the shape model is
C     C        provided by DSK data. This makes it easy to 
C     C        verify the solar incidence angle when
C     C        the target is modeled as an  ellipsoid.
C     C
C              CALL SUBSLR ( SUBMTH(I),  TARGET,  ET,
C          .                 FIXREF,     ABCORR,  OBSRVR,
C          .                 SSOLPT,     TRGEPC,  SRFVEC  )
C     C
C     C        Now find the sub-spacecraft point.
C     C
C              CALL SUBPNT ( SUBMTH(I),  TARGET,  ET,
C          .                 FIXREF,     ABCORR,  OBSRVR,
C          .                 SSCPT,      TRGEPC,  SRFVEC )
C     C
C     C        Find the phase, solar incidence, and emission
C     C        angles at the sub-solar point on Mars as
C     C        seen from MGS at time ET.
C     C
C              CALL ILUMIN ( ILUMTH(I), TARGET, 
C          .                 ET,        FIXREF,  ABCORR,
C          .                 OBSRVR,    SSOLPT,  TRGEPC,
C          .                 SRFVEC,    SSLPHS,  SSLSOL,
C          .                 SSLEMI                      )
C     C
C     C        Do the same for the sub-spacecraft point.
C     C
C              CALL ILUMIN ( ILUMTH(I), TARGET,  
C          .                 ET,        FIXREF,  ABCORR,
C          .                 OBSRVR,    SSCPT,   TRGEPC,
C          .                 SRFVEC,    SSCPHS,  SSCSOL,
C          .                 SSCEMI                      )
C     C
C     C        Convert the angles to degrees and write them out.
C     C
C              SSLPHS = DPR() * SSLPHS
C              SSLSOL = DPR() * SSLSOL
C              SSLEMI = DPR() * SSLEMI
C
C              SSCPHS = DPR() * SSCPHS
C              SSCSOL = DPR() * SSCSOL
C              SSCEMI = DPR() * SSCEMI
C
C              WRITE (*,F2) ' '
C              WRITE (*,F2) '   ILUMIN method: '//ILUMTH(I)
C              WRITE (*,F2) '   SUBPNT method: '//SUBMTH(I)
C              WRITE (*,F2) '   SUBSLR method: '//SUBMTH(I)
C              WRITE (*,F2) ' '
C              WRITE (*,F2) '      Illumination angles at the '
C          .   //           'sub-solar point:'
C              WRITE (*,F2) ' '
C
C              WRITE (*,F1) '      Phase angle           (deg.): ',
C          .                SSLPHS
C              WRITE (*,F1) '      Solar incidence angle (deg.): ',
C          .                SSLSOL
C              WRITE (*,F1) '      Emission angle        (deg.): ',
C          .                SSLEMI
C              WRITE (*,F2) ' '
C
C              IF ( I .EQ. 1 ) THEN
C                 WRITE (*,F2) '        The solar incidence angle '
C          .      //           'should be 0.'
C                 WRITE (*,F2) '        The emission and phase '
C          .      //           'angles should be equal.'
C                 WRITE (*,F2) ' '
C              END IF
C
C
C              WRITE (*,F2) '      Illumination angles at the '
C          .   //          'sub-s/c point:'
C              WRITE (*,F2) ' '
C              WRITE (*,F1) '      Phase angle           (deg.): ',
C          .               SSCPHS
C              WRITE (*,F1) '      Solar incidence angle (deg.): ',
C          .               SSCSOL
C              WRITE (*,F1) '      Emission angle        (deg.): ',
C          .               SSCEMI
C              WRITE (*,F2) ' '
C
C              IF ( I .EQ. 1 ) THEN
C                 WRITE (*,F2) '        The emission angle '
C          .      //           'should be 0.'
C                 WRITE (*,F2) '        The solar incidence '
C          .      //           'and phase angles should be equal.'
C              END IF
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
C        UTC epoch is 2003 OCT 13 06:00:00 UTC
C
C           ILUMIN method: Ellipsoid
C           SUBPNT method: Near Point/Ellipsoid
C           SUBSLR method: Near Point/Ellipsoid
C
C              Illumination angles at the sub-solar point:
C
C              Phase angle           (deg.):   138.370270685
C              Solar incidence angle (deg.):     0.000000000
C              Emission angle        (deg.):   138.370270685
C
C                The solar incidence angle should be 0.
C                The emission and phase angles should be equal.
C
C              Illumination angles at the sub-s/c point:
C
C              Phase angle           (deg.):   101.439331040
C              Solar incidence angle (deg.):   101.439331041
C              Emission angle        (deg.):     0.000000002
C
C                The emission angle should be 0.
C                The solar incidence and phase angles should be equal.
C
C           ILUMIN method: DSK/Unprioritized
C           SUBPNT method: DSK/Nadir/Unprioritized
C           SUBSLR method: DSK/Nadir/Unprioritized
C
C              Illumination angles at the sub-solar point:
C
C              Phase angle           (deg.):   138.387071677
C              Solar incidence angle (deg.):     0.967122745
C              Emission angle        (deg.):   137.621480599
C
C              Illumination angles at the sub-s/c point:
C
C              Phase angle           (deg.):   101.439331359
C              Solar incidence angle (deg.):   101.555993667
C              Emission angle        (deg.):     0.117861156 
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
C     N.J. Bachman   (JPL)
C     S.C. Krening   (JPL)
C
C$ Version
C
C-    SPICELIB Version 2.0.0, 04-APR-2017 (NJB)
C
C        Fixed some header comment typos.
C
C     15-AUG-2016 (NJB) 
C
C        Now supports DSK usage. No longer includes zzabcorr.inc.
C        String 'SUN' passed to ILLUMG has been changed to '10'.
C
C        Now supports transmission aberration corrections.
C
C     
C-    SPICELIB Version 1.2.0, 04-APR-2011 (NJB) (SCK)
C
C        The routine has been completely re-implemented:
C        it now calls ILLUMG.
C
C        The meta-kernel used for the header example program
C        has been updated. The example program outputs have
C        been updated as well.
C
C        References to the new PXFRM2 routine were added
C        to the Detailed Output section.  
C
C-    SPICELIB Version 1.1.0, 17-MAY-2010 (NJB) 
C
C        Bug fix: ILUMIN now returns immediately if a target
C        radius lookup fails.
C
C-    SPICELIB Version 1.0.1, 06-FEB-2009 (NJB) 
C
C        Typo correction: changed FIXFRM to FIXREF in header
C        documentation. Meta-kernel name suffix was changed to
C        ".tm" in header code example.
C
C-    SPICELIB Version 1.0.0, 02-MAR-2008 (NJB)
C
C-&
 
C$ Index_Entries
C
C     illumination angles
C     lighting angles
C     phase angle
C     solar incidence angle
C     emission angle
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
      LOGICAL               RETURN
 
C
C     Standard SPICE error handling.
C
      IF ( RETURN () ) THEN
         RETURN
      END IF

      CALL CHKIN ( 'ILUMIN' )

      CALL ILLUMG ( METHOD,  TARGET,  '10',    ET,     
     .              FIXREF,  ABCORR,  OBSRVR,  SPOINT,
     .              TRGEPC,  SRFVEC,  PHASE,   SOLAR,   EMISSN )
  
 
      CALL CHKOUT ( 'ILUMIN' )
      RETURN
      END
